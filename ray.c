/*
Outputs from time command:

	Baseline (100 frames, 2048x1536):
		real    10m6.024s
		user    2m37.799s
		sys     0m10.268s

	Multithreaded Rendering (100 frames, 2048x1536):
		real    7m22.293s
		user    3m0.609s
		sys     0m10.756s
*/

#include <stdio.h>
#include <sys/ioctl.h>
#include <errno.h>
#include <semaphore.h>
#include <pthread.h>

#include "ray_ast.h"
#include "ray_math.h"
#include "ray.yacc.generated_h"
#include "ray.lex.generated_h"

#include "ray_render.h"
#include "ray_bmp.h"
#include "ray_physics.h"
#include "ray_console.h"

// This should always be 2, I just didn't want to have any out of place constants
#define NUM_FRAME_BUFFERS 2

// Defines the number of threads used for rendering the pixels in the scene itself
#define NUM_RENDER_THREADS 5

// Render arguments passed to the functions instantiated in each thread in the thread pool
typedef struct render_args {
	// Scene context
	struct context *ctx;
	// Index of the render thread in the thread pool
	int rt;
	// The minimum pixel on the x axis to start rendering at
	int xmin;
	// The maximum pixel on the x axis to finish rendering at
	int xmax;
} render_args;

// Mutexes denoting which threads have data that is ready to process
sem_t *render_ready_mutexes;
// Mutexes denoting which threads in the rendering thread pool are done with their data
sem_t *render_done_mutexes;

// Mutex for when the frame buffer is being swapped
sem_t buf_swap_mutex;

// Global definition of the frame buffer being actively rendered to
struct framebuffer_pt4 *render_fb;
// Global definition of the frame buffer being saved to the output file
struct framebuffer_pt4 *save_fb;

// Global variable for the filepath of the current save buffer
char current_filepath[128];

#define CHECK(x)	do { if (!(x)) { fprintf(stderr, "%s:%d CHECK failed: %s, errno %d %s\n", __FILE__, __LINE__, #x, errno, strerror(errno)); abort(); } } while(0)

// Function definition for the render scene function saving to the render frame buffer
void *render_scene(void *args);
// Function definition for the save scene function saving to a file from the save frame buffer
void *save_scene(void *args);

int main(int argc, char **argv) {
	int rc;

	struct context *ctx = new_context();
	FILE *finput = NULL;
	yyscan_t scanner;
	yylex_init(&scanner);

	// Read from a script. By default this is stdin.
	if (argc > 1) {
		// If a file is specified as a command line argument, read from that instead of stdin.
		const char *source = argv[1];
		finput = fopen(source, "rb");
		if (finput == NULL) {
			fprintf(stderr, "Could not open '%s' for reading, errno %d (%s)\n", source, errno, strerror(errno));
			return 1;
		}
		yyset_in(finput, scanner);
	}
	// Parse the input file and run the parsed script if parsing was successful.
	if ((rc = yyparse(ctx, scanner)) != 0) {
		fprintf(stderr, "Parse failure for script\n");
		goto out;		
	}

	// Calculate framebuffer size. If we're outputting into a bmp file, use a high resolution.
	// If we're rendering to the active console, use ioctls to find the window size.
	// Initialize a framebuffer with the chosen resolution.
	//
	// For section 2.2, we will need to replace the single framebuffer here with an array containing
	// two framebuffers, one for even frames, one for odd.
	struct framebuffer_pt4 **fbs = malloc(sizeof(struct framebuffer_pt4*) * NUM_FRAME_BUFFERS);
	int render_to_console = 1;
	if (argc > 2) {
		render_to_console = 0;
		///////////////////////////////////////////////////////////////////////////////////////
		// HINT: changing the resolutions here will alter the performance. If you want bmps  //
		// but faster, try lowering the resolution here.                                     //
		///////////////////////////////////////////////////////////////////////////////////////
		// const int x_res = 267;
		// const int y_res = 200;
		// Using the default resolutions
		const int x_res = 2048;
		const int y_res = 1536;

		// Initialize the frame buffers for file writing (2 of them)
		for(int i = 0; i < NUM_FRAME_BUFFERS; i++)
			fbs[i] = new_framebuffer_pt4(x_res, y_res);
	} else {
		struct winsize ws;

		if (isatty(STDIN_FILENO)) {
			if ((rc = ioctl(0, TIOCGWINSZ, &ws)) < 0) {
				fprintf(stderr, "Failed to get window size: %d %s\n", errno, strerror(errno));
				return 1;
			}
			printf ("cols (x) %d lines (y) %d\n", ws.ws_col, ws.ws_row);
		} else {
			ws.ws_row = 128;
			ws.ws_col = 128;
		}

		for(int i = 0; i < NUM_FRAME_BUFFERS; i++)
			fbs[i] = new_framebuffer_pt4(ws.ws_col, ws.ws_row);
	}

	///////////////////////////////////////////////////////////////////////////////////////
	// Now we have a framebuffer and a scene graph.                                      //
	// Alternate render and physics passes.                                              //
	// However: we can parallelize the output here, as long as we are not corrupting the //
	// framebuffer whilst outputting.                                                    //
	// TODO: section 2: instead of one framebuffer, use two, and write file output for   //
	//                  frame N-1 whilst drawing frame N.                                //
	///////////////////////////////////////////////////////////////////////////////////////

	// TODO: Free these
	// Thread pool for setting up the rendering threads
	pthread_t *render_threads = malloc(sizeof(pthread_t) * NUM_RENDER_THREADS);
	// A single thread for saving the file, this is not part of a pool
	pthread_t file_thread;

	// This determines the number of pixels per render thread on the x axis based on the horizontal
	// resolution and the number of render threads in the pool
	const int x_pixels_per_thread = fbs[0]->width / NUM_RENDER_THREADS;

	// The frame buffer used for rendering is initialized as the second one so the program has
	// something to save during the first loop
	render_fb = fbs[1];

	// Initialize the buffer swap mutex to zero because they should only be swapped after rendering is done
	sem_init(&buf_swap_mutex, 0, 0);

	// Initialize the render mutex arrays
	render_ready_mutexes = (sem_t*)malloc(sizeof(sem_t) * NUM_RENDER_THREADS);
	render_done_mutexes = (sem_t*)malloc(sizeof(sem_t) * NUM_RENDER_THREADS);
	for(int rt = 0; rt < NUM_RENDER_THREADS; rt++) {
		// Initialize the mutex for render data being ready on a render thread in the pool and check for an error
		if(sem_init(&(render_ready_mutexes[rt]), 0, 0) != 0) {
			printf("[ray.c -> main] Unable to initialize render ready mutex %d: %d\n", rt, errno);
			goto out;
		}

		// Initialize the mutex for render data being processed on a render thread in the pool and check for an error
		if(sem_init(&(render_done_mutexes[rt]), 0, 0) != 0) {
			printf("[ray.c -> main] Unable to initialize render done mutex %d: %d\n", rt, errno);
			goto out;
		}

		// The minimum x pixel value is the index of the render thread times the number of x pixels per thread
		const int xmin = rt * x_pixels_per_thread;
		// The maximum x pixel value is the index of the render thread plus one times the number of x pixels per thread
		// If the next index is the last index, then the maximum x pixel should be the x resolution of the window
		const int xmax = rt + 1 == NUM_RENDER_THREADS ? fbs[0]->width : rt * x_pixels_per_thread + x_pixels_per_thread;

		// Create a new render argument structure for a particular render thread. This is constant to the thread in the pool
		// so the thread renders the same pixels every time.
		render_args *r_args = (render_args*)malloc(sizeof(render_args));
		r_args->ctx = ctx;
		r_args->rt = rt;
		r_args->xmin = xmin;
		r_args->xmax = xmax;

		// Create a thread in the pool and give it the render_scene function to call with the relevant arguments
		int c_err;
		if((c_err = pthread_create(&(render_threads[rt]), NULL, &render_scene, (void*)r_args))) {
			printf("[ray.c -> main] Unable to create render thread: %d\n", c_err);
			goto out;
		}
	}

	// These for loops start the rendering process in the second frame buffer so there is data to save on the first loop
	// Posts to each render ready mutex in the thread pool to inform the render threads that the scene positions are modified and ready to render
	for(int rt = 0; rt < NUM_RENDER_THREADS; rt++) {
		if(sem_post(&(render_ready_mutexes[rt])) != 0) {
			printf("[ray.c -> main] Sem post error (setup frame) on render ready mutex %d: %d\n", rt, errno);
			goto out;
		}
	}

	// Velocities can be calculated while rendering because it doesn't overwrite positions
	calc_velocities(ctx);
	
	// Waits on render data to be done in each of the threads in the pool before updating positions
	for(int rt = 0; rt < NUM_RENDER_THREADS; rt++) {
		if(sem_wait(&(render_done_mutexes[rt])) != 0) {
			printf("[ray.c -> main] Sem wait error (setup frame) on render done mutex %d: %d\n", rt, errno);
			goto out;
		}
	}

	// Velocities need to be updated after the rendering is done because it modifies positions and frame would be torn if these
	// were updated in parallel
	update_positions(ctx);

	// Loops through all frames in the animation
	for (int frame = 0; frame < 100; frame++) {
		// Swap the save frame buffer for the render frame buffer
		save_fb = render_fb;

		if(frame % 2 == 0) {
			// On an even frame, the render buffer is the first frame buffer
			render_fb = fbs[0];

			// Indicate that the buffers are swapped
			if(sem_post(&buf_swap_mutex) != 0)
				printf("[ray.c -> main] Even frame sem post error on buffer swap mutex: %d\n", errno);
		} else {
			// On an odd frame, the render buffer is the second frame buffer
			render_fb = fbs[1];

			// Indicate that the buffers are swapped
			if(sem_post(&buf_swap_mutex) != 0)
				printf("[ray.c -> main] Odd frame sem post error on buffer swap mutex: %d\n", errno);
		}

		// Wait for the frame buffers to swap before doing anything else because this would mess things up
		if(sem_wait(&buf_swap_mutex) != 0)
			printf("[ray.c -> main] Sem wait error on buffer swap mutex: %d\n", errno);

		if (render_to_console) {
			render_console(fbs[frame % 2]);
		} else {
			// Write the file name to the global file path so it can be used by the save function
			snprintf(current_filepath, sizeof(current_filepath)-1, "%s-%05d.bmp", argv[2], frame);
			
			// The file thread is spawned each frame because it's only one thread and we can wait for it to
			// finish execution with join later on (makes my life so much easier lol). It calls the save_scene
			// function with no arguments
			int c_err;
			if((c_err = pthread_create(&file_thread, NULL, &save_scene, NULL))) {
				printf("[ray.c -> main] Unable to create file thread: %d\n", c_err);
				goto out;
			}
		}

		// Posts to each render ready mutex in the thread pool to inform the render threads that the scene positions are modified and ready to render
		for(int rt = 0; rt < NUM_RENDER_THREADS; rt++) {
			if(sem_post(&(render_ready_mutexes[rt])) != 0) {
				printf("[ray.c -> main] Sem post error (setup frame) on render ready mutex %d: %d\n", rt, errno);
				goto out;
			}
		}

		// Velocities can be calculated while rendering because it doesn't overwrite positions
		calc_velocities(ctx);

		// Waits on render data to be done in each of the threads in the pool before updating positions
		for(int rt = 0; rt < NUM_RENDER_THREADS; rt++) {
			if(sem_wait(&(render_done_mutexes[rt])) != 0) {
				printf("[ray.c -> main] Sem wait error (setup frame) on render done mutex %d: %d\n", rt, errno);
				goto out;
			}
		}

		// Velocities need to be updated after the rendering is done because it modifies positions and frame would be torn if these
		// were updated in parallel
		update_positions(ctx);

		// Wait for the file thread to finish execution and then join it onto the main thread before continuing the loop
		int j_err;
		if((j_err = pthread_join(file_thread, NULL))) {
			printf("[ray.c -> main] Unable to join file thread: %d\n", j_err);
		}
	}

	// Since the render threads are constantly executing, they need to be canceled after their work is done
	int c_err;
	for(int rt = 0; rt < NUM_RENDER_THREADS; rt++) {
		if((c_err = pthread_cancel(render_threads[rt]))) {
			printf("[ray.c -> main] Unable to cancel render thread: %d\n", c_err);
			goto out;
		}
	}

out:
	yylex_destroy(scanner);
	if (finput) fclose(finput);
	free_context(ctx);
	// Free both frame buffers
	for(int i = 0; i < NUM_FRAME_BUFFERS; i++) {
		if (fbs[i]) free_framebuffer_pt4(fbs[i]);
	}
	if (fbs) free(fbs);
	// There will be NO memory leaks here :)
	if (render_ready_mutexes) free(render_ready_mutexes);
	if (render_done_mutexes) free(render_done_mutexes);
	if (render_threads) free(render_threads);

	return 0;
}

// Render scene file set up for thread execution with void pointer arguments
void *render_scene(void *args) {
	// Void pointer arguments are cast to a render arguments struct
	render_args *r_args = (render_args*)args;
	
	// Extract the scene context from the render arguments
	struct context *ctx = r_args->ctx;

	// These are used for the raytracing calculations, not for rendering pixels
	const int xmax = render_fb->width;
	const int ymax = render_fb->height;

	// These are the pixel bounds for the slice rendered by this render thread
	// Variable minimum x pixel rendered on this thread
	const int pix_x_min = r_args->xmin;
	// Variable maximum x pixel rendered on this thread
	const int pix_x_max = r_args->xmax;
	// Y min and max are always constant
	const int pix_y_min = 0;
	const int pix_y_max = ymax;

	// The index of this render thread in the thread pool
	const int rt = r_args->rt;
	
	// Argument structure is freed because all the variables are saved
	free(r_args);

	// Permanent loop that runs work when it's available
	while(1) {
		// Wait for render work to be ready on this thread
		if(sem_wait(&(render_ready_mutexes[rt])) != 0)
			printf("[ray.c -> render_scene] Sem wait error on render ready mutex %d: %d\n", rt, errno);

		double left_right_angle;
		double up_down_angle;
		if (xmax > ymax) {
			left_right_angle = M_PI / 3;
			up_down_angle = left_right_angle / xmax * ymax;
		} else {
			up_down_angle = M_PI / 3;
			left_right_angle = up_down_angle / ymax * xmax;
		}

		double left_right_start = -left_right_angle / 2;
		double left_right_step = left_right_angle / (xmax - 1);
		double up_down_start = up_down_angle / 2;
		double up_down_step = up_down_angle / (ymax - 1);

		//printf("left_right_start %lf left_right_step %lf up_down_start %lf up_down_step %lf\n", left_right_start, left_right_step, up_down_start, up_down_step);

		// The x pixel bounds are the minimum and maximum assigned to this thread so it only renders a
		// smaller section for performance and parallelization
		for (int x = pix_x_min; x < pix_x_max; x++) {
			double xangle = -(left_right_start + left_right_step * x);
			for (int y = pix_y_min; y < pix_y_max; y++) {
				double yangle = up_down_start - up_down_step * y;

				// I'm 99% sure this is wrong but it looks okay
				pt3 direction = {{sin(xangle), sin(yangle), cos(yangle) * cos(xangle)}};
				pt3_normalize_mut(&direction);
				ray r = {{{0, 0, -20}}, direction};
				//printf("ray: %lf %lf %lf\n", direction.v[0], direction.v[1], direction.v[2]);
				pt4 px_color = {0};
				raytrace(ctx, &r, &px_color, 3);
				// Rendered to the global render frame buffer
				framebuffer_pt4_set(render_fb, x, y, px_color);
			}
		}

		// Broadcast that render work is done on this thread
		if(sem_post(&(render_done_mutexes[rt])) != 0)
			printf("[ray.c -> render_scene] Sem post error on render ready mutex %d: %d\n", rt, errno);
	}
	
	return NULL;
}

// Save scene file set up for thread execution with void pointer arguments
void *save_scene(void *args) {
	// The filepath is read from the global variable
	char *filepath = malloc(sizeof(current_filepath));
	strcpy(filepath, current_filepath);

	// The save file buffer is written into the bitmap file for this frame
	CHECK(render_bmp(save_fb, filepath) == 0);

	// The thread-specific filepath is free
	free(filepath);

	return NULL;
}