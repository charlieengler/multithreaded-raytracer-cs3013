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

#define NUM_FRAME_BUFFERS 2

// TODO: Comments
#define NUM_RENDER_THREADS 20

// TODO: Comments
typedef struct render_args {
	struct framebuffer_pt4 *fb;
	struct context *ctx;
	int index;
	int xmin, xmax;
} render_args;

// TODO: Free these monstrosities
sem_t ***fb_ready_mutexes;
sem_t ***fb_done_mutexes;

#define CHECK(x)	do { if (!(x)) { fprintf(stderr, "%s:%d CHECK failed: %s, errno %d %s\n", __FILE__, __LINE__, #x, errno, strerror(errno)); abort(); } } while(0)

// TODO: Comments
void *render_scene(void *args);

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
		const int x_res = 267;
		const int y_res = 200;
		// const int x_res = 2048;
		// const int y_res = 1536;

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

	// TODO: Comments
	pthread_t render_threads[NUM_RENDER_THREADS];

	fb_ready_mutexes = (sem_t***)malloc(sizeof(sem_t**) * NUM_FRAME_BUFFERS);
	fb_done_mutexes = (sem_t***)malloc(sizeof(sem_t**) * NUM_FRAME_BUFFERS);
	for(int i = 0; i < NUM_FRAME_BUFFERS; i++) {
		sem_t **ready_mutexes = (sem_t**)malloc(sizeof(sem_t*) * NUM_RENDER_THREADS);
		fb_ready_mutexes[i] = ready_mutexes;

		sem_t **done_mutexes = (sem_t**)malloc(sizeof(sem_t*) * NUM_RENDER_THREADS);
		fb_done_mutexes[i] = done_mutexes;

		for(int j = 0; j < NUM_RENDER_THREADS; j++) {
			sem_t *ready_mutex = malloc(sizeof(sem_t));
			sem_t *done_mutex = malloc(sizeof(sem_t));

			fb_ready_mutexes[i][j] = ready_mutex;
			fb_done_mutexes[i][j] = done_mutex;

			// TODO: Error checking
			sem_init(ready_mutex, 0, 0);
			sem_init(done_mutex, 0, 0);
		}
	}

	const int x_pixels_per_thread = fbs[0]->width / NUM_RENDER_THREADS;

	for(int rt = 0; rt < NUM_RENDER_THREADS; rt++) {
		const int xmin = rt * x_pixels_per_thread;
		const int xmax = rt + 1 == NUM_RENDER_THREADS ? fbs[0]->width : rt * x_pixels_per_thread + x_pixels_per_thread;

		render_args *r_args = (render_args*)malloc(sizeof(render_args));
		r_args->fb = fbs[0];
		r_args->ctx = ctx;
		r_args->index = rt;
		r_args->xmin = xmin;
		r_args->xmax = xmax;

		int c_err;
		if((c_err = pthread_create(&render_threads[rt], NULL, &render_scene, (void*)r_args))) {
			printf("[ray.c -> main] Unable to create render thread: %d\n", c_err);
			goto out;
		}
	}

	// TODO: Comments for all of this
	for (int frame = 0; frame < 100; frame++) {
		for(int rt = 0; rt < NUM_RENDER_THREADS; rt++) {
			// TODO: Error checking
			sem_post(fb_ready_mutexes[0][rt]);
		}

		for(int rt = 0; rt < NUM_RENDER_THREADS; rt++) {
			// TODO: Error checking
			sem_wait(fb_done_mutexes[0][rt]);
		}

		calc_velocities(ctx);
		update_positions(ctx);

		if (render_to_console) {
			render_console(fbs[0]);
		} else {
			char filepath[128];
			snprintf(filepath, sizeof(filepath)-1, "%s-%05d.bmp", argv[2], frame);
			CHECK(render_bmp(fbs[0], filepath) == 0);
		}
	}

	for(int rt = 0; rt < NUM_RENDER_THREADS; rt++) {
		int c_err;
		if((c_err = pthread_cancel(render_threads[rt]))) {
			printf("[ray.c -> main] Unable to cancel render thread: %d\n", c_err);
			goto out;
		}
	}

out:
	yylex_destroy(scanner);
	if (finput) fclose(finput);
	free_context(ctx);
	// TODO: Comments
	for(int i = 0; i < NUM_FRAME_BUFFERS; i++) {
		if (fbs[i]) free_framebuffer_pt4(fbs[i]);
	}
	if (fbs) free(fbs);

	return 0;
}


void *render_scene(void *args) {
	render_args *r_args = (render_args*)args;

	struct framebuffer_pt4 *fb = r_args->fb;
	struct context *ctx = r_args->ctx;

	const int xmax = fb->width;
	const int ymax = fb->height;

	const int pix_x_min = r_args->xmin;
	const int pix_x_max = r_args->xmax;
	const int pix_y_min = 0;
	const int pix_y_max = ymax;

	const int index = r_args->index;

	free(r_args);

	while(1) {
		// TODO: Error checking
		sem_wait(fb_ready_mutexes[0][index]);

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
				framebuffer_pt4_set(fb, x, y, px_color);
			}
		}

		// TODO: Error checking
		sem_post(fb_done_mutexes[0][index]);
	}
	
	return NULL;
}