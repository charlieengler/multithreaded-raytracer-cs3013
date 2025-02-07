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

#define CHECK(x)	do { if (!(x)) { fprintf(stderr, "%s:%d CHECK failed: %s, errno %d %s\n", __FILE__, __LINE__, #x, errno, strerror(errno)); abort(); } } while(0)

void render_scene(struct framebuffer_pt4 *fb, const struct context* ctx);

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
	struct framebuffer_pt4 *fb = NULL;
	int render_to_console = 1;
	if (argc > 2) {
		render_to_console = 0;
		///////////////////////////////////////////////////////////////////////////////////////
		// HINT: changing the resolutions here will alter the performance. If you want bmps  //
		// but faster, try lowering the resolution here.                                     //
		///////////////////////////////////////////////////////////////////////////////////////
		//fb = new_framebuffer_pt4(2560, 1440);
		fb = new_framebuffer_pt4(2048, 1536);
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

		fb = new_framebuffer_pt4(ws.ws_col, ws.ws_row);
	}

	///////////////////////////////////////////////////////////////////////////////////////
	// Now we have a framebuffer and a scene graph.                                      //
	// Alternate render and physics passes.                                              //
	// However: we can parallelize the output here, as long as we are not corrupting the //
	// framebuffer whilst outputting.                                                    //
	// TODO: section 2: instead of one framebuffer, use two, and write file output for   //
	//                  frame N-1 whilst drawing frame N.                                //
	///////////////////////////////////////////////////////////////////////////////////////

	for (int frame = 0; frame < 100; frame++) {
		render_scene(fb, ctx);
		step_physics(ctx);
		if (render_to_console) {
			render_console(fb);
		} else {
			char filepath[128];
			snprintf(filepath, sizeof(filepath)-1, "%s-%05d.bmp", argv[2], frame);
			CHECK(render_bmp(fb, filepath) == 0);
		}
	}

out:
	yylex_destroy(scanner);
	if (finput) fclose(finput);
	free_context(ctx);
	if (fb) free_framebuffer_pt4(fb);

	return 0;
}


void render_scene(struct framebuffer_pt4 *fb, const struct context *ctx) {
	double left_right_angle;
	double up_down_angle;
	int xmax = fb->width;
	int ymax = fb->height;
	if (xmax > ymax) {
		left_right_angle = M_PI / 3;
		up_down_angle = left_right_angle / xmax * ymax;
	} else {
		up_down_angle = M_PI / 3;
		left_right_angle = up_down_angle / ymax * xmax;
	}

	double left_right_start = - left_right_angle / 2;
	double left_right_step = left_right_angle / (xmax - 1);
	double up_down_start = up_down_angle / 2;
	double up_down_step = up_down_angle / (ymax - 1);

	//printf("left_right_start %lf left_right_step %lf up_down_start %lf up_down_step %lf\n", left_right_start, left_right_step, up_down_start, up_down_step);

	for (int x = 0; x < xmax; x++) {
		double xangle = -(left_right_start + left_right_step * x);
		for (int y = 0; y < ymax; y++) {
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

}

