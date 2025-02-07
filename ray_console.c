#include <stdio.h>

#include "ray_render.h"

void render_console(struct framebuffer_pt4 *fb) {
	for (int y = 0; y < fb->height; y++) {
		for (int x = 0; x < fb->width; x++) {
			const pt4 *p = framebuffer_pt4_get(fb, x, y);
			//printf("x %d y %d is %lf %lf %lf\n", x, y, p->v[0], p->v[1], p->v[2]);
			printf("\033[48;2;%u;%u;%um ",
					color_double_to_u8(p->v[0]),
					color_double_to_u8(p->v[1]),
					color_double_to_u8(p->v[2]));
		}
		printf("\n");
	}
}
