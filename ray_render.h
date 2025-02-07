#ifndef RAY_RENDER_H__
#define RAY_RENDER_H__

#include <inttypes.h>
#include <stdlib.h>
#include <string.h>

#include "ray_ast.h"

#define DEFINE_FRAMEBUFFER(t, n)								\
struct framebuffer_##t {									\
	int width;										\
	int height;										\
	int channels;										\
	t *pixels;										\
};												\
												\
static inline struct framebuffer_##t *new_framebuffer_##t(int width, int height) {		\
	struct framebuffer_##t *ret = malloc(sizeof(*ret) + sizeof(t) * n * width * height);	\
	ret->width = width;									\
	ret->height = height;									\
	ret->channels = n;									\
	ret->pixels = (void*)(ret + 1);								\
	return ret;										\
}												\
												\
static inline void free_framebuffer_##t(struct framebuffer_##t *fb) {				\
	free(fb);										\
}												\
												\
static inline int framebuffer_##t##_pixel_index(struct framebuffer_##t *fb, int x, int y) {	\
	return n * (y * fb->width + x);								\
}												\
												\
static inline void framebuffer_##t##_set(struct framebuffer_##t *fb, int x, int y, t px) {	\
	fb->pixels[framebuffer_##t##_pixel_index(fb, x, y)] = px;				\
}												\
												\
static inline t* framebuffer_##t##_get(struct framebuffer_##t *fb, int x, int y) {		\
	return &fb->pixels[framebuffer_##t##_pixel_index(fb, x, y)];				\
}

//DEFINE_FRAMEBUFFER(uint8_t)	// for the rgb data before export
DEFINE_FRAMEBUFFER(pt4, 1)	// for the actual rendering passes

static inline uint8_t color_double_to_u8(double d) {
	if (d < 0) return 0;
	if (d >= 1) return 255;
	return (uint8_t)(d * 255);
}

int raytrace(const struct context *ctx, const ray *r, pt4 *ret, int depth);

#endif	// RAY_RENDER_H__

