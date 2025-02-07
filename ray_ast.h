#ifndef RAY_AST_H_
#define RAY_AST_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

static inline double dot(const double *a, const double *b, int n) {
	double ret = 0;
	for (int i = 0; i < n; i++)
		ret += a[i] * b[i];
	return ret;
}

typedef struct pt3 {
	double v[3];
} pt3;

typedef struct pt4 {
	double v[4];
} pt4;

static inline double pt3_dot(const pt3 *a, const pt3 *b) {
	return dot(a->v, b->v, 3);
}

static inline double pt3_dotv(const pt3 a, const pt3 b) {
	return dot(a.v, b.v, 3);
}

static inline double pt4_dot(const pt4 *a, const pt4 *b) {
	return dot(a->v, b->v, 4);
}

static inline double mag(const pt3 *a) {
	return sqrt(pt3_dot(a, a));
}

static inline void pt3_normalize_mut(pt3 *a) {
	double m = mag(a);
	if (m > 0)
		for (int i = 0; i < 3; i++)
			a->v[i] /= m;
}

static inline pt3 pt3_normalize(const pt3 *a) {
	pt3 ret = *a;
	double m = mag(a);
	if (m > 0)
		for (int i = 0; i < 3; i++)
			ret.v[i] /= m;
	return ret;
}

static inline pt3 pt3_addv(const pt3 a, const pt3 b) {
	pt3 ret;
	for (int i = 0; i < 3; i++)
		ret.v[i] = a.v[i] + b.v[i];
	return ret;
}

static inline pt3 pt3_add(const pt3 *a, const pt3 *b) {
	pt3 ret;
	for (int i = 0; i < 3; i++)
		ret.v[i] = a->v[i] + b->v[i];
	return ret;
}

static inline pt3 pt3_sub(const pt3 *a, const pt3 *b) {
	pt3 ret;
	for (int i = 0; i < 3; i++)
		ret.v[i] = a->v[i] - b->v[i];
	return ret;
}

static inline pt3 pt3_mul(const pt3 *a, double t) {
	pt3 ret;
	for (int i = 0; i < 3; i++)
		ret.v[i] = a->v[i] * t;
	return ret;
}

static inline pt4 pt4_mul(const pt4 *a, double t) {
	pt4 ret;
	for (int i = 0; i < 4; i++)
		ret.v[i] = a->v[i] * t;
	return ret;
}

static inline void pt4_add_mut(pt4 *a, const pt4 *b) {
	for (int i = 0; i < 4; i++)
		a->v[i] += b->v[i];
}

static inline pt4 pt4_mul_ptwise(const pt4 *a, const pt4 *b) {
	pt4 ret = {{ a->v[0] * b->v[0], a->v[1] * b->v[1], a->v[2] * b->v[2], a->v[3] * b->v[3] }};
	return ret;
}

static inline double pt3_pt3_dist(const pt3 *a, const pt3 *b) {
	pt3 diff = pt3_sub(a, b);
	return mag(&diff);
}

typedef struct ray {
	pt3 origin;
	pt3 direction;
} ray;

static inline pt3 ray_project(const ray *r, double t) {
	pt3 m = pt3_mul(&r->direction, t);
	return pt3_add(&r->origin, &m);
}

typedef struct mat3 {
	double m[9];
} mat3;

static inline pt3 mat3_pt3_mul(const mat3 *m, const pt3 *p) {
	pt3 ret = {
		.v = {
			m->m[0] * p->v[0] + m->m[1] * p->v[1] + m->m[2] * p->v[2],
			m->m[3] * p->v[0] + m->m[4] * p->v[1] + m->m[5] * p->v[2],
			m->m[6] * p->v[0] + m->m[7] * p->v[1] + m->m[8] * p->v[2]
		},
	};
	return ret;
}

static inline mat3 mat3_mat3_mul(const mat3 *a, const mat3 *b) {
	mat3 ret = {
		.m = {
			a->m[0] * b->m[0] + a->m[3] * b->m[1] + a->m[6] * b->m[2],
			a->m[1] * b->m[0] + a->m[4] * b->m[1] + a->m[7] * b->m[2],
			a->m[2] * b->m[0] + a->m[5] * b->m[1] + a->m[8] * b->m[2],
			a->m[0] * b->m[3] + a->m[3] * b->m[4] + a->m[6] * b->m[5],
			a->m[1] * b->m[3] + a->m[4] * b->m[4] + a->m[7] * b->m[5],
			a->m[2] * b->m[3] + a->m[5] * b->m[4] + a->m[8] * b->m[5],
			a->m[0] * b->m[6] + a->m[3] * b->m[7] + a->m[6] * b->m[8],
			a->m[1] * b->m[6] + a->m[4] * b->m[7] + a->m[7] * b->m[8],
			a->m[2] * b->m[6] + a->m[5] * b->m[7] + a->m[8] * b->m[8],
		},
	};
	return ret;
}

static inline mat3 rotation(double yaw, double pitch, double roll) {
	double cosa = cos(yaw);
	double cosb = cos(pitch);
	double cosy = cos(roll);
	double sina = sin(yaw);
	double sinb = sin(pitch);
	double siny = sin(roll);

	mat3 ret = {
		.m = {
			cosa * cosb,	cosa * sinb * siny - sina * cosy,	cosa * sinb * cosy + sina * siny,
			sina * cosb,	sina * sinb * siny + cosa * cosy,	sina * sinb * cosy - cosa * siny,
			-sinb,		cosb * siny,				cosb * cosy,
		},
	};
	return ret;
}

typedef struct color {
	pt4 rgba;
	double reflectance;
} color;

typedef struct sphere {
	color color;
	pt3 position;
	pt3 velocity;
	double mass;
	double radius;
} sphere;

typedef struct plane {
	color color;
	pt3 position;
	pt3 normal;
} plane;

typedef struct light {
	color color;
	pt3 position;
} light;

typedef struct spherespec {
	color *color;
	pt3 *position;
	pt3 *velocity;
	double mass;
	double radius;
	struct spherespec *next;
} spherespec;

typedef struct planespec {
	color *color;
	pt3 *position;
	pt3 *normal;
	struct planespec *next;
} planespec;

typedef struct colorspec {
	pt4 *rgba;
	double reflectance;
	struct colorspec *next;
} colorspec;

typedef struct lightspec {
	color *color;
	pt3 *position;
	struct lightspec *next;
} lightspec;

typedef struct spherespecs {
	spherespec *first;
	spherespec *last;
} spherespecs;

typedef struct planespecs {
	planespec *first;
	planespec *last;
} planespecs;

typedef struct lightspecs {
	lightspec *first;
	lightspec *last;
} lightspecs;

typedef struct colorspecs {
	colorspec *first;
	colorspec *last;
} colorspecs;

struct context {
	int num_spheres;
	int num_planes;
	int num_lights;
	sphere *spheres;	
	plane *planes;
	light *lights;
};

#define __CONTEXT_APPEND(ctx, x, v)	do {					\
	if (ctx->x == NULL) {							\
		ctx->x = malloc(sizeof(*ctx->x));				\
	} else {								\
		ctx->x = realloc(ctx->x, sizeof(*ctx->x) * (ctx->num_##x + 1));		\
	}									\
	ctx->x[ctx->num_##x++] = v;						\
} while(0)

static inline void context_add_sphere(struct context *ctx, sphere s) {
	__CONTEXT_APPEND(ctx, spheres, s);
}
static inline void context_add_plane(struct context *ctx, plane s) {
	__CONTEXT_APPEND(ctx, planes, s);
}
static inline void context_add_light(struct context *ctx, light s) {
	__CONTEXT_APPEND(ctx, lights, s);
}	

static inline void apply_planespec(plane* p, planespec *spec) {
	if (spec->color != NULL) p->color = *spec->color;
	if (spec->position != NULL) p->position = *spec->position;
	if (spec->normal != NULL) p->normal = *spec->normal;
}

static inline void apply_spherespec(sphere* p, spherespec *spec) {
	if (spec->color != NULL) p->color = *spec->color;
	if (spec->position != NULL) p->position = *spec->position;
	if (spec->velocity != NULL) p->velocity = *spec->velocity;
	if (spec->mass > 0) p->mass = spec->mass;
	if (spec->radius > 0) p->radius = spec->radius;
}

static inline void apply_lightspec(light* p, lightspec *spec) {
	if (spec->color != NULL) p->color = *spec->color;
	if (spec->position != NULL) p->position = *spec->position;
}

static inline void apply_colorspec(color* p, colorspec *spec) {
	if (spec->rgba != NULL) p->rgba = *spec->rgba;
	if (spec->reflectance > 0) p->reflectance = spec->reflectance;
}

#define append_ll(a, b)		do { if (a->first == NULL) { a->first = a->last = b; } else { a->last->next = b; a->last = b; b->next = NULL; } } while(0)
#define prepend_ll(a, b)	do { if (a->first == NULL) { a->first = a->last = b; } else { b->next = a->first; a->first = b; } } while(0)
#define FREE_LL(t, x)		do { struct t *p = x->first; while (p) { struct t *next = p->next; free_##t(p); p = next; } } while(0)

#define CREATE_NEW_FN(x)	static inline struct x *new_##x() { struct x *p = malloc(sizeof(struct x)); memset(p, 0, sizeof(struct x)); return p; }
CREATE_NEW_FN(pt3)
CREATE_NEW_FN(pt4)
CREATE_NEW_FN(spherespec)
CREATE_NEW_FN(spherespecs)
CREATE_NEW_FN(planespec)
CREATE_NEW_FN(planespecs)
CREATE_NEW_FN(lightspec)
CREATE_NEW_FN(lightspecs)
CREATE_NEW_FN(color)
CREATE_NEW_FN(colorspec)
CREATE_NEW_FN(colorspecs)
CREATE_NEW_FN(context)

static inline void free_pt3(pt3 *p) { free(p); }
static inline void free_pt4(pt4 *p) { free(p); }
static inline void free_color(color *p) { free(p); }
static inline void free_spherespec(spherespec *p) {
	if (p->color != NULL) free(p->color);
	if (p->position != NULL) free(p->position);
	if (p->velocity != NULL) free(p->velocity);
	free(p);
}
static inline void free_planespec(planespec *p) {
	if (p->color != NULL) free(p->color);
	if (p->position != NULL) free(p->position);
	if (p->normal != NULL) free(p->normal);
	free(p);
}
static inline void free_lightspec(lightspec *p) {
	if (p->color != NULL) free(p->color);
	if (p->position != NULL) free(p->position);
	free(p);
}
static inline void free_colorspec(colorspec *p) {
	if (p->rgba != NULL) free(p->rgba);
	free(p);
}
static inline void free_spherespecs(spherespecs *ps) {
	FREE_LL(spherespec, ps);
	free(ps);
}
static inline void free_planespecs(planespecs *ps) {
	FREE_LL(planespec, ps);
	free(ps);
}
static inline void free_lightspecs(lightspecs *ls) {
	FREE_LL(lightspec, ls);
	free(ls);
}
static inline void free_colorspecs(colorspecs *ps) {
	FREE_LL(colorspec, ps);
	free(ps);
}
static inline void free_context(struct context *ctx) {
	if (ctx->lights) free(ctx->lights);
	if (ctx->planes) free(ctx->planes);
	if (ctx->spheres) free(ctx->spheres);
	free(ctx);
}

// Hacks here because the lexer and parser are co-dependent for type definitions.
#define YY_TYPEDEF_YY_SCANNER_T
typedef void * yyscan_t;

#endif	// RAY_AST_H_

