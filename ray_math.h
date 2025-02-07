#ifndef RAY_MATH_H__
#define RAY_MATH_H__

#include "ray_ast.h"

int intersect_ray_sphere(const ray *r, const sphere *s, double *t, pt3 *q);
int intersect_ray_plane(const ray *r, const plane *p, double *t, pt3 *q);
int intersect_sphere_sphere(const sphere *a, const sphere *b, pt3 *pos, pt3 *normal);
int intersect_sphere_plane(const sphere *a, const plane *b);

#endif	// RAY_MATH_H__

