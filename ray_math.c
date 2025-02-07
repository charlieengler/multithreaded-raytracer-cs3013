#include "ray_math.h"

// Adapted from https://gamedev.stackexchange.com/questions/96459/fast-ray-sphere-collision-code
// Intersects ray r = p + td, |d| = 1, with sphere s and, if intersecting, 
// returns t value of intersection and intersection point q 
int intersect_ray_sphere(const ray *r, const sphere *s, double *t, pt3 *q)
{
	pt3 m = pt3_sub(&r->origin, &s->position);
	double b = pt3_dot(&m, &r->direction);
	double c = pt3_dot(&m, &m) - s->radius * s->radius;

	// Exit if r’s origin outside s (c > 0) and r pointing away from s (b > 0) 
	if (c > 0.0f && b > 0.0f) return 0; 
	double discr = b*b - c; 

	// A negative discriminant corresponds to ray missing sphere
	if (discr < 0.0f) return 0; 

	// Ray now found to intersect sphere, compute smallest t value of intersection
	double t1 = -b - sqrt(discr);

	// If t is negative, ray started inside sphere so clamp t to zero 
	if (t1 < 0.0f) t1 = 0.0f; 
	*t = t1;
	*q = ray_project(r, t1);

	return 1;
}

int intersect_ray_plane(const ray *r, const plane *p, double *t, pt3 *q) {
	double d = pt3_dot(&r->direction, &p->normal);
	if (d > 0.000001 || d < -0.000001) {
		pt3 origin_diff = pt3_sub(&p->position, &r->origin);
		*t = pt3_dot(&origin_diff, &p->normal) / d;
		if (*t > 0) {
			*q = ray_project(r, *t);
			return 1;
		}
	}

	return 0;
}

int intersect_sphere_sphere(const sphere *a, const sphere *b, pt3 *pos, pt3 *normal) {
	double min_dist = a->radius + b->radius;
	double dist = pt3_pt3_dist(&a->position, &b->position);
	if (dist < min_dist) {
		// collision!
		pt3 direction = pt3_sub(&a->position, &b->position);
		pt3_normalize_mut(&direction);
		*normal = direction;
		ray r = {a->position, direction};
		*pos = ray_project(&r, a->radius);
		return 1;
	}
	return 0;
}

int intersect_sphere_plane(const sphere *s, const plane *p) {
	double distance = pt3_dotv(pt3_sub(&s->position, &p->position), p->normal);
	return distance < (s->radius + 0.001);
}
