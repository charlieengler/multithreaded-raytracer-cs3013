
#include "ray_physics.h"
#include "ray_math.h"

static const double framerate = 25;
static const pt3 gravity = {{0, -9.8 / framerate, 0}};

void step_physics(struct context *ctx) {
	// First we loop through all the spheres in the world, collision detect them with everything else,
	// and modify the velocities _in place_. Note that velocities may be manipulated during the render
	// phase as they don't affect position until they are applied.
	for (int i = 0; i < ctx->num_spheres; i++) {
		sphere *si = &ctx->spheres[i];
		pt3 *pi = &si->position;
		pt3 *vi = &si->velocity;
		for (int j = i+1; j < ctx->num_spheres; j++) {
			sphere *sj = &ctx->spheres[j];
			pt3 *pj = &sj->position;
			pt3 *vj = &sj->velocity;
			pt3 pdiff = pt3_sub(pi, pj);
			pt3_normalize_mut(&pdiff);
			pt3 vdiff = pt3_sub(vi, vj);

			// if they're heading away from each other already, this is an unlucky frame
			int approaching = pt3_dot(&vdiff, &pdiff) < 0;
			if (!approaching) continue;

			// check for a collision (position distance less than sum of radii)
			pt3 hit;
			pt3 normal;
			if (!intersect_sphere_sphere(si, sj, &hit, &normal)) continue;

			// find projections of existing velocities onto the collision normal
			// use the radii to find comparative masses. we assume uniform density for simplicity.
			double ri = si->radius;
			double rj = sj->radius;
			double mi = 4.0 / 3.0 * M_PI * ri * ri * ri;
			double mj = 4.0 / 3.0 * M_PI * rj * rj * rj;

			double a = pt3_dotv(pt3_mul(&normal, 2), vdiff) / (1/mi + 1/mj);

			if (0) fprintf(stderr, "collision pre  i %d j %d, vi %lf %lf %lf vj %lf %lf %lf ri %lf rj %lf mi %lf mj %lf a %lf\n", i, j,
					vi->v[0], vi->v[1], vi->v[2],
					vj->v[0], vj->v[1], vj->v[2],
					ri, rj, mi, mj,
					a);

			*vi = pt3_addv(*vi, pt3_mul(&normal, -a / mi));
			*vj = pt3_addv(*vj, pt3_mul(&normal, a / mj));

			if (0) fprintf(stderr, "collision post i %d j %d, vi %lf %lf %lf vj %lf %lf %lf\n", i, j,
					vi->v[0], vi->v[1], vi->v[2],
					vj->v[0], vj->v[1], vj->v[2]);
		}
		for (int j = 0; j < ctx->num_planes; j++) {
			plane *p = &ctx->planes[j];
			if (intersect_sphere_plane(si, p)) {
				double d = pt3_dot(vi, &p->normal);
				if (d < 0) {
					*vi = pt3_addv(*vi, pt3_mul(&p->normal, -1.99 * d));
				}
			} else {
				// Apply a velocity step to all the spheres; discretized acceleration.
				// But not if we're touching the plane, otherwise we sink into it.
				si->velocity = pt3_addv(si->velocity, gravity);
			}
		}
	}

	// Now that velocities are altered in the previous loop, apply the velocities (and gravity) to everything.
	// This can't be run in parallel with rendering as it occurs now because the raytracer is dependent on the
	// positions of each sphere to know how to draw it.
	for (int i = 0; i < ctx->num_spheres; i++) {
		ctx->spheres[i].position = pt3_addv(ctx->spheres[i].position, pt3_mul(&ctx->spheres[i].velocity, 1.0/24));
	}
}

