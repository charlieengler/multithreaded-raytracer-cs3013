
#include "ray_render.h"
#include "ray_ast.h"
#include "ray_math.h"

static const pt4 ambient_light = {{0.2, 0.2, 0.2, 1.0}};

// This function should not need to be changed, unless you want to play with the rendering.
int raytrace(const struct context *ctx, const ray *r, pt4 *ret, int depth) {
	double best_t = -1;
	pt3 hit = {0};
	pt3 normal = {0};
	int sphere_hit_index = -1;
	int plane_hit_index = -1;
	for (int i = 0; i < ctx->num_spheres; i++) {
		const struct sphere *const s = &ctx->spheres[i];
		double t;
		pt3 q;
		if (intersect_ray_sphere(r, s, &t, &q)) {
			//printf("sphere t %lf best_t %lf\n", t, best_t);
			if (best_t < 0 || t < best_t) {
				best_t = t;
				hit = q;
				normal = pt3_sub(&q, &s->position);
				pt3_normalize_mut(&normal); 
				sphere_hit_index = i;
			}
			if (0)
			printf("ray: %lf %lf %lf hit sphere %d at %lf %lf %lf\n", r->direction.v[0], r->direction.v[1], r->direction.v[2],
				i, q.v[0], q.v[1], q.v[2] );
			// pr
		}	
	}

	for (int i = 0; i < ctx->num_planes; i++) {
		const struct plane *const p = &ctx->planes[i];
		double t;
		pt3 q;
		if (intersect_ray_plane(r, p, &t, &q)) {
			//printf("plane t %lf best_t %lf\n", t, best_t);
			if (best_t < 0 || t < best_t) {
				best_t = t;
				hit = q;
				normal = p->normal;
				pt3_normalize_mut(&normal); 
				sphere_hit_index = -1;
				plane_hit_index = i;
			}
		}
	}

	// if we have an intersection, apply the color, any maybe recurse.
	const color *c;
	if (sphere_hit_index >= 0) {
		c = &ctx->spheres[sphere_hit_index].color;
	} else if (plane_hit_index >= 0) {
		c = &ctx->planes[plane_hit_index].color;
	} else {
		return 0;
	}

	if (!ret) {
		// this is a light reachability test.
		return 1;
	}

	// add ambient term
	pt4 ambient = pt4_mul_ptwise(&ambient_light, &c->rgba);
	pt4_add_mut(ret, &ambient);
	
	// add a tiny bit to 'hit' to stop zfighting at edges
	pt3 normal_out_bump = pt3_mul(&normal, 0.00001);
	pt3 hit_out_bump = pt3_add(&hit, &normal_out_bump);

	// dot product of the direction to the surface hit normal.
	// negative means we're entering a shape, positive means we're leaving a shape (at least in spheres)
	double nd = pt3_dot(&normal, &r->direction);

	// fire a ray towards light sources.
	if (depth > 0) {
		for (int i = 0; i < ctx->num_lights; i++) {
			const light *const l = &ctx->lights[i];
			pt3 lightdir = pt3_sub(&l->position, &hit);
			pt3_normalize_mut(&lightdir);
			ray rlight = {hit_out_bump, lightdir};
			//printf("light %d to thing\n", i);
			if (!raytrace(ctx, &rlight, NULL, 0)) {
				// add diffuse & maybe specular term for this light
				double light_directness = pt3_dot(&normal, &lightdir);
				//printf("light %d to thing unobstructed, dot %lf\n", i, light_directness);
				if (light_directness > 0) {
					//printf("light directness!\n");
					double light_distance = pt3_pt3_dist(&l->position, &hit);
					double ldist_inv_square = 1.0 / light_distance * light_distance;
					pt4 light_diffuse = pt4_mul(&l->color.rgba, ldist_inv_square * light_directness);
					pt4 surface_diffuse = pt4_mul_ptwise(&light_diffuse, &c->rgba);
					pt4_add_mut(ret, &surface_diffuse);
				}	
			}
		}
	}

	// reflection, if there.
	if (depth > 0 && c->reflectance > 0 && nd < 0) {
		pt3 bounce_normal = pt3_mul(&normal, nd * -2);
		pt3 bounced = pt3_add(&r->direction, &bounce_normal);
		pt3_normalize_mut(&bounced);
		ray bounce_ray = {hit_out_bump, bounced};
		pt4 bounce_color = {0};
		if (raytrace(ctx, &bounce_ray, &bounce_color, depth - 1)) {
			pt4 bounce_color_scaled = pt4_mul(&bounce_color, c->reflectance);
			pt4_add_mut(ret, &bounce_color_scaled);
		}
	}

	// refraction, if there. Alpha channel is see-thru-ness.
	// I couldn't get this working in time.
	if (0 && depth > 0 && c->rgba.v[3] < 1.0) {
		pt3 normal_in_bump = pt3_mul(&normal, -0.00001);
		pt3 hit_in_bump = pt3_add(&hit, &normal_in_bump);

		static const double refractive_proportinality_of_glass = 1.5;	// glass
		double refractive_proportionality = refractive_proportinality_of_glass;
		double angle_of_incidence;
		if (nd > 0) {
			angle_of_incidence = acos(-nd);
			refractive_proportionality = 1.0 / refractive_proportionality;
 		} else {
			angle_of_incidence = acos(nd);
		}
		double angle_of_refraction = angle_of_incidence * refractive_proportionality;	// roughly, for glass
		if (angle_of_refraction < M_PI/2) {
		//printf("nd %lf angle_of_incidence %lf angle_of_refraction %lf\n", nd, angle_of_incidence, angle_of_refraction);
			// lets refract!
			pt3 refract_normal = pt3_mul(&normal, nd * refractive_proportionality);	// this isn't snells law but whatever
			pt3 refracted = pt3_add(&r->direction, &refract_normal);
			pt3_normalize_mut(&refracted);
			ray refract_ray = {hit_in_bump, refracted};
			pt4 refract_color = {0};
			if (raytrace(ctx, &refract_ray, &refract_color, depth - 1)) {
				pt4 refract_color_scaled = pt4_mul(&refract_color, 1.0 - c->rgba.v[3]);
				pt4_add_mut(ret, &refract_color_scaled);
			}		
		}
	}

	return 1;
}


