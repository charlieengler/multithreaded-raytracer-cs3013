#ifndef RAY_PHYSICS_H__
#define RAY_PHYSICS_H__

#include "ray_ast.h"
#include "ray_math.h"

void step_physics(struct context *ctx);
void calc_velocities(struct context *ctx);
void update_positions(struct context *ctx);

#endif	// RAY_PHYSICS_H__
