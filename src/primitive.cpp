#include "primitive.h"
#include <stdlib.h>
#include <rw/rw_th.h>

#include "ray.h"
#include "global.h"
#include "metrics.h"
#include "mesh.h"

bool solve_quadratic(const float a, const float b, const float c, float *x0, float *x1) {
  float discriminant = b * b - 4.0 * a * c;
  if (discriminant < 0) return false;
  else if (discriminant == 0.0) {
    *x0 = 0.5f * b/a;
    *x1 = *x0;
  } else {
    float q = (b > 0) ?
      -0.5f * (b + rwm_sqrt(discriminant)) :
      -0.5f * (b - rwm_sqrt(discriminant));
    *x0 = q / a;
    *x1 = c / q;
  }

  // We want the min to be at x0
  if (*x0 > *x1) std::swap(*x0, *x1);

  return true;
}

Sphere sphere_create(Vec3 world_pos, float radius) {
  Sphere result;
  result.radius = radius;
  result.position = world_pos;
  return result;
}

bool sphere_intersect(Sphere *s, Ray *r, IntersectInfo *out_ii) {
  rwth_atomic_add_i64((int64_t volatile *) &mtr_num_sphere_tests, 1);
  Vec3 l = r->origin - s->position;
  float a = rwm_v3_inner(r->dir, r->dir);
  float b = 2.0f * rwm_v3_inner(r->dir, l);
  float c = rwm_v3_inner(l, l) - SQUARE(s->radius);
  float t0, t1;
  if (!solve_quadratic(a, b, c, &t0, &t1)) {
    return false;
  }

  if (t0 < 0 && t1 < 0) return false;

  // We want to return the closest t
  t0 = MIN(t0, t1);
  if (t0 < 0 || r->at_t < t0) return false;
  r->at_t = t0;

#if 0
  Point3 hit = ray_at_t(r);
  float phi = atan2(hit.y, hit.x);
  float u = phi / rwm_to_radians(360.0f);
  float theta = acos(rwm_clamp(-1.0, hit.z/s->radius, 1.0));
  float v = (theta - acos(-1.0f)) / (acos(1.0f) - acos(-1.0f));
  out_ii->tex_coord = rwm_v2_init(u, v);
#endif

  out_ii->normal = rwm_v3_normalize(ray_at_t(r) - s->position);
  out_ii->hit_point = ray_at_t(r);
  return true;
}

Rect3 sphere_get_bounds(Sphere *s) {
  Rect3 result = rwm_r3_init_v3(
    s->position - rwm_v3_init(s->radius, s->radius, s->radius),
    s->position + rwm_v3_init(s->radius, s->radius, s->radius)
  );
  return result;
}

bool triangle_intersect(Triangle *tri, Ray *r, IntersectInfo *out_ii) {
  rwth_atomic_add_i64((int64_t volatile *) &mtr_num_triangle_tests, 1);
  // Moller-Trumbore algorithm
  Vec3 e0 = rwm_v3_subtract(tri->v1, tri->v0);
  Vec3 e1 = rwm_v3_subtract(tri->v2, tri->v0);
  Vec3 pv = rwm_v3_cross(r->dir, e1);
  float det = rwm_v3_inner(pv, e0);

#if defined(CULL_BACK)
  // If det is negative then the triangle is backfacing
  // If it is close to 0, then it misses the triangle
  if (det < EPSILON) return false;
#else
  // ray and triangle are parallel
  if (ABS(det) < EPSILON) return false;
#endif

  float inv_det = 1.0f/det;
  Vec3 tv = r->origin - tri->v0;
  tri->u = rwm_v3_inner(tv, pv) * inv_det;
  if (tri->u < 0 || tri->u > 1) return false;

  Vec3 qv = rwm_v3_cross(tv, e0);
  tri->v = rwm_v3_inner(r->dir, qv) * inv_det;
  if (tri->v < 0 || tri->u + tri->v > 1) return false;

  tri->w = 1 - tri->u - tri->v;

  float t = rwm_v3_inner(e1, qv) * inv_det;
  if (r->at_t < t || t < 0) return false;
  r->at_t = t;

  // DEBUG
  out_ii->color = rwm_v3_init(tri->u, tri->v, tri->w);

  return true;
}

bool plane_intersect(Plane *plane, Ray *r, IntersectInfo *out_ii) {
  rwth_atomic_add_i64((int64_t volatile *) &mtr_num_plane_tests, 1);
  float den = rwm_v3_inner(plane->n, r->dir);
  if (ABS(den) > EPSILON) {
    float t = rwm_v3_inner(rwm_v3_subtract(plane->p, r->origin), plane->n)/den;
    // Need to check if we're intersecting the plane from behind
    // Imagine extending the ray behind the camera
    if (t <= 0 || r->at_t < t) return false;
    r->at_t = t;
    out_ii->normal = plane->n;
    return true;
  }
  return false;
}

bool disk_intersect(Plane *plane, float radius, Ray *r, IntersectInfo *out_ii) {
  // Store the previous t in case the ray hits the plane but not inside the disk
  float prev_t = r->at_t;
  float den = rwm_v3_inner(r->dir, plane->n);
  // Check if the ray hits the plane
  if (plane_intersect(plane, r, out_ii)) {
    // Check if it is inside the disk
    if (rwm_v3_length(rwm_v3_subtract(ray_at_t(r), plane->p)) <= radius) {
      return true;
    }
  }
  // We might have incorrectly set t by the previous intersection test
  // So set it back to the last valid value
  r->at_t = prev_t;
  return false;
}

bool mesh_intersect(Mesh *mesh, Ray *r, IntersectInfo *out_ii) {
  bool did_intersect = false;
  for (int j = 0; j < mesh->f.size(); j++) {
    Triangle triangle;
    triangle.v0 = mesh->v[mesh->v_idx[j * 3]];
    triangle.v1 = mesh->v[mesh->v_idx[j * 3 + 1]];
    triangle.v2 = mesh->v[mesh->v_idx[j * 3 + 2]];
    if (triangle_intersect(&triangle, r, out_ii)) {
      rwth_atomic_add_i64((int64_t volatile *) &mtr_num_triangle_isect, 1);
      // Get surface properties: texture coordinates, vertex normal
      Vec3 col = {triangle.u , triangle.v, triangle.w};
      out_ii->color = col;
      Vec2 uv0 = mesh->uv[mesh->uv_idx[j * 3]];
      Vec2 uv1 = mesh->uv[mesh->uv_idx[j * 3 + 1]];
      Vec2 uv2 = mesh->uv[mesh->uv_idx[j * 3 + 2]];
      out_ii->tex_coord = (triangle.w * uv0) + (triangle.u * uv1) + (triangle.v * uv2);
      Vec3 n0 = mesh->n[mesh->n_idx[j * 3]];
      Vec3 n1 = mesh->n[mesh->n_idx[j * 3 + 1]];
      Vec3 n2 = mesh->n[mesh->n_idx[j * 3 + 2]];
      out_ii->normal = (triangle.w * n0) + (triangle.u * n1) + (triangle.v * n2);

      did_intersect = true;
    }
  }
  return did_intersect;
}

bool bb_intersect(Rect3 *bb, Ray *r, IntersectInfo *out_ii) {
  #if 0
  // TODO(ray): The optimized version doesn't work properly...
  // Something is wrong. Come back to this later.
  float t_min = (bb->p[r->sign[0]].x - r->origin.x) * r->inv_dir.x;
  float t_max = (bb->p[1-r->sign[0]].x - r->origin.x) * r->inv_dir.x;
  float t_ymin = (bb->p[r->sign[1]].y - r->origin.y) * r->inv_dir.y;
  float t_ymax = (bb->p[1-r->sign[1]].y - r->origin.y) * r->inv_dir.y;
  if (t_min > t_ymin || t_ymin > t_max) return false;
  if (t_ymin > t_min) t_min = t_ymin;
  if (t_ymax < t_max) t_max = t_ymax;
  float t_zmin = (bb->p[r->sign[2]].z - r->origin.z) * r->inv_dir.z;
  float t_zmax = (bb->p[1-r->sign[2]].z - r->origin.z) * r->inv_dir.z;
  if (t_min > t_zmax || t_zmin > t_max) return false;
  if (t_zmin > t_min) t_min = t_zmin;
  if (t_zmax < t_max) t_max = t_zmax;
  return (t_min < FLT_MAX) && (t_max > 0);
#else
  float t0 = 0;
  float t1 = FLT_MAX;
  for (int i = 0; i < 3; i++) {
    float t_near = (bb->min_p.e[i] - r->origin.e[i]) * r->inv_dir.e[i];
    float t_far = (bb->max_p.e[i] - r->origin.e[i]) * r->inv_dir.e[i];
    if (t_near > t_far) std::swap(t_near, t_far);
    t0 = t_near > t0 ? t_near : t0;
    t1 = t_far < t1 ? t_far : t1;
    if (t0 > t1) return false;
  }
#endif
  return true;
}
