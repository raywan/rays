#include "ray.h"
#include <math.h>
#include <stdlib.h>
#include <algorithm>

Ray ray_init(Vec3 origin, Vec3 dir) {
  Ray result;
  result.origin = origin;
  result.dir = dir;
  result.inv_dir.x = 1/result.dir.x;
  result.inv_dir.y = 1/result.dir.y;
  result.inv_dir.z = 1/result.dir.z;
  result.sign[0] = result.dir.x < 0;
  result.sign[1] = result.dir.y < 0;
  result.sign[2] = result.dir.z < 0;
  result.at_t = FLT_MAX;
  return result;
}

Vec3 ray_at_t(Ray *r) {
  Vec3 result = rwm_v3_add(r->origin, rwm_v3_scalar_mult(r->at_t, r->dir));
  return result;
}

Vec3 reflect(Vec3 incident, Vec3 normal) {
  return incident - 2 * rwm_v3_inner(incident, normal) * normal;
}

Vec3 refract(Vec3 incident, Vec3 normal, float ior) {
  float cos_i = rwm_clamp(-1.0f, rwm_v3_inner(incident, normal), 1.0f);
  float eta_i = 1;
  float eta_t = ior;
  Vec3 n = normal;
  if (cos_i < 0) {
    cos_i = -cos_i;
  } else {
    std::swap(eta_i, eta_t);
    n = -normal;
  }
  float eta = eta_i/eta_t;
  float k = 1 - eta * eta * (1 - cos_i * cos_i);
  Vec3 result = k < 0 ? rwm_v3_zero() : eta * incident + (eta * cos_i - sqrtf(k)) * n;
  return result;
}

// Returns the amount reflected, kr
// To get the amount transmitted, kt = 1 - kr
float fresnel(Vec3 incident, Vec3 normal, float ior) {
  float cos_i = rwm_clamp(-1.0f, rwm_v3_inner(incident, normal), 1.0f);
  float eta_i = 1;
  float eta_t = ior;
  if (cos_i > 0) std::swap(eta_i, eta_t);
  float sin_t = eta_i/eta_t * rwm_sqrt(MAX(0.0f, 1 - SQUARE(cos_i)));
  if (sin_t >= 1.0f) {
    return 1.0f;
  }

  float cos_t = rwm_sqrt(MAX(0.0f, 1 - SQUARE(sin_t)));
  cos_i = fabsf(cos_i);
  float r_para = ((eta_t * cos_i) - (eta_i * cos_t)) / ((eta_t * cos_i) + (eta_i * cos_t));
  float r_perp = ((eta_i * cos_i) - (eta_t * cos_t)) / ((eta_i * cos_i) + (eta_t * cos_t));
  return (SQUARE(r_para) + SQUARE(r_perp)) / 2;
}
