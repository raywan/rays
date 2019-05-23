#ifndef __RAY_H__
#define __RAY_H__

#include <rw/rw_math.h>
#include <float.h>

#define SHADOW_ISECT_PRINTF(message, ...) if (r->type == RT_SHADOW) printf(message, ##__VA_ARGS__)

enum RayType {
  RT_CAMERA,
  RT_SHADOW,
  RT_REFLECT,
  RT_REFRACT,
  RT_GI,
};

struct Ray {
  RayType type;
  float at_t;
  Vec3 origin;
  Vec3 dir;

  Vec3 inv_dir;
  int sign[3];
};

Ray ray_init(Vec3 origin, Vec3 dir);
Vec3 ray_at_t(Ray *r);

Vec3 reflect(Vec3 incident, Vec3 normal);
Vec3 refract(Vec3 incident, Vec3 normal, float ior);

// Returns the amount reflected, kr
// To get the amount transmitted, kt = 1 - kr
float fresnel(Vec3 incident, Vec3 normal, float ior);


#endif
