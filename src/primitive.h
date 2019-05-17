#ifndef __PRIMITIVE_H__
#define __PRIMITIVE_H__

#include <rw/rw_math.h>
#include "ray.h"
#include "mesh.h"

struct Sphere {
  Vec3 position; // world space position
  float radius;
};

struct Triangle {
  Vec3 v0, v1, v2;
  float u,v,w;
};

struct Plane {
  Vec3 p;
  Vec3 n;
  float w, h;
};

enum PrimitiveType {
  PT_SPHERE,
  PT_TRIANGLE,
  PT_PLANE,
  PT_MESH,
};

// TODO(ray): Move this elsewhere
enum MaterialType {
  M_DIFFUSE = 1 << 0,
  M_REFLECT = 1 << 1,
  M_REFRACT = 1 << 2,
  M_RR = 1 << 3,
};

struct Material {
  MaterialType type;
  Vec3 albedo;
  float ior;
  bool use_texture;
};

struct IntersectInfo {
  int id;
  Vec3 normal;
  Vec3 hit_point;
  Vec3 color;
  Vec2 tex_coord;
  Material *material;
};

Sphere sphere_create(Vec3 world_pos, float radius);
bool sphere_intersect(Sphere *s, Ray *r, IntersectInfo *out_ii);
Rect3 sphere_get_bounds(Sphere *s);
bool triangle_intersect(Triangle *tri, Ray *r, IntersectInfo *out_ii);
bool plane_intersect(Plane *plane, Ray *r, IntersectInfo *out_ii);
bool disk_intersect(Plane *plane, float radius, Ray *r, IntersectInfo *out_ii);
bool mesh_intersect(Mesh *mesh, Ray *r, IntersectInfo *out_ii);
bool bb_intersect(Rect3 *bb, Ray *r, IntersectInfo *out_ii);


#endif
