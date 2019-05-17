#ifndef __WORLD_H__
#define __WORLD_H__

#include <vector>
#include <rw/rw_math.h>
#include "primitive.h"
#include "bvh.h"

struct Mesh;

enum LightType {
  LT_SPHERE,
  LT_DIRECTION,
};

struct Light {
  LightType type;
  Vec3 position;
  Vec3 color;
  float intensity;
  // For directional light
  Vec3 dir;
};

struct World {
  std::vector<Sphere> spheres;
  std::vector<Material> sphere_materials;

  std::vector<Plane> planes;
  std::vector<Material> plane_materials;

  std::vector<Mesh *> meshes;
  std::vector<Material> mesh_materials;

  std::vector<Triangle> tris;

  std::vector<Light> lights;
  std::vector<BVHPrimitive> bvh_prims;
};

void create_world(World *w);

#endif
