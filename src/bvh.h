#ifndef __BVH_H__
#define __BVH_H__

#include <rw/rw_math.h>
#include "primitive.h"

struct World;

struct BVHNode {
  int id;
  Rect3 bounds;
  int prim_idx;
  int split_axis;
  BVHNode *left;
  BVHNode *right;
  bool leaf;
};

struct BVHPrimitive {
  int idx;
  int p_idx;
  PrimitiveType type;
  Rect3 bounds;
  Point3 centroid;
  // For meshes
  int mesh_idx;
  int f_idx; 
};

std::vector<BVHPrimitive> bvh_preprocess_world(World *w);
BVHNode *bvh_build(World *world);
bool bvh_intersect(World *world, BVHNode *root, Ray *orig_ray, IntersectInfo *out_ii, Ray *out_r);

#endif
