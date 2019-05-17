#ifndef __BVH_H__
#define __BVH_H__

#include <rw/rw_math.h>
#include "primitive.h"

struct World;

struct BVHNode {
  Rect3 bounds;
  int prim_idx;
  int split_axis;
  BVHNode *left;
  BVHNode *right;
  bool leaf;
};

struct BVHPrimitive {
  int idx;
  int f_idx; // for meshes
  PrimitiveType type;
  Rect3 bounds;
  Point3 centroid;
};

std::vector<BVHPrimitive> bvh_preprocess_world(World *w);
BVHNode *bvh_build(World *world);
bool bvh_intersect(Ray *r, IntersectInfo *ii);

#endif
