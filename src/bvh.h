#ifndef __BVH_H__
#define __BVH_H__

#include <rw/rw_math.h>

struct BVHNode {
  Rect3 bounds;
  int split_axis;
  int left;
  int right;
  bool leaf;
};

enum PrimitiveType {
  PT_SPHERE,
  PT_TRIANGLE,
  PT_PLANE,
  PT_MESH,
};


struct BVHPrimitive {
  int idx;
  int f_idx; // for meshes
  PrimitiveType type;
  Rect3 bounds;
  Point3 centroid;
};

std::vector<BVHPrimitive> bvh_preprocess_world(World *w);

#endif
