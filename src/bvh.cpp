#include "bvh.h"
#include "mesh.h"
#include "world.h"

#define BOUND_MESH 1

uint32_t xor_shift_u32(uint32_t *state) {
  uint32_t x = *state;
  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 15;
  *state = x;
  return x;
}

std::vector<BVHPrimitive> bvh_preprocess_world(World *w) {
  std::vector<BVHPrimitive> result;
  // Spheres
  for (int i = 0; i < w->spheres.size(); i++) {
    BVHPrimitive prim;
    prim.idx = i;
    prim.type = PT_SPHERE;
    prim.bounds = sphere_get_bounds(&w->spheres[i]);
    prim.centroid = 0.5 * prim.bounds.min_p + 0.5f * prim.bounds.max_p;
   result.push_back(prim);
  }
  // Meshes
  for (int i = 0; i < w->meshes.size(); i++) {
    // Loop through each triangle in the mesh
    int cur_v_idx = 0;
    Mesh *cur_mesh = w->meshes[i];
    BVHPrimitive prim;
    prim.idx = i;
    for (int j = 0; j < cur_mesh->f.size(); j++) {
      prim.f_idx = j;
      // Intialize the bounds with one point of the triangle
      prim.bounds = rwm_r3_init_p(cur_mesh->v[cur_mesh->v_idx[cur_v_idx++]]);
      for (int k = 1; k < cur_mesh->f[j]; k++) {
        // Extend the bounding box using the rest of the points
        prim.bounds = rwm_r3_union_p(prim.bounds, cur_mesh->v[cur_mesh->v_idx[cur_v_idx++]]);
      }
      prim.type = PT_TRIANGLE;
      prim.centroid = 0.5 * prim.bounds.min_p + 0.5f * prim.bounds.max_p;
      result.push_back(prim);
    }
  }

  return result;
}

BVHNode *bvh_recursive_build(std::vector<BVHPrimitive> *prims, int lo, int hi) {
  // Calculate total bounds of the primitives
  BVHNode *node = (BVHNode *) malloc(sizeof(BVHNode));
  // NOTE(ray): I don't think i need to do this computation here
  Rect3 total_bounds = (*prims)[lo].bounds;
  for (int i = lo + 1; i < hi; i++) {
    total_bounds = rwm_r3_union(total_bounds, (*prims)[i].bounds);
  }

  int n_primitives = hi - lo;

  if (n_primitives == 1) {
    node->leaf = true;
    node->left = NULL;
    node->right = NULL;
    node->bounds = total_bounds;
    node->prim_idx = lo;
  } else {
    // Partition based on the centroid
    Rect3 centroid_bounds = rwm_r3_init_p((*prims)[lo].centroid);
    for (int i = lo + 1; i < hi; i++) {
      centroid_bounds = rwm_r3_union_p(centroid_bounds, (*prims)[i].centroid);
    }
    int axis = rwm_r3_max_extent(centroid_bounds);
    int mid = (hi - lo)/2;
    if (centroid_bounds.max_p.e[axis] == centroid_bounds.min_p.e[axis]) {
      node->leaf = true;
      node->left = NULL;
      node->right = NULL;
      node->bounds = total_bounds;
      node->prim_idx = lo;
    } else {
      node->left = bvh_recursive_build(prims, lo, mid);
      node->right = bvh_recursive_build(prims, mid, hi);
      node->bounds = rwm_r3_union(node->left->bounds, node->right->bounds);
      node->split_axis = axis;
    }
  }
  return node;
}

BVHNode *bvh_build(World *world) {
  std::vector<BVHPrimitive> prims = bvh_preprocess_world(world);
  BVHNode *root = bvh_recursive_build(&prims, 0, prims.size());
  return root;
}

bool bvh_intersect(Ray *r, IntersectInfo *ii) {
  return false;
}
