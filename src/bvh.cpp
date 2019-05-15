#include "bvh.h"
#include "mesh.h"

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
    prim.type = PT_SPHERE;
    prim.idx = i;
    prim.bounds = sphere_get_bounds(&w->spheres[i]);
    prim.centroid = 0.5 * prim.bounds.min_p + 0.5f * prim.bounds.max_p;
    result.push_back(prim);
  }
  // Meshes
  for (int i = 0; i < w->meshes.size(); i++) {
#if BOUND_MESH
    BVHPrimitive prim;
    prim.idx = i;
    prim.type = PT_MESH;
    prim.bounds = mesh_get_bounds(w->meshes[i]);
    prim.centroid = 0.5 * prim.bounds.min_p + 0.5f * prim.bounds.max_p;
    result.push_back(prim);
#else
    // Loop through each triangle in the mesh
    int cur_v_idx = 0;
    Mesh *cur_mesh = w->meshes[i];
    for (int j = 0; j < cur_mesh->f.size(); j++) {
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
#endif

  return result;
}

BVHNode *bvh_recursive_build() {
}

BVHNode *bvh_build(World *world) {
  // Process world, loop over all primitives and meshes to create BVHPrimitives
  // bvh_recursive_build();
}

bool bvh_intersect(Ray *r, IntersectInfo *ii) {

}
