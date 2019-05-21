#include "bvh.h"
#include "mesh.h"
#include "world.h"

#define BOUND_MESH 1

#define DEBUG 1

uint32_t xor_shift_u32(uint32_t *state) {
  uint32_t x = *state;
  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 15;
  *state = x;
  return x;
}

BVHNode *bvh_leaf_node(int prim_idx, int axis, Rect3 bounds) {
  BVHNode *node = (BVHNode *) malloc(sizeof(BVHNode));
  node->split_axis = axis;
  node->prim_idx = prim_idx;
  node->bounds = bounds;
  node->left = NULL;
  node->right = NULL;
  node->leaf = true;
  return node;
}

std::vector<BVHPrimitive> bvh_preprocess_world(World *w) {
  std::vector<BVHPrimitive> result;
  int cur_idx = 0;
  // Spheres
  for (int i = 0; i < w->spheres.size(); i++) {
    BVHPrimitive prim;
    prim.idx = cur_idx++;
    prim.p_idx = i;
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
    prim.mesh_idx = i;
    for (int j = 0; j < cur_mesh->f.size(); j++) {
      prim.f_idx = j;
      // Intialize the bounds with one point of the triangle
      prim.bounds = rwm_r3_init_p(cur_mesh->v[cur_mesh->v_idx[cur_v_idx++]]);
      for (int k = 1; k < cur_mesh->f[j]; k++) {
        // Extend the bounding box using the rest of the points
        prim.bounds = rwm_r3_union_p(prim.bounds, cur_mesh->v[cur_mesh->v_idx[cur_v_idx++]]);
      }
      prim.idx = cur_idx++;
      prim.type = PT_TRIANGLE;
      prim.centroid = 0.5 * prim.bounds.min_p + 0.5f * prim.bounds.max_p;
      result.push_back(prim);
    }
  }

  return result;
}

void print_prims(BVHPrimitive *prims, int n) {
#if DEBUG
  for (int i = 0; i < n; i++) {
    printf("idx %d\n", prims[i].idx);
    printf("\tp_idx %d\n", prims[i].p_idx);
    printf("\tmesh_idx %d\n", prims[i].mesh_idx);
    printf("\tf_idx %d\n", prims[i].f_idx);
    printf("\ttype: %d\n", prims[i].type);
    rwm_v3_printf("\tbounds.min_p", &prims[i].bounds.min_p);
    rwm_v3_printf("\tbounds.max_p", &prims[i].bounds.max_p);
    rwm_v3_printf("\tcentroid", &prims[i].centroid);
    puts("");
  }
#endif
}

void bvhn_print(BVHNode *node) {
#if DEBUG
  printf("node.leaf: %d\n", node->leaf);
  if (node->leaf) {
    printf("node.prim_idx: %d\n", node->prim_idx);
  } else {
    printf("node.left: %p\n", node->left);
    printf("node.right: %p\n", node->right);
  }
  rwm_v3_printf("node.bounds.min_p", &node->bounds.min_p);
  rwm_v3_printf("node.bounds.max_p", &node->bounds.max_p);
  puts("");
#endif
}

int x_axis_comp(const void *a, const void *b) {
  BVHPrimitive *left = (BVHPrimitive *) a;
  BVHPrimitive *right = (BVHPrimitive *) b;
  if (left->bounds.min_p.x - right->bounds.min_p.x < 0.0) return -1;
  return 1;
}

int y_axis_comp(const void *a, const void *b) {
  BVHPrimitive *left = (BVHPrimitive *) a;
  BVHPrimitive *right = (BVHPrimitive *) b;
  if (left->bounds.min_p.y - right->bounds.min_p.y < 0.0) return -1;
  return 1;
}

int z_axis_comp(const void *a, const void *b) {
  BVHPrimitive *left = (BVHPrimitive *) a;
  BVHPrimitive *right = (BVHPrimitive *) b;
  if (left->bounds.min_p.z - right->bounds.min_p.z < 0.0) return -1;
  return 1;
}


BVHNode *bvh_recursive_build(BVHPrimitive *prims, int n, int *total_nodes) {
  // Calculate total bounds of the primitives
  // printf("n == %d\n", n);
  BVHNode *node = (BVHNode *) malloc(sizeof(BVHNode));
  node->leaf = false;
  (*total_nodes)++;
  static uint32_t rng_state = 4;
  int axis = xor_shift_u32(&rng_state) % 3;
  // printf("Sorting on axis: %d\n", axis);
  if (axis == 0) {
    qsort(prims, n, sizeof(BVHPrimitive), x_axis_comp);
  } else if (axis == 1) {
    qsort(prims, n, sizeof(BVHPrimitive), y_axis_comp); 
  } else {
    qsort(prims, n, sizeof(BVHPrimitive), z_axis_comp);
  }
  #if 1
  if (n == 1) {
    node->left = bvh_leaf_node(prims[0].idx, axis, prims[0].bounds);
    node->right = NULL;
    (*total_nodes)++;
    bvhn_print(node->left);
  } else if (n == 2) {
    node->left = bvh_leaf_node(prims[0].idx, axis, prims[0].bounds);
    bvhn_print(node->left);
    node->right = bvh_leaf_node(prims[1].idx, axis, prims[1].bounds);
    bvhn_print(node->right);
    *total_nodes += 2;
  } else {
    node->left = bvh_recursive_build(prims, n/2, total_nodes);
    node->right = bvh_recursive_build(prims+n/2, n-n/2, total_nodes);
  }
  node->bounds = rwm_r3_union(node->left->bounds, node->right->bounds);
  bvhn_print(node);
  #endif
  return node;
}

BVHNode *bvh_build(World *world) {
  world->bvh_prims = bvh_preprocess_world(world);
  std::vector<BVHPrimitive> bvh_prims_work_copy = world->bvh_prims;
  printf("total prims: %llu\n", world->bvh_prims.size());
  print_prims(world->bvh_prims.data(), world->bvh_prims.size());
  int total_nodes = 0;
  BVHNode *root = bvh_recursive_build(bvh_prims_work_copy.data(), world->bvh_prims.size(), &total_nodes);
  printf("total nodes: %d\n", total_nodes);
  return root;
}

// Need to return the leaf node of the primitive we intersect with
// Check bb
// Check if children are leaves, if yes, then check them and return
// Else. recurse on children if they exist
bool bvh_intersect(BVHNode *root, Ray *r, BVHNode *out_node) {
  if (!root) return false;
  IntersectInfo ii;
  if (bb_intersect(&root->bounds, r, &ii)) {
    if (root->leaf) {
      out_node = root;
      return true;
    }
    return bvh_intersect(root->left, r, out_node) || 
      bvh_intersect(root->right, r, out_node);
  }
  return false;
}
