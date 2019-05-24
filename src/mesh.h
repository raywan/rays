#ifndef __MESH_H__
#define __MESH_H__

#include <vector>
#include <rw/rw_math.h>
#include <rw/rw_transform.h>

#define NUM_PACKED_ELEMENTS 8

//  TODO(ray): Remove any use of std::vector!!!!!!!!!!!!!!!!

enum class FaceFormat {
  UNDEFINED,
  V,
  VVT,
  VVN,
  VVTVN
};

typedef struct Mesh {
  std::vector<int> f;
  std::vector<int> v_idx;
  std::vector<int> uv_idx;
  std::vector<int> n_idx;
  std::vector<int> t_idx;
  // Raw data from the OBJ file
  std::vector<Vec3> v;
  std::vector<Vec2> uv;
  std::vector<Vec3> n;
  std::vector<Vec3> t;
  // Processed face data, index aligned and for GL_ARRAY_BUFFER rendering
  // NOTE(ray): This is useful for debugging, but we'll be using the packed data
  // most of the time
  std::vector<Vec3> buf_v;
  std::vector<Vec2> buf_uv;
  std::vector<Vec3> buf_n;
  std::vector<Vec3> buf_t;

  // For EBO
  float *packed;

  FaceFormat format;

  Transform obj_world_tr;
} Mesh;

uint8_t load_obj(Mesh *out_mesh, const char *obj_path, bool pack_data);

Rect3 mesh_get_bounds(Mesh *m);

Mesh mesh_make_plane();

void mesh_attach_transform(Mesh *mesh, Transform obj_world_tr);

#endif
