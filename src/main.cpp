#include <stdio.h>
#include <pthread.h>
#include <unistd.h>
#include <stdlib.h>
#include <vector>
#include <queue>
#include <random>
#include <float.h>
#include <rw/rw_types.h>
#define RWM_IMPLEMENTATION
#include <rw/rw_math.h>
#define RWTM_IMPLEMENTATION
#include <rw/rw_time.h>
#define RWTR_IMPLEMENTATION
#include <rw/rw_transform.h>
#define RWTH_IMPLEMENTATION
#include <rw/rw_th.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image_write.h>

#include "mesh.h"
// #include "bvh.h"

#define SHADOW_ISECT_PRINTF(message, ...) if (r->type == RT_SHADOW) printf(message, ##__VA_ARGS__)

#define WIDTH 640
#define HEIGHT 480
#define ASPECT ((float)WIDTH/(float)HEIGHT)

#define EPSILON 0.00001f
#define BIAS 0.0001f
#define REFRACT_BIAS 0.00001f

#define MAX_DEPTH 2
#define SAMPLES_PER_PIXEL 4

#define USE_GLOBAL_ILLUMINATION 1
#define NUM_PT_SAMPLES 32

#define TILE_X 16
#define TILE_Y 16
#define NUM_THREADS 8

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

struct BVHNode {
  Rect3 bounds;
  int prim_idx;
  int split_axis;
  BVHNode *left;
  BVHNode *right;
  bool leaf;
};

std::default_random_engine generator;
std::uniform_real_distribution<float> distribution(0, 1);
std::uniform_real_distribution<float> spp_distribution(0, 1);
pthread_mutex_t jq_mutex;

// Metrics
static char *output_name = NULL;
static uint64_t num_primary_rays = 0;
static uint64_t num_sphere_tests = 0;
static uint64_t num_sphere_isect = 0;
static uint64_t num_plane_tests = 0;
static uint64_t num_plane_isect = 0;
static uint64_t num_triangle_tests = 0;
static uint64_t num_triangle_isect = 0;

struct Tile {
  Vec2 top_right;
};

enum RayType {
  RT_CAMERA,
  RT_SHADOW,
  RT_REFLECT,
  RT_REFRACT,
};

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

struct Ray {
  RayType type;
  float at_t;
  Vec3 origin;
  Vec3 dir;

  Vec3 inv_dir;
  int sign[3];
};

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

struct IntersectInfo {
  int id;
  Vec3 normal;
  Vec3 hit_point;
  Vec3 color;
  Vec2 tex_coord;
  Material *material;
};

struct Camera {
  Vec3 position;
  Vec3 target;
  Vec3 up;
  float fov;

  // Computed
  Vec3 lower_left;
  Vec3 horizontal;
  Vec3 vertical;
  float lens_radius;
  float focus_dist;
};

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

struct WorkerData {
  int tid;
  std::queue<Tile> *job_queue;
  World *world;
  Camera *camera;
  int *film;
};

void create_coordinate_system(Vec3 normal, Vec3 *out_nt, Vec3 *out_nb) {
  if (std::fabs(normal.x) > std::fabs(normal.y)) {
    *out_nt = rwm_v3_scalar_div(
      rwm_v3_init(normal.z, 0, -normal.x),
      rwm_sqrt(SQUARE(normal.x) + SQUARE(normal.z))
    );
  } else {
    *out_nt = rwm_v3_scalar_div(
      rwm_v3_init(0, -normal.z, normal.y),
      rwm_sqrt(SQUARE(normal.y) + SQUARE(normal.z))
    );
  }
  *out_nb = rwm_v3_cross(normal, *out_nt);
}

Vec3 uniform_sample_hemisphere(float r1, float r2) {
  float sin_theta = rwm_sqrt(1 - SQUARE(r1));
  float phi = 2 * PI * r2;
  float x = sin_theta * cosf(phi);
  float z = sin_theta * sinf(phi);
  return rwm_v3_init(x, r1, z);
}

Vec3 reflect(Vec3 incident, Vec3 normal) {
  return incident - 2 * rwm_v3_inner(incident, normal) * normal;
}

Vec3 refract(Vec3 incident, Vec3 normal, float ior) {
  float cos_i = rwm_clamp(-1.0f, rwm_v3_inner(incident, normal), 1.0f);
  float eta_i = 1;
  float eta_t = ior;
  Vec3 n = normal;
  if (cos_i < 0) {
    cos_i = -cos_i;
  } else {
    std::swap(eta_i, eta_t);
    n = -normal;
  }
  float eta = eta_i/eta_t;
  float k = 1 - eta * eta * (1 - cos_i * cos_i);
  Vec3 result = k < 0 ? rwm_v3_zero() : eta * incident + (eta * cos_i - sqrtf(k)) * n;
  return result;
}

// Returns the amount reflected, kr
// To get the amount transmitted, kt = 1 - kr
float fresnel(Vec3 incident, Vec3 normal, float ior) {
  float cos_i = rwm_clamp(-1.0f, rwm_v3_inner(incident, normal), 1.0f);
  float eta_i = 1;
  float eta_t = ior;
  if (cos_i > 0) std::swap(eta_i, eta_t);
  float sin_t = eta_i/eta_t * rwm_sqrt(MAX(0.0f, 1 - SQUARE(cos_i)));
  if (sin_t >= 1.0f) {
    return 1.0f;
  }

  float cos_t = rwm_sqrt(MAX(0.0f, 1 - SQUARE(sin_t)));
  cos_i = fabsf(cos_i);
  float r_para = ((eta_t * cos_i) - (eta_i * cos_t)) / ((eta_t * cos_i) + (eta_i * cos_t));
  float r_perp = ((eta_i * cos_i) - (eta_t * cos_t)) / ((eta_i * cos_i) + (eta_t * cos_t));
  return (SQUARE(r_para) + SQUARE(r_perp)) / 2;
}

Ray ray_init(Vec3 origin, Vec3 dir) {
  Ray result;
  result.origin = origin;
  result.dir = dir;
  result.inv_dir.x = 1/result.dir.x;
  result.inv_dir.y = 1/result.dir.y;
  result.inv_dir.z = 1/result.dir.z;
  result.sign[0] = result.dir.x < 0;
  result.sign[1] = result.dir.y < 0;
  result.sign[2] = result.dir.z < 0;
  result.at_t = FLT_MAX;
  return result;
}

Camera camera_init(Vec3 position, Vec3 target, Vec3 up, float fov, float aperature) {
  Camera result;
  result.position = position;
  result.target = target;
  result.up = up;
  result.lens_radius = aperature/2.0f;
  // result.focus_dist = rwm_v3_length(position - target);
  // TODO(ray): Add depth of field later
  float focus_dist = 1.0f;
  Vec3 dir = rwm_v3_normalize(position - target);
  Vec3 right = rwm_v3_normalize(rwm_v3_cross(up, dir));
  Vec3 new_up = rwm_v3_cross(dir, right);

  float theta = rwm_to_radians(fov);
  float half_height = tan(theta/2);
  float half_width = ASPECT * half_height;
  // result.lower_left = position - (half_width * right) - (half_height * new_up) - dir;
  result.lower_left = position - (half_width * focus_dist * right) - (half_height * focus_dist * new_up) - dir * focus_dist;
  result.horizontal = 2 * half_width * focus_dist * right;
  result.vertical = 2 * half_height * focus_dist * new_up;
  // result.position = rwm_v3_zero();
  // result.lower_left = rwm_v3_init(-2.0, -1.0, -1.0);
  // result.horizontal = rwm_v3_init(4.0, 0.0, 0.0);
  // result.vertical = rwm_v3_init(0.0, 2.0, 0.0);
  return result;
}

Camera camera_init_default() {
  return camera_init(
    rwm_v3_init(0.0, 0.0, 0.0), // position
    rwm_v3_init(0.0, 0.0, -1.0), // target
    rwm_v3_init(0.0, 1.0, 0.0), // up
    90.0f, // fov
    2.0f // aperature
  );
}

// For simplicity, instead of using PBRT camera, use the RTiaW camera
Ray camera_get_ray(Camera *camera, float res_x, float res_y) {
  // Normalize raster space coordinates
  float u = (float) res_x/(float) WIDTH;
  float v = (float) res_y/(float) HEIGHT;

  Ray result = ray_init(
    camera->position,
    camera->lower_left + (u*camera->horizontal) + (v*camera->vertical) - camera->position
  );
  result.type = RT_CAMERA;
  return result;
}


Vec3 ray_at_t(Ray *r) {
  Vec3 result = rwm_v3_add(r->origin, rwm_v3_scalar_mult(r->at_t, r->dir));
  return result;
}

int create_png_pixel(Vec3 color) {
	int r = (int) (255.99f * rwm_clamp01(color.r));
	int g = (int) (255.99f * rwm_clamp01(color.g));
	int b = (int) (255.99f * rwm_clamp01(color.b));
	return (0xFF << 24) | (b << 16) | (g << 8) | (r << 0);
}

int write_png(const char *filename, int *data, int w, int h) {
  int result = stbi_write_png(filename, w, h, 4, (void *) data, w * 4);
  return result;
}

float get_checkerboard(Vec2 uv) {
  float checker = (fmod(uv.u * 10, 1.0) > 0.5) ^ (fmod(uv.v * 10, 1.0) < 0.5);
  float c = 0.3 * (1 - checker) + 0.7 * checker;
  return c;
  // return checker;
}

bool solve_quadratic(const float a, const float b, const float c, float *x0, float *x1) {
  float discriminant = b * b - 4.0 * a * c;
  if (discriminant < 0) return false;
  else if (discriminant == 0.0) {
    *x0 = 0.5f * b/a;
    *x1 = *x0;
  } else {
    float q = (b > 0) ?
      -0.5f * (b + rwm_sqrt(discriminant)) :
      -0.5f * (b - rwm_sqrt(discriminant));
    *x0 = q / a;
    *x1 = c / q;
  }

  // We want the min to be at x0
  if (*x0 > *x1) std::swap(*x0, *x1);

  return true;
}

Sphere sphere_create(float radius, Vec3 world_pos) {
  Sphere result;
  result.radius = radius;
  result.position = world_pos;
  return result;
}

bool sphere_intersect(Sphere *s, Ray *r, IntersectInfo *out_ii) {
  rwth_atomic_add_i64((int64_t volatile *) &num_sphere_tests, 1);
  Vec3 l = r->origin - s->position;
  float a = rwm_v3_inner(r->dir, r->dir);
  float b = 2.0f * rwm_v3_inner(r->dir, l);
  float c = rwm_v3_inner(l, l) - SQUARE(s->radius);
  float t0, t1;
  if (!solve_quadratic(a, b, c, &t0, &t1)) {
    return false;
  }

  if (t0 < 0 && t1 < 0) return false;

  // We want to return the closest t
  t0 = MIN(t0, t1);
  if (t0 < 0 || r->at_t < t0) return false;
  r->at_t = t0;

#if 0
  Point3 hit = ray_at_t(r);
  float phi = atan2(hit.y, hit.x);
  float u = phi / rwm_to_radians(360.0f);
  float theta = acos(rwm_clamp(-1.0, hit.z/s->radius, 1.0));
  float v = (theta - acos(-1.0f)) / (acos(1.0f) - acos(-1.0f));
  out_ii->tex_coord = rwm_v2_init(u, v);
#endif

  out_ii->normal = rwm_v3_normalize(ray_at_t(r) - s->position);
  out_ii->hit_point = ray_at_t(r);
  return true;
}

Rect3 sphere_get_bounds(Sphere *s) {
  Rect3 result = rwm_r3_init_v3(
    s->position - rwm_v3_init(s->radius, s->radius, s->radius),
    s->position + rwm_v3_init(s->radius, s->radius, s->radius)
  );
  return result;
}

bool triangle_intersect(Triangle *tri, Ray *r, IntersectInfo *out_ii) {
  rwth_atomic_add_i64((int64_t volatile *) &num_triangle_tests, 1);
  // Moller-Trumbore algorithm
  Vec3 e0 = rwm_v3_subtract(tri->v1, tri->v0);
  Vec3 e1 = rwm_v3_subtract(tri->v2, tri->v0);
  Vec3 pv = rwm_v3_cross(r->dir, e1);
  float det = rwm_v3_inner(pv, e0);

#if defined(CULL_BACK)
  // If det is negative then the triangle is backfacing
  // If it is close to 0, then it misses the triangle
  if (det < EPSILON) return false;
#else
  // ray and triangle are parallel
  if (ABS(det) < EPSILON) return false;
#endif

  float inv_det = 1.0f/det;
  Vec3 tv = r->origin - tri->v0;
  tri->u = rwm_v3_inner(tv, pv) * inv_det;
  if (tri->u < 0 || tri->u > 1) return false;

  Vec3 qv = rwm_v3_cross(tv, e0);
  tri->v = rwm_v3_inner(r->dir, qv) * inv_det;
  if (tri->v < 0 || tri->u + tri->v > 1) return false;

  tri->w = 1 - tri->u - tri->v;

  float t = rwm_v3_inner(e1, qv) * inv_det;
  if (r->at_t < t || t < 0) return false;
  r->at_t = t;

  // DEBUG
  out_ii->color = rwm_v3_init(tri->u, tri->v, tri->w);

  return true;
}

bool plane_intersect(Plane *plane, Ray *r, IntersectInfo *out_ii) {
  rwth_atomic_add_i64((int64_t volatile *) &num_plane_tests, 1);
  float den = rwm_v3_inner(plane->n, r->dir);
  if (ABS(den) > EPSILON) {
    float t = rwm_v3_inner(rwm_v3_subtract(plane->p, r->origin), plane->n)/den;
    // Need to check if we're intersecting the plane from behind
    // Imagine extending the ray behind the camera
    if (t <= 0 || r->at_t < t) return false;
    r->at_t = t;
    out_ii->normal = plane->n;
    return true;
  }
  return false;
}

bool disk_intersect(Plane *plane, float radius, Ray *r, IntersectInfo *out_ii) {
  // Store the previous t in case the ray hits the plane but not inside the disk
  float prev_t = r->at_t;
  float den = rwm_v3_inner(r->dir, plane->n);
  // Check if the ray hits the plane
  if (plane_intersect(plane, r, out_ii)) {
    // Check if it is inside the disk
    if (rwm_v3_length(rwm_v3_subtract(ray_at_t(r), plane->p)) <= radius) {
      return true;
    }
  }
  // We might have incorrectly set t by the previous intersection test
  // So set it back to the last valid value
  r->at_t = prev_t;
  return false;
}

bool mesh_intersect(Mesh *mesh, Ray *r, IntersectInfo *out_ii) {
  bool did_intersect = false;
  for (int j = 0; j < mesh->f.size(); j++) {
    Triangle triangle;
    triangle.v0 = mesh->v[mesh->v_idx[j * 3]];
    triangle.v1 = mesh->v[mesh->v_idx[j * 3 + 1]];
    triangle.v2 = mesh->v[mesh->v_idx[j * 3 + 2]];
    if (triangle_intersect(&triangle, r, out_ii)) {
      rwth_atomic_add_i64((int64_t volatile *) &num_triangle_isect, 1);
      // Get surface properties: texture coordinates, vertex normal
      Vec3 col = {triangle.u , triangle.v, triangle.w};
      out_ii->color = col;
      Vec2 uv0 = mesh->uv[mesh->uv_idx[j * 3]];
      Vec2 uv1 = mesh->uv[mesh->uv_idx[j * 3 + 1]];
      Vec2 uv2 = mesh->uv[mesh->uv_idx[j * 3 + 2]];
      out_ii->tex_coord = (triangle.w * uv0) + (triangle.u * uv1) + (triangle.v * uv2);
      Vec3 n0 = mesh->n[mesh->n_idx[j * 3]];
      Vec3 n1 = mesh->n[mesh->n_idx[j * 3 + 1]];
      Vec3 n2 = mesh->n[mesh->n_idx[j * 3 + 2]];
      out_ii->normal = (triangle.w * n0) + (triangle.u * n1) + (triangle.v * n2);

      did_intersect = true;
    }
  }
  return did_intersect;
}

bool bb_intersect(Rect3 *bb, Ray *r, IntersectInfo *out_ii) {
  float t_min = (bb->p[r->sign[0]].x - r->origin.x) * r->inv_dir.x;
  float t_max = (bb->p[1-r->sign[0]].x - r->origin.x) * r->inv_dir.x;
  float t_ymin = (bb->p[r->sign[1]].y - r->origin.y) * r->inv_dir.y;
  float t_ymax = (bb->p[1-r->sign[1]].y - r->origin.y) * r->inv_dir.y;
  if (t_min > t_ymax || t_ymin > t_max) return false;
  t_min = MAX(t_min, t_ymin);
  t_max = MIN(t_max, t_ymax);
  float t_zmin = (bb->p[r->sign[2]].z - r->origin.z) * r->inv_dir.z;
  float t_zmax = (bb->p[1-r->sign[2]].z - r->origin.z) * r->inv_dir.z;
  if (t_min > t_zmax || t_zmin > t_max) return false;
  t_min = MAX(t_min, t_zmin);
  t_max = MIN(t_max, t_zmax);
  r->at_t = t_min;

  return true;
}

bool trace(World *world, Ray *r, IntersectInfo *out_ii) {
  bool did_intersect = false;
#if 1
  Ray bb_ray = *r;
  for (int i = 0; i < world->bvh_prims.size(); i++) {
    if (bb_intersect(&world->bvh_prims[i].bounds, &bb_ray, out_ii)) {
      switch (world->bvh_prims[i].type) {
        case PT_SPHERE: {
          if (sphere_intersect(&(world->spheres[world->bvh_prims[i].idx]), r, out_ii)) {
            rwth_atomic_add_i64((int64_t volatile *) &num_sphere_isect, 1);
            out_ii->material = &world->sphere_materials[world->bvh_prims[i].idx];
            did_intersect = true;
          }
        } break;
        case PT_TRIANGLE: {
          int f_idx = world->bvh_prims[i].f_idx;
          Mesh *cur_mesh  = world->meshes[world->bvh_prims[i].idx];
          Triangle triangle;
          triangle.v0 = cur_mesh->v[cur_mesh->v_idx[f_idx * 3]];
          triangle.v1 = cur_mesh->v[cur_mesh->v_idx[f_idx * 3 + 1]];
          triangle.v2 = cur_mesh->v[cur_mesh->v_idx[f_idx * 3 + 2]];
          if (triangle_intersect(&triangle, r, out_ii)) {
            rwth_atomic_add_i64((int64_t volatile *) &num_triangle_isect, 1);
            // Get surface properties: texture coordinates, vertex normal
            Vec3 col = {triangle.u , triangle.v, triangle.w};
            out_ii->color = col;
            Vec2 uv0 = cur_mesh->uv[cur_mesh->uv_idx[f_idx * 3]];
            Vec2 uv1 = cur_mesh->uv[cur_mesh->uv_idx[f_idx * 3 + 1]];
            Vec2 uv2 = cur_mesh->uv[cur_mesh->uv_idx[f_idx * 3 + 2]];
            out_ii->tex_coord = (triangle.w * uv0) + (triangle.u * uv1) + (triangle.v * uv2);
            Vec3 n0 = cur_mesh->n[cur_mesh->n_idx[f_idx * 3]];
            Vec3 n1 = cur_mesh->n[cur_mesh->n_idx[f_idx * 3 + 1]];
            Vec3 n2 = cur_mesh->n[cur_mesh->n_idx[f_idx * 3 + 2]];
            out_ii->normal = (triangle.w * n0) + (triangle.u * n1) + (triangle.v * n2);
            out_ii->material = &world->mesh_materials[world->bvh_prims[i].idx];
            did_intersect = true;
          }
        } break;
        // case PT_MESH: {
          // if (mesh_intersect(world->meshes[world->bvh_prims[i].idx], r, out_ii)) {
            // out_ii->material = &world->mesh_materials[world->bvh_prims[i].idx];
            // did_intersect = true;
          // }
        // } break;
        default:
          break;
      }
    }
  }
#else
  for (int i = 0; i < world->spheres.size(); i++) {
    if (sphere_intersect(&(world->spheres[i]), r, out_ii)) {
      rwth_atomic_add_i64((int64_t volatile *) &num_sphere_isect, 1);
      out_ii->material = &world->sphere_materials[i];
      did_intersect = true;
      // SHADOW_ISECT_PRINTF("shadow ray isect sphere\n");
    }
  }

  for (int i = 0; i < world->planes.size(); i++) {
    if (plane_intersect(&(world->planes[i]), r, out_ii)) {
      rwth_atomic_add_i64((int64_t volatile *) &num_plane_isect, 1);
    // if (disk_intersect(&(world->planes[i]), 2, r, out_ii)) {
      out_ii->material = &world->plane_materials[i];
      did_intersect = true;
    }
  }

  for (int i = 0; i < world->tris.size(); i++) {
    if (triangle_intersect(&(world->tris[i]), r, out_ii)) {
      rwth_atomic_add_i64((int64_t volatile *) &num_triangle_isect, 1);
      did_intersect = true;
    }
  }

  for (int i = 0; i < world->meshes.size(); i++) {
    if (mesh_intersect(world->meshes[i], r, out_ii)) {
      out_ii->material = &world->mesh_materials[i];
      did_intersect = true;
      // SHADOW_ISECT_PRINTF("shadow ray isect mesh\n");
      // SHADOW_ISECT_PRINTF("shadow ray t: %f\n", r->at_t);
    }
  }
#endif
  out_ii->hit_point = ray_at_t(r);
  return did_intersect;
}

void calculate_light_contribution(Light *l, IntersectInfo *ii, Vec3 *out_light_dir, Vec3 *out_light_intensity, float *out_dist) {
  Vec3 light_dir;
  Vec3 light_intensity;
  float dist_from_light;
  if (l->type == LT_SPHERE) {
    light_dir = ii->hit_point - l->position;
    float r2 = rwm_v3_length_squared(light_dir);
    dist_from_light = rwm_sqrt(r2);
    light_dir = rwm_v3_normalize(light_dir);
    light_intensity = l->intensity/(4.0 * PI * r2) * l->color;
  } else if (l->type == LT_DIRECTION) {
    light_dir = l->dir;
    light_intensity = l->intensity * l->color;
    dist_from_light = FLT_MAX;
  }
  *out_light_dir = light_dir;
  *out_light_intensity = light_intensity;
  *out_dist = dist_from_light;
};

Vec3 cast_ray(World *world, Ray *r, int cur_depth) {
  if (cur_depth > MAX_DEPTH) return rwm_v3_zero();

  Vec3 color = rwm_v3_zero();
  IntersectInfo ii;
  if (trace(world, r, &ii)) {
    // Trace visibility ray to each light source
    switch (ii.material->type) {
      case M_DIFFUSE: {
#if 1
        Vec3 direct_lighting = rwm_v3_zero();
        for (int i = 0; i < world->lights.size(); i++) {
          float dist_from_light;
          Vec3 light_dir;
          Vec3 light_intensity;
          calculate_light_contribution(&world->lights[i], &ii, &light_dir, &light_intensity, &dist_from_light);

          Ray shadow_ray = ray_init(
            ii.hit_point + ii.normal * BIAS,
            -light_dir
          );
          shadow_ray.at_t = dist_from_light;
          shadow_ray.type = RT_SHADOW;

          // Trace back to the light source to check visibility
          IntersectInfo shading_ii;
          if (!trace(world, &shadow_ray, &shading_ii)) {
            float NdV = MAX(0.0f, rwm_v3_inner(ii.normal, -light_dir));
            direct_lighting += light_intensity * NdV;
          }
        }

        Vec3 indirect_lighting = rwm_v3_zero();
#if USE_GLOBAL_ILLUMINATION
        Vec3 Nt, Nb;
        create_coordinate_system(ii.normal, &Nt, &Nb);
        float pdf = 1/(2*PI);

        for (int n = 0; n < NUM_PT_SAMPLES; n++) {
          float r1 = distribution(generator);
          float r2 = distribution(generator);
          Vec3 sample = uniform_sample_hemisphere(r1, r2);
          // Transform the sample into the intersection normal coordinate space
          Vec3 sample_normal_space = rwm_v3_init(
            (sample.x * Nb.x) + (sample.y * ii.normal.x) + (sample.z * Nt.x),
            (sample.x * Nb.y) + (sample.y * ii.normal.y) + (sample.z * Nt.y),
            (sample.x * Nb.z) + (sample.y * ii.normal.z) + (sample.z * Nt.z)
          );
          Ray sample_ray = ray_init(
            ii.hit_point + sample_normal_space * BIAS,
            sample_normal_space
          );
          indirect_lighting += r1 * rwm_v3_scalar_div(cast_ray(world, &sample_ray, cur_depth+1), pdf);
        }
        indirect_lighting = rwm_v3_scalar_div(indirect_lighting, (float) NUM_PT_SAMPLES);
#endif
        color = (rwm_v3_scalar_div(direct_lighting, PI) + 2 * indirect_lighting) * ii.material->albedo;
        if (ii.material->use_texture) {
          color *= get_checkerboard(ii.tex_coord);
          // color = (rwm_v3_scalar_div(direct_lighting, PI) + 2 * indirect_lighting) * get_checkerboard(ii.tex_coord);
        } else {
          // color = (rwm_v3_scalar_div(direct_lighting, PI) + 2 * indirect_lighting) * ii.material->albedo;
        }

#else
        float NdV = MAX(0.0f, rwm_v3_inner(ii.normal, -r->dir));
        // color = get_checkerboard(ii.tex_coord) * ii.normal * NdV;
        color = ii.normal * NdV;
#endif
      } break;
      case M_REFLECT: {
        Vec3 reflect_dir = rwm_v3_normalize(reflect(r->dir, ii.normal));
        Ray reflect_ray = ray_init(
          ii.hit_point + ii.normal * BIAS,
          reflect_dir
        );
        reflect_ray.type = RT_REFLECT;
        color += cast_ray(world, &reflect_ray, cur_depth + 1);
      } break;
      case M_REFRACT: {
        bool outside = rwm_v3_inner(r->dir, ii.normal) < 0;
        Vec3 refract_origin = outside ?
          ii.hit_point - ii.normal * REFRACT_BIAS :
          ii.hit_point + ii.normal * REFRACT_BIAS;
        Vec3 refract_dir = rwm_v3_normalize(refract(r->dir, ii.normal, ii.material->ior));
        Ray refraction_ray = ray_init(refract_origin, refract_dir);
        refraction_ray.type = RT_REFRACT;
        color += cast_ray(world, &refraction_ray, cur_depth + 1);
      } break;
#if 1
      case M_RR: {
        Vec3 refract_color = rwm_v3_zero();
        float kr =  fresnel(r->dir, ii.normal, ii.material->ior);
        bool outside = rwm_v3_inner(r->dir, ii.normal) < 0;
        if (kr < 1) {
          Vec3 refract_origin = outside ?
            ii.hit_point - ii.normal * REFRACT_BIAS :
            ii.hit_point + ii.normal * REFRACT_BIAS;
          Vec3 refract_dir = rwm_v3_normalize(refract(r->dir, ii.normal, ii.material->ior));
          Ray refraction_ray = ray_init(refract_origin, refract_dir);
          refraction_ray.type = RT_REFRACT;
          refract_color += cast_ray(world, &refraction_ray, cur_depth + 1);
        }
        Vec3 reflect_color = rwm_v3_zero();
        Vec3 reflect_origin = outside ?
          ii.hit_point + ii.normal * BIAS :
          ii.hit_point - ii.normal * BIAS;
        Ray reflect_ray = ray_init(
          reflect_origin,
          reflect(r->dir, ii.normal)
        );
        reflect_ray.type = RT_REFLECT;
        reflect_color = cast_ray(world, &reflect_ray, cur_depth + 1);
        color += reflect_color * kr + refract_color * (1 - kr);
      } break;
#endif
      default:
        break;
    }
  // } else if (cur_depth == 0) {
    // color = rwm_v3_init(0.5, 0.7, 1.0);
  } else {
    float t = 0.5f*(r->dir.y + 1.0f);
    color = ((1.0f - t)*rwm_v3_init(1.0f, 1.0f, 1.0f) + t * rwm_v3_init(0.5f, 0.7f, 1.0f)) * 0.5f;
    // color = rwm_v3_init(1.0f, 1.0f, 1.0f);
  }
  return color;
}

void print_run_info() {
  puts("Hello Rays.");
  puts("========================================================================");
  printf("NUM THREADS:              %d\n", NUM_THREADS);
  printf("WIDTH:                    %d\n", WIDTH);
  printf("HEIGHT:                   %d\n", HEIGHT);
  printf("MAX DEPTH:                %d\n", MAX_DEPTH);
  printf("SAMPLES PER PIXEL:        %d\n", SAMPLES_PER_PIXEL);
  printf("GLOBAL ILLUMINATION:      %d\n", USE_GLOBAL_ILLUMINATION);
  printf("NUM PATH TRACING SAMPLES: %d\n", NUM_PT_SAMPLES);
  if (output_name != NULL) {
    printf("OUTPUT FILE:              %s\n", output_name);
  } else {
    printf("OUTPUT FILE:              render.png\n");
  }
  puts("========================================================================");
}

void construct_tiles(std::queue<Tile> *jq) {
  int num_tiles_w = WIDTH/TILE_X;
  int num_tiles_h = HEIGHT/TILE_Y;
  for (int i = 0; i < num_tiles_w; i++) {
    for (int j = 0; j < num_tiles_h; j++) {
      Tile t;
      t.top_right = rwm_v2_init(i*16, j*16);
      jq->push(t);
    }
  }
}

void *worker(void *t_arg) {
  WorkerData *data = (WorkerData *) t_arg;
  while (1) {
    // Try to take from the job queue
    pthread_mutex_lock(&jq_mutex);
    if (data->job_queue->size() == 0) {
      pthread_mutex_unlock(&jq_mutex);
      break;
    }
    Tile t = data->job_queue->front();
    data->job_queue->pop();
    pthread_mutex_unlock(&jq_mutex);

    // Do work
    for (int i = 0; i < TILE_Y; i++) {
      for (int j = 0; j < TILE_X; j++) {
        int x = (int) t.top_right.x + j;
        int y = (int) t.top_right.y + TILE_Y - 1 - i;
        Vec3 color = rwm_v3_zero();
        for (int s = 0; s < SAMPLES_PER_PIXEL; s++) {
          rwth_atomic_add_i64((int64_t volatile *) &num_primary_rays, 1);
          float s1 = SAMPLES_PER_PIXEL > 1 ? spp_distribution(generator) : 0;
          float s2 = SAMPLES_PER_PIXEL > 1 ? spp_distribution(generator) : 0;
          Ray r = camera_get_ray(data->camera, x+s1, y+s2);
          color += cast_ray(data->world, &r, 0);
        }
        color = rwm_v3_scalar_div(color, SAMPLES_PER_PIXEL);
        data->film[(HEIGHT - y) * WIDTH + x] = create_png_pixel(color);
      }
    }
  }
  pthread_exit(NULL);
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

int main(int argc, char *argv[]) {
  rwtm_init();
  uint64_t start_time = rwtm_now();

  // Initialize the camera
  puts("Initializing camera...");
  int camera_shot = 0;
  if (argc > 1) {
    camera_shot = atoi(argv[1]);
    if (argc == 3) {
      output_name = argv[2];
    }
  }

  print_run_info();

  Camera camera;
  if (camera_shot == 0) {
#if 0
    camera = camera_init_default();
  } else if (camera_shot == 1) {
#endif
    camera = camera_init(
      rwm_v3_init(0.0, 1.0, 1.0), // position
      rwm_v3_init(0.0, 0.0, -1.0), // target
      rwm_v3_init(0.0, 1.0, 0.0), // up
      90.0f, // fov
      2.0f // aperature
    );
  } else if (camera_shot == 2) {
    camera = camera_init(
      rwm_v3_init(-1.0, 1.0, 1.0), // position
      rwm_v3_init(0.0, 0.0, 0.0), // target
      rwm_v3_init(0.0, 1.0, 0.0), // up
      90.0f, // fov
      2.0f // aperature
    );
  }

  // Intialize scene
  printf("Loading OBJ files:\n\t");
  Mesh suz;
  load_obj(&suz, "suzanne.obj", false);
  Mesh plane = mesh_make_plane();

  puts("Creating world...");
  World world;
  // world.spheres.push_back((Sphere) {rwm_v3_init(0,0.25,-0.8), 0.5});
  // world.spheres.push_back((Sphere) {rwm_v3_init(-1,0.25,-0.5), 0.5});
  world.spheres.push_back((Sphere) {rwm_v3_init(0.6,0.12,0.0), 0.25});
  world.spheres.push_back((Sphere) {rwm_v3_init(0,0.25,-0.5), 0.5});
  // world.spheres.push_back((Sphere) {rwm_v3_init(1,0.25,-1), 0.5});
  // world.spheres.push_back((Sphere) {rwm_v3_init(-1,0.25,-1), 0.5});
  // world.spheres.push_back((Sphere) {rwm_v3_init(-1,0.25,-2), 1.5});
  world.sphere_materials.push_back((Material) {M_REFLECT, rwm_v3_init(1.0, 1.0, 1.0)});
  // world.sphere_materials.push_back((Material) {M_RR, rwm_v3_init(1.0, 1.0, 1.0), 1.55});
  world.sphere_materials.push_back((Material) {M_DIFFUSE, rwm_v3_init(1.0, 0.0, 0.0)});
  // world.sphere_materials.push_back((Material) {M_DIFFUSE, rwm_v3_init(1.0, 1.0, 1.0)});

  // world.planes.push_back((Plane) {rwm_v3_init(0,-0.5,-1), rwm_v3_init(0,1.0,0)});
  // world.plane_materials.push_back((Material) {M_DIFFUSE, rwm_v3_init(1.0, 1.0, 1.0)});
  // world.planes.push_back((Plane) {rwm_v3_init(0,0,-5), rwm_v3_init(0,0,1.0)});
  // world.plane_materials.push_back((Material) {M_DIFFUSE, rwm_v3_init(1.0, 1.0, 1.0)});
  // world.tris.push_back((Triangle) {
    // rwm_v3_init(-1,-1.0,-5),
    // rwm_v3_init(1,-1.0,-5),
    // rwm_v3_init(0,1.0,-5)
  // });
  // world.meshes.push_back(&suz);
  world.meshes.push_back(&plane);
  // world.mesh_materials.push_back((Material) {M_REFLECT, rwm_v3_zero()});
  world.mesh_materials.push_back((Material) {M_DIFFUSE, rwm_v3_init(1.0, 1.0, 1.0), 0, true});
  puts("Adding lights");

  world.lights.push_back((Light) {
      LT_SPHERE,
      rwm_v3_init(2.0, 2.0, 1.0), // position
      rwm_v3_init(1.0f, 1.0f, 1.0f), // color
      300.0f, // intensity
  });

  // world.lights.push_back((Light) {
      // LT_SPHERE,
      // rwm_v3_init(-2.0, 2.0, -3.0), // position
      // rwm_v3_init(0.0f, 0.0f, 1.0f), // color
      // 600.0f, // intensity
  // });

#if 1
  world.lights.push_back((Light) {
      LT_DIRECTION,
      rwm_v3_zero(),
      rwm_v3_init(1.0f, 1.0f, 1.0f), // color
      1.0f, // intensity
      rwm_v3_init(0.0, -1.0, 0.0), // position
  });
#endif

  world.bvh_prims = bvh_preprocess_world(&world);

  int *data = (int *) malloc(WIDTH * HEIGHT * sizeof(int));
  int *cur_data = data;

#if 1
  puts("Begin multithreaded tracing...");
  pthread_mutex_init(&jq_mutex, NULL);
  std::queue<Tile> job_queue;
  construct_tiles(&job_queue);
  pthread_t threads[NUM_THREADS];
  WorkerData worker_data[NUM_THREADS];
  int rc;
  // Create workers
  for (int i = 0; i < NUM_THREADS; i++) {
    worker_data[i].tid = i;
    worker_data[i].job_queue = &job_queue;
    worker_data[i].film = data;
    worker_data[i].world = &world;
    worker_data[i].camera = &camera;
    rc = pthread_create(&threads[i], NULL, worker, (void *) &worker_data[i]);
  }
  for (int i = 0; i < NUM_THREADS; i++) {
    pthread_join(threads[i], NULL);
  }
#else
	// Begin tracing
  puts("Begin tracing...");
  for (int i = 0; i < HEIGHT; i++) {
    for (int j = 0; j < WIDTH; j++) {
      int x = j;
      int y = HEIGHT - 1 - i;
      Vec3 color = rwm_v3_zero();
			for (int s = 0; s < SAMPLES_PER_PIXEL; s++) {
        float s1 = SAMPLES_PER_PIXEL > 1 ? spp_distribution(generator) : 0;
        float s2 = SAMPLES_PER_PIXEL > 1 ? spp_distribution(generator) : 0;
				Ray r = camera_get_ray(&camera, x+s1, y+s2);
        color += cast_ray(&world, &r, 0);
      }
      color = rwm_v3_scalar_div(color, SAMPLES_PER_PIXEL);
      cur_data[WIDTH * i  + j] = create_png_pixel(color);
    }
  }

  puts("Writing results to PNG...");
#endif

  int result;
  if (output_name != NULL) {
    result = write_png(output_name, data, WIDTH, HEIGHT);
  } else {
    result = write_png("render.png", data, WIDTH, HEIGHT);
  }
  assert(result == 1);
  puts("");

  printf("Render time (rwtm):                     %fs\n", rwtm_to_sec(rwtm_since(start_time)));
  printf("Total number of primary rays:           %llu\n", num_primary_rays);
  printf("Total number of sphere tests:           %llu\n", num_sphere_tests);
  printf("Total number of sphere intersections:   %llu\n", num_sphere_isect);
  printf("Total number of plane tests:            %llu\n", num_plane_tests);
  printf("Total number of plane intersections:    %llu\n", num_plane_isect);
  printf("Total number of triangle tests:         %llu\n", num_triangle_tests);
  printf("Total number of triangle intersections: %llu\n", num_triangle_isect);
  puts("Done.");

  return 0;
}
