#include "render.h"
#include <random>
#include <math.h>

#include <rw/rw_math.h>
#define RWTH_IMPLEMENTATION
#include <rw/rw_th.h>

#include "global.h"
#include "texture.h"
#include "utils.h"
#include "metrics.h"

std::default_random_engine generator;
std::uniform_real_distribution<float> distribution(0, 1);
std::uniform_real_distribution<float> spp_distribution(0, 1);

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

bool trace(World *world, Ray *r, IntersectInfo *out_ii) {
  bool did_intersect = false;
#if 1
  Ray bb_ray = *r;
  for (int i = 0; i < world->bvh_prims.size(); i++) {
    if (bb_intersect(&world->bvh_prims[i].bounds, &bb_ray, out_ii)) {
      switch (world->bvh_prims[i].type) {
        case PT_SPHERE: {
          if (sphere_intersect(&(world->spheres[world->bvh_prims[i].idx]), r, out_ii)) {
            rwth_atomic_add_i64((int64_t volatile *) &mtr_num_sphere_isect, 1);
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
            rwth_atomic_add_i64((int64_t volatile *) &mtr_num_triangle_isect, 1);
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
    // "Skybox", just a gradient
    float t = 0.5f*(r->dir.y + 1.0f);
    color = ((1.0f - t)*rwm_v3_init(1.0f, 1.0f, 1.0f) + t * rwm_v3_init(0.5f, 0.7f, 1.0f)) * 0.5f;
    // color = rwm_v3_init(1.0f, 1.0f, 1.0f);
  }
  return color;
}

void create_render_args() {

}

void render(WorkerData *data, Tile t) {
  for (int i = 0; i < TILE_Y; i++) {
    for (int j = 0; j < TILE_X; j++) {
      int x = (int) t.top_right.x + j;
      int y = (int) t.top_right.y + TILE_Y - 1 - i;
      Vec3 color = rwm_v3_zero();
      for (int s = 0; s < SAMPLES_PER_PIXEL; s++) {
        rwth_atomic_add_i64((int64_t volatile *) &mtr_num_primary_rays, 1);
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

void render(RenderArgs *ra) {
  int *cur_data = ra->film;
  for (int i = 0; i < HEIGHT; i++) {
    for (int j = 0; j < WIDTH; j++) {
      int x = j;
      int y = HEIGHT - 1 - i;
      Vec3 color = rwm_v3_zero();
			for (int s = 0; s < SAMPLES_PER_PIXEL; s++) {
        float s1 = SAMPLES_PER_PIXEL > 1 ? spp_distribution(generator) : 0;
        float s2 = SAMPLES_PER_PIXEL > 1 ? spp_distribution(generator) : 0;
				Ray r = camera_get_ray(ra->camera, x+s1, y+s2);
        color += cast_ray(ra->world, &r, 0);
      }
      color = rwm_v3_scalar_div(color, SAMPLES_PER_PIXEL);
      cur_data[WIDTH * i  + j] = create_png_pixel(color);
    }
  }
}
