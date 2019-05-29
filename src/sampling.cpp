#include "sampling.h"
#include <rw/rw_math.h>
#include <math.h>
#include <random>
#include <stdint.h>

#define ONE_MINUS_EPSILON (0x1.fffffep-1)
static std::default_random_engine generator;
static std::uniform_real_distribution<float> unif_01(0.0, 1.0);

// u is a 2D uniform random sample
// Samples uniformly from a hemisphere where the z-axis is up
Vec3 uniform_sample_hemisphere(Point2 u) {
  float z = u.e[0];
  float r = rwm_sqrt(MAX(0.0, 1.0 - SQUARE(z)));
  float phi = 2 * PI * u.e[1];
  return rwm_v3_init(r * cosf(phi), r * sinf(phi), z);
}

float uniform_hemisphere_pdf() {
  return INV_2PI;
}

Vec3 cosine_sample_hemisphere(Point2 u) {
  Point2 d = concentric_sample_disk(u);
  float z = rwm_sqrt(MAX(0.0, 1.0 - SQUARE(d.x) - SQUARE(d.y)));
  return rwm_v3_init(d.x, d.y, z);
}

float cosine_hemisphere_pdf(float cos_theta) {
  return cos_theta * INV_PI;
}

Vec3 uniform_sample_sphere(Point2 u) {
  float z = 1 - 2*u.e[0];
  float r = rwm_sqrt(MAX(0.0, 1.0 - SQUARE(z)));
  float phi = 2 * PI * u.e[1];
  return rwm_v3_init(r * cosf(phi), r * sinf(phi), z);
}

float uniform_sphere_pdf() {
  return INV_4PI;
}

Point2 uniform_sample_disk(Point2 u) {
  float r = rwm_sqrt(u.e[0]);
  float theta = 2 * PI * u.e[1];
  return rwm_v2_init(r * cosf(theta), r * sinf(theta));
}

Point2 concentric_sample_disk(Point2 u) {
  // We want the random variable to be [-1, 1]
  Point2 offset = 2 * u - rwm_v2_init(1, 1);
  if (offset.x == 0.0 && offset.y == 0.0) return rwm_v2_init(0, 0);
  float theta;
  float r;
  if (ABS(offset.x) > ABS(offset.y)) {
    r = offset.x;
    theta = (PI/4.0f) * (offset.y / offset.x);
  } else {
    r = offset.y;
    theta = (PI/2.0f) - (PI/4.0f) * (offset.x / offset.y);
  }
  return r * rwm_v2_init(cosf(theta), sinf(theta));
}

Vec3 uniform_sample_cone(Point2 u, float cos_theta_max) {
  float cos_theta = (1 - u.e[0]) + u.e[0] * cos_theta_max;
  float sin_theta = rwm_sqrt(1 - SQUARE(cos_theta));
  float phi = u.e[1] * 2 * PI;
  return rwm_v3_init(cosf(phi) * sin_theta, sinf(phi) * sin_theta, cos_theta);
}

float uniform_cone_pdf(float cos_theta_max) {
  return 1/(2 * PI * (1 - cos_theta_max));
}

void stratified_sample_1d(float *out_samples, int n_samples, bool jitter) {
  float inv_n_samples = 1.0f/n_samples;
  for (int i = 0; i < n_samples; i++) {
    float delta = jitter ? unif_01(generator) : 0.5f;
    out_samples[i] = MIN((i + delta) * inv_n_samples, ONE_MINUS_EPSILON);
  }
}

void stratified_sample_2d(Point2 *out_samples, int nx, int ny, bool jitter) {
  float dx = 1.0f/nx;
  float dy = 1.0f/ny;
  for (int y = 0; y < ny; y++) {
    for (int x = 0; x < nx; x++) {
      float jx = jitter ? unif_01(generator) : 0.5f;
      float jy = jitter ? unif_01(generator) : 0.5f;
      out_samples->x = MIN((x + jx) * dx, ONE_MINUS_EPSILON);
      out_samples->y = MIN((y + jy) * dy, ONE_MINUS_EPSILON);
      ++out_samples;
    }
  }
}

void latin_hypercube(float *samples, int n_samples, int n_dim) {
  float inv_n_samples = 1.0 / n_samples;
  // Generate LHS samples along the diagonal
  for (int i = 0; i < n_samples; i++) {
    for (int j = 0; j < n_dim; j++) {
      float sj = (i + unif_01(generator)) * inv_n_samples;
      samples[n_dim * i + j] = MIN(sj, ONE_MINUS_EPSILON);
    }
  }

  // Perumute
  std::default_random_engine g2;
  for (int i = 0; i < n_samples; i++) {
    for (int j = 0; j < n_dim; j++) {
      std::uniform_int_distribution<int32_t> unif_u32(0, n_samples-j);
      int other = j + unif_u32(g2);
      std::swap(samples[n_dim * j + i], samples[n_dim * other + i]);
    }
  }
}