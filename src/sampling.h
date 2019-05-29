#ifndef __SAMPLING_H__
#define __SAMPLING_H__

#include <rw/rw_math.h>

// u is a 2D uniform random sample
Vec3 uniform_sample_hemisphere(Point2 u);
float uniform_hemisphere_pdf();

Vec3 uniform_sample_sphere(Point2 u);
float uniform_sphere_pdf();

Vec3 cosine_sample_hemisphere(Point2 u);
float cosine_hemisphere_pdf(float cos_theta);

Point2 uniform_sample_disk(Point2 u);
Point2 concentric_sample_disk(Point2 u);

Vec3 uniform_sample_cone(Point2 u, float cos_theta_max);
float uniform_cone_pdf(float cos_theta_max);

void stratified_sample_1d(float *out_samples, int n_samples, bool jitter);
void stratified_sample_2d(Point2 *out_samples, int nx, int ny, bool jitter);

// For Multiple Importance Sampling two distributions
float mis_balance_heuristic(int nf, float f_pdf, int ng, float g_pdf);
float mis_power_heuristic(int nf, float f_pdf, int ng, float g_pdf);

#endif
