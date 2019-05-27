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

void stratified_sample_1d(float *out_samples, int n_samples, bool jitter);
void stratified_sample_2d(Point2 *out_samples, int nx, int ny, bool jitter);

#endif