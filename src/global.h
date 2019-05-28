#ifndef __GLOBAL_H__
#define __GLOBAL_H__

#define WIDTH 640
#define HEIGHT 480
#define ASPECT ((float)WIDTH/(float)HEIGHT)

#define EPSILON 0.00001f
#define BIAS 0.0001f
#define REFRACT_BIAS 0.00001f

#define MAX_DEPTH 3
#define SAMPLES_PER_PIXEL 2048

#define USE_GLOBAL_ILLUMINATION 1
#define NUM_PT_SAMPLES 1

#define TILE_X 1
#define TILE_Y 1
#define NUM_THREADS 8

#endif
