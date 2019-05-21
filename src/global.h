#ifndef __GLOBAL_H__
#define __GLOBAL_H__

#define WIDTH 640
#define HEIGHT 480
#define ASPECT ((float)WIDTH/(float)HEIGHT)

#define EPSILON 0.00001f
#define BIAS 0.0001f
#define REFRACT_BIAS 0.00001f

#define MAX_DEPTH 2
#define SAMPLES_PER_PIXEL 4

#define USE_GLOBAL_ILLUMINATION 0
#define NUM_PT_SAMPLES 32

#define TILE_X 16
#define TILE_Y 16
#define NUM_THREADS 1

#endif
