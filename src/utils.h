#ifndef __UTILS_H__
#define __UTILS_H__

#include <rw/rw_math.h>

int create_png_pixel(Vec3 color);
int write_png(const char *filename, int *data, int w, int h);

#endif
