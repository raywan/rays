#ifndef __TEXTURE_H__
#define __TEXTURE_H__

#include <rw/rw_math.h>
#include <math.h>

static inline float get_checkerboard(Vec2 uv) {
  float checker = (fmod(uv.u * 10, 1.0) > 0.5) ^ (fmod(uv.v * 10, 1.0) < 0.5);
  float c = 0.3 * (1 - checker) + 0.7 * checker;
  return c;
  // return checker;
}

#endif
