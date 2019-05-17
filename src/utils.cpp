#include "utils.h"
#include <rw/rw_math.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image_write.h>

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
