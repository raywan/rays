#include "camera.h"
#include <math.h>
#include <float.h>

#include "global.h"

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
    43.0f, // fov
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
