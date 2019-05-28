#ifndef __CAMERA_H__
#define __CAMERA_H__

#include <rw/rw_math.h>
#include "ray.h"

struct Camera {
  Vec3 position;
  Vec3 target;
  Vec3 up;
  float fov;

  // Computed
  Vec3 lower_left;
  Vec3 horizontal;
  Vec3 vertical;
  float lens_radius;
  float focus_dist;
};

Camera camera_init(Vec3 position, Vec3 target, Vec3 up, float fov, float aperature);
Camera camera_init_default();
Camera create_scene_camera(int camera_shot);

// For simplicity, instead of using PBRT camera, use the RTiaW camera
Ray camera_get_ray(Camera *camera, float res_x, float res_y);


#endif
