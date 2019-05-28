#include "world.h"
#include <rw/rw_transform.h>

#include "bvh.h"

void create_cornell_box(World *world) {

  Mesh *bottom_plane = new Mesh;
  *bottom_plane = mesh_make_plane();
  mesh_attach_transform(bottom_plane, rwtr_trs(
    rwm_v3_init(0.0,0.0,0.0),
    rwm_v3_init(2.0,2.0,2.0),
    -1, 0)
  );
  Mesh *top_plane = new Mesh;
  *top_plane = mesh_make_plane();
  mesh_attach_transform(top_plane, rwtr_trs(
    rwm_v3_init(0.0,2.0,0.0),
    rwm_v3_init(2.0,2.0,2.0),
    0, 180)
  );
  Mesh *back_plane = new Mesh;
  *back_plane = mesh_make_plane();
  mesh_attach_transform(back_plane, rwtr_trs(
    rwm_v3_init(0.0,0.0,-2.0),
    rwm_v3_init(2.0,2.0,2.0),
    0, 90)
  );
  Mesh *left_plane = new Mesh;
  *left_plane = mesh_make_plane();
  mesh_attach_transform(left_plane, rwtr_trs(
    rwm_v3_init(-1.5,0.0,0.0),
    rwm_v3_init(2.0,2.0,2.0),
    2, -90)
  );
  Mesh *right_plane = new Mesh;
  *right_plane = mesh_make_plane();
  mesh_attach_transform(right_plane, rwtr_trs(
    rwm_v3_init(1.5,0.0,0.0),
    rwm_v3_init(2.0,2.0,2.0),
    2, 90)
  );
  world->meshes.push_back(bottom_plane);
  world->meshes.push_back(top_plane);
  world->meshes.push_back(back_plane);
  world->meshes.push_back(left_plane);
  world->meshes.push_back(right_plane);

  world->mesh_materials.push_back({M_DIFFUSE, rwm_v3_init(1.0, 1.0, 1.0), 0, false});
  world->mesh_materials.push_back({M_DIFFUSE, rwm_v3_init(1.0, 1.0, 1.0), 0, false});
  world->mesh_materials.push_back({M_DIFFUSE, rwm_v3_init(1.0, 1.0, 1.0), 0, false});
  world->mesh_materials.push_back({M_DIFFUSE, rwm_v3_init(1.0, 0.0, 0.0), 0, false});
  world->mesh_materials.push_back({M_DIFFUSE, rwm_v3_init(0.0, 1.0, 0.0), 0, false});

  world->lights.push_back({
    LT_SPHERE,
    rwm_v3_init(0.0, 1.9, -1.0), // position
    rwm_v3_init(1.0f, 1.0f, 1.0f), // color
    5.0f, // intensity
  });
}

void create_world(World *world) {
  printf("Loading OBJ files:\n\t");

  Mesh *suz = new Mesh;
  load_obj(suz, "suzanne.obj", false);
  mesh_attach_transform(suz, rwtr_trs(
    rwm_v3_init(-1.0,0.75,-0.80),
    rwm_v3_init(0.5, 0.5, 0.5),
    1, 45.0)
  );
  Mesh *plane = new Mesh;
  *plane = mesh_make_plane();
  mesh_attach_transform(plane, rwtr_trs(
    rwm_v3_init(0.0,0.0,0.0),
    rwm_v3_init(2.0,2.0,2.0),
    0, 0)
  );

  puts("Creating world...");
  // world->spheres.push_back({rwm_v3_init(0,0.25,-0.8), 0.5});
  // world->spheres.push_back({rwm_v3_init(-1,0.25,-0.5), 0.5});
  world->spheres.push_back(sphere_create(rwm_v3_init(0.6,0.12,0.0), 0.25));
  world->spheres.push_back(sphere_create(rwm_v3_init(0,0.25,-0.5), 0.5));
  world->spheres.push_back({rwm_v3_init(1,0.25,-1), 0.5});
  // world->spheres.push_back({rwm_v3_init(-1,0.25,-1), 0.5});
  // world->spheres.push_back({rwm_v3_init(-1,0.25,-2), 1.5});
  world->sphere_materials.push_back({M_RR, rwm_v3_init(1.0, 1.0, 1.0), 1.55});
  world->sphere_materials.push_back({M_DIFFUSE, rwm_v3_init(1.0, 0.0, 0.0)});
  world->sphere_materials.push_back({M_REFLECT, rwm_v3_init(1.0, 1.0, 1.0)});
  // world->sphere_materials.push_back({M_DIFFUSE, rwm_v3_init(1.0, 1.0, 1.0)});

  // world->planes.push_back({rwm_v3_init(0,-0.5,-1), rwm_v3_init(0,1.0,0)});
  // world->plane_materials.push_back({M_DIFFUSE, rwm_v3_init(1.0, 1.0, 1.0)});
  // world->planes.push_back({rwm_v3_init(0,0,-5), rwm_v3_init(0,0,1.0)});
  // world->plane_materials.push_back({M_DIFFUSE, rwm_v3_init(1.0, 1.0, 1.0)});
  // world->tris.push_back({
    // rwm_v3_init(-1,-1.0,-5),
    // rwm_v3_init(1,-1.0,-5),
    // rwm_v3_init(0,1.0,-5)
  // });
  world->meshes.push_back(suz);
  world->meshes.push_back(plane);
  // world->mesh_materials.push_back({M_REFLECT, rwm_v3_zero()});
  world->mesh_materials.push_back({M_DIFFUSE, rwm_v3_init(0.8, 0.8, 0.8), 0, false});
  world->mesh_materials.push_back({M_DIFFUSE, rwm_v3_init(1.0, 1.0, 1.0), 0, true});

  // create_cornell_box(world);

  puts("Adding lights");

  world->lights.push_back({
    LT_SPHERE,
    rwm_v3_init(2.0, 2.0, 1.0), // position
    rwm_v3_init(1.0f, 1.0f, 1.0f), // color
    100.0f, // intensity
  });

  // world->lights.push_back({
  //   LT_SPHERE,
  //   rwm_v3_init(-2.0, 2.0, -3.0), // position
  //   rwm_v3_init(0.0f, 0.0f, 1.0f), // color
  //   600.0f, // intensity
  // });

#if 1
  world->lights.push_back({
      LT_DIRECTION,
      rwm_v3_zero(),
      rwm_v3_init(1.0f, 1.0f, 1.0f), // color
      1.0f, // intensity
      rwm_v3_init(0.0, -1.0, 0.0), // position
  });
#endif

  world->bvh_root = bvh_build(world);
}
