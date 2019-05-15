#ifndef __WORLD_H__
#define __WORLD_H__

struct Mesh;

struct World {
  std::vector<Sphere> spheres;
  std::vector<Material> sphere_materials;

  std::vector<Plane> planes;
  std::vector<Material> plane_materials;

  std::vector<Mesh *> meshes;
  std::vector<Material> mesh_materials;

  std::vector<Triangle> tris;

  std::vector<Light> lights;
  std::vector<BVHPrimitive> bvh_prims;
};

#endif
