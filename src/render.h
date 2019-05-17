#ifndef __RENDER_H__
#define __RENDER_H__

#include <queue>
#include "camera.h"
#include "world.h"

struct RenderArgs {
  World *world;
  Camera *camera;
  int *film;
};

struct Tile {
  Vec2 top_right;
};

struct WorkerData {
  int tid;
  std::queue<Tile> *job_queue;
  RenderArgs render_args;
  World *world;
  Camera *camera;
  int *film;
};

void construct_tiles(std::queue<Tile> *jq);
void render(WorkerData *data, Tile t);
void render(RenderArgs *ra);

#endif
