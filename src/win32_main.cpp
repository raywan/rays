#include <windows.h>
#include <assert.h>
#include <stdio.h>
#include <queue>

#include <rw/rw_time.h>
#define RWM_IMPLEMENTATION
#include <rw/rw_math.h>
#define RWTR_IMPLEMENTATION
#include <rw/rw_transform.h>

#include "render.h"
#include "mesh.h"
#include "global.h"
#include "ray.h"
#include "camera.h"
#include "primitive.h"
#include "world.h"
#include "bvh.h"
#include "utils.h"
#include "metrics.h"

HANDLE jq_mutex;

unsigned int worker(void *t_arg) {
  WorkerData *data = (WorkerData *) t_arg;
  while (1) {
    switch (WaitForSingleObject(jq_mutex, 0xFFFFFFFF)) {
      case WAIT_OBJECT_0: {
        if (data->job_queue->size() == 0) {
          ReleaseMutex(jq_mutex);
          return 1;
        }
        Tile t = data->job_queue->front();
        data->job_queue->pop();
        ReleaseMutex(jq_mutex);
        // Do work
        render(data, t);
      } break;
      case WAIT_ABANDONED:
        return 0;
    }
  }
  return 1;
}

int main(int argc, char *argv[]) {
  print_run_info();

  rwtm_init();
  mtr_start_time = rwtm_now();

  // Initialize the camera
  puts("Initializing camera...");
  int camera_shot = 0;
  if (argc > 1) {
    camera_shot = atoi(argv[1]);
    if (argc == 3) {
      output_name = argv[2];
    }
  }

  Camera camera;
  if (camera_shot == 0) {
#if 0
    camera = camera_init_default();
  } else if (camera_shot == 1) {
#endif
    camera = camera_init(
      rwm_v3_init(0.0, 1.0, 1.0), // position
      rwm_v3_init(0.0, 0.0, -1.0), // target
      rwm_v3_init(0.0, 1.0, 0.0), // up
      90.0f, // fov
      2.0f // aperature
    );
  } else if (camera_shot == 2) {
    camera = camera_init(
      rwm_v3_init(-1.0, 1.0, 1.0), // position
      rwm_v3_init(0.0, 0.0, 0.0), // target
      rwm_v3_init(0.0, 1.0, 0.0), // up
      90.0f, // fov
      2.0f // aperature
    );
  }

  // Intialize scene
  World world;
  create_world(&world);
  world.bvh_root = bvh_build(&world);

  int *data = (int *) malloc(WIDTH * HEIGHT * sizeof(int));
  int *cur_data = data;

#if 1
#if NUM_THREADS != 1
  puts("Begin multithreaded tracing...");

  std::queue<Tile> job_queue;
  construct_tiles(&job_queue);

  HANDLE threads[NUM_THREADS];
  WorkerData worker_data[NUM_THREADS];

  jq_mutex = CreateMutex(NULL, false, NULL);
  if (jq_mutex == NULL) {
    puts("Create jq_mutex failed");
    return 1;
  }

  unsigned long tid;
  for (int i = 0; i < NUM_THREADS; i++) {
    worker_data[i].tid = i;
    worker_data[i].job_queue = &job_queue;
    worker_data[i].film = data;
    worker_data[i].world = &world;
    worker_data[i].camera = &camera;
    threads[i] = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE) worker, &worker_data[i], 0, &tid);
  }
  WaitForMultipleObjects(NUM_THREADS, threads, TRUE, 0xFFFFFFFF);
  for (int i = 0; i < NUM_THREADS; i++) {
    CloseHandle(threads[i]);
  }
  CloseHandle(jq_mutex);
#else
  // Begin tracing
  puts("Begin tracing...");
  RenderArgs ra;
  ra.world = &world;
  ra.camera = &camera;
  ra.film = data;
  render(&ra);

#endif

  puts("Writing results to PNG...");
  int result;
  if (output_name != NULL) {
    result = write_png(output_name, data, WIDTH, HEIGHT);
  } else {
    result = write_png("render.png", data, WIDTH, HEIGHT);
  }
  assert(result == 1);
  puts("");

  print_post_run_metrics();
  puts("Done.");
#endif

  return 0;
}
