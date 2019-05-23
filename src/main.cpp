#include <stdio.h>
#include <pthread.h>
#include <queue>

#include <rw/rw_time.h>
#define RWM_IMPLEMENTATION
#include <rw/rw_math.h>

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

pthread_mutex_t jq_mutex;

void *worker(void *t_arg) {
  WorkerData *data = (WorkerData *) t_arg;
  while (1) {
    // Try to take from the job queue
    pthread_mutex_lock(&jq_mutex);
    if (data->job_queue->size() == 0) {
      pthread_mutex_unlock(&jq_mutex);
      break;
    }
    Tile t = data->job_queue->front();
    data->job_queue->pop();
    pthread_mutex_unlock(&jq_mutex);

    // Do work
    render(data, t);
  }
  pthread_exit(NULL);
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
  world.bvh_prims = bvh_preprocess_world(&world);

  int *data = (int *) malloc(WIDTH * HEIGHT * sizeof(int));
  int *cur_data = data;

#if NUM_THREADS != 1
  puts("Begin multithreaded tracing...");
  pthread_mutex_init(&jq_mutex, NULL);
  std::queue<Tile> job_queue;
  construct_tiles(&job_queue);
  pthread_t threads[NUM_THREADS];
  WorkerData worker_data[NUM_THREADS];
  int rc;
  // Create workers
  for (int i = 0; i < NUM_THREADS; i++) {
    worker_data[i].tid = i;
    worker_data[i].job_queue = &job_queue;
    worker_data[i].film = data;
    worker_data[i].world = &world;
    worker_data[i].camera = &camera;
    rc = pthread_create(&threads[i], NULL, worker, (void *) &worker_data[i]);
  }
  for (int i = 0; i < NUM_THREADS; i++) {
    pthread_join(threads[i], NULL);
  }
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

  return 0;
}
