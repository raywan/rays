#include "metrics.h"
#include <stdio.h>
#include <stdint.h>
#include "global.h"

#define RWTM_IMPLEMENTATION
#include <rw/rw_time.h>

// Metrics
char *output_name;
uint64_t mtr_start_time;
uint64_t mtr_num_primary_rays;
uint64_t mtr_num_sphere_tests;
uint64_t mtr_num_sphere_isect;
uint64_t mtr_num_plane_tests;
uint64_t mtr_num_plane_isect;
uint64_t mtr_num_triangle_tests;
uint64_t mtr_num_triangle_isect;

void print_run_info() {
  puts("Hello Rays.");
  puts("========================================================================");
  printf("NUM THREADS:              %d\n", NUM_THREADS);
  printf("WIDTH:                    %d\n", WIDTH);
  printf("HEIGHT:                   %d\n", HEIGHT);
  printf("MAX DEPTH:                %d\n", MAX_DEPTH);
  printf("SAMPLES PER PIXEL:        %d\n", SAMPLES_PER_PIXEL);
  printf("GLOBAL ILLUMINATION:      %d\n", USE_GLOBAL_ILLUMINATION);
  printf("NUM PATH TRACING SAMPLES: %d\n", NUM_PT_SAMPLES);
  if (output_name != NULL) {
    printf("OUTPUT FILE:              %s\n", output_name);
  } else {
    printf("OUTPUT FILE:              render.png\n");
  }
  puts("========================================================================");
}

void print_post_run_metrics() {
  printf("Render time (rwtm):                     %fs\n", rwtm_to_sec(rwtm_since(mtr_start_time)));
  printf("Total number of primary rays:           %llu\n", mtr_num_primary_rays);
  printf("Total number of sphere tests:           %llu\n", mtr_num_sphere_tests);
  printf("Total number of sphere intersections:   %llu\n", mtr_num_sphere_isect);
  printf("Total number of plane tests:            %llu\n", mtr_num_plane_tests);
  printf("Total number of plane intersections:    %llu\n", mtr_num_plane_isect);
  printf("Total number of triangle tests:         %llu\n", mtr_num_triangle_tests);
  printf("Total number of triangle intersections: %llu\n", mtr_num_triangle_isect);
}
