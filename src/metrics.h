#ifndef __METRICS_H__
#define __METRICS_H__

#include <stdint.h>

extern char *output_name;
extern uint64_t mtr_start_time;
extern uint64_t mtr_num_primary_rays;
extern uint64_t mtr_num_sphere_tests;
extern uint64_t mtr_num_sphere_isect;
extern uint64_t mtr_num_plane_tests;
extern uint64_t mtr_num_plane_isect;
extern uint64_t mtr_num_triangle_tests;
extern uint64_t mtr_num_triangle_isect;

void print_run_info();
void print_post_run_metrics();

#endif
