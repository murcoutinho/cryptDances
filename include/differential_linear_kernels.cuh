#ifndef DIFFKH
#define DIFFKH 

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <string.h>
#include <curand.h>
#include <curand_kernel.h>
#include <mpi.h>
#include "util.cuh"
#include "gputimer.cuh"
#include "algorithms.cuh"
#include "types.cuh"
#include "cryptanalysis_types.cuh"

#define TYPE_DIFFERENTIAL 0
#define TYPE_LINEAR 1

void compute_all_single_bit_differential_correlation(int alg_type, int subrounds, int last_subround, 
uint64_t number_of_trials, const char *out_file_name);

void compute_differential_or_linear_correlation(diff_lin_t *diff_lin, int type);
void search_until_find_correlation(diff_lin_t *diff_lin, int type);
#endif