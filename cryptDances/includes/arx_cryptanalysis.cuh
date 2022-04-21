#ifndef ARX_CRYPTANALYSIS
#define ARX_CRYPTANALYSIS

#include "constants.cuh"
#include "config.cuh"
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
#include "differential_linear_kernels.cuh"
#include "pnb_kernels.cuh"
#include <string.h>
#include "types.cuh"
#include <math.h>
#include "automatic_linear_expansions.cuh"
#include <time.h>
#include "pnb_kernels.cuh"
#include "cryptanalysis_types.cuh"

extern int my_rank, num_procs;

#endif