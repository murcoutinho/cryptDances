#ifndef PNB_KERNELS_CUH
#define PNB_KERNELS_CUH

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
/*
Alg of sec 3.3 of New features of latin dances.
*/
__global__ void compute_neutrality_kernel(unsigned long long seed, uint32_t enc_rounds, uint32_t dec_rounds,
int IDW, int IDB, int ODW, int ODB, 
int ntest, unsigned long long *d_result, int alg_type, 
int r1, int r2, int r3);


void compute_correlation_of_g(pnb_t *pnb);

void compute_neutrality_vector(pnb_t *pnb, uint64_t number_of_trials);

void get_pnb_list(pnb_t *pnb);

//compute complexity of pnb attack
void compute_complexity_of_the_attack(pnb_t *pnb);

void pnb_attack_for_single_bit_differential(
    int idw, 
    int idb, 
    int odw, 
    int odb, 
    int subrounds,
    int differential_part_subrounds,
    int linear_part_subrounds, 
    double threshold,
    uint64_t number_of_trials_for_neutrality,
    uint64_t number_of_trials_for_bias_of_g,
    int alg_type,
    FILE *output_file
    );

#endif