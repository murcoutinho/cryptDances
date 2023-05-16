#ifndef UTIL_H
#define UTIL_H

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <string.h>
#include <curand.h>
#include <curand_kernel.h>
#include <mpi.h>
#include "types.cuh"
#include "arx_cryptanalysis.cuh"
#include "cryptanalysis_types.cuh"

#define XOR3(x,y,z) ((x) ^ (y) ^ (z))
#define EQ(x,y,z) ((~(x) ^ (y)) & (~(x) ^ (z)))

#define GENERATE_KEY_NONCE_CTR(k, nonce, ctr, rng)\
    k[0] = curand(&rng);\
    k[1] = curand(&rng);\
    k[2] = curand(&rng);\
    k[3] = curand(&rng);\
    k[4] = curand(&rng);\
    k[5] = curand(&rng);\
    k[6] = curand(&rng);\
    k[7] = curand(&rng);\
    nonce[0] = curand(&rng);\
    nonce[1] = curand(&rng);\
    ctr[0] = curand(&rng);\
    ctr[1] = curand(&rng);\

#define GENERATE_RANDOM_STATE(state,size)               \
    do {                                           \
        for(int i = 0; i < size; i++) {              \
            state[i] = curand(&rng);                \
        }                                          \
    } while(0)

#define is_empty(...) ( sizeof( (char[]){#__VA_ARGS__} ) == 1 )

#define dprintf(fp, ...)\
    {\
            printf (__VA_ARGS__);\
            if(fp!=NULL) fprintf (fp, __VA_ARGS__);\
    }\

__host__ __device__ void get_difference_from_single_bit(uint32_t diff[MAXIMUM_STATE_SIZE], int word, int bit);
__host__ __device__ void print_state(uint32_t state[MAXIMUM_STATE_SIZE]);
__host__ __device__ void xor_array(uint32_t *z, uint32_t *x, uint32_t *y, int size);
__host__ __device__ void and_array(uint32_t *z, uint32_t *x, uint32_t *y, int size);
__host__ __device__ uint8_t xor_bits_of_state(uint32_t state[MAXIMUM_STATE_SIZE], int size);
__host__ __device__ void transform_state_to_bits(uint32_t state[MAXIMUM_STATE_SIZE], uint8_t bits[MAXIMUM_STATE_SIZE]);
__host__ __device__ void update_result(int result[MAXIMUM_STATE_SIZE_IN_BITS], uint8_t bits[MAXIMUM_STATE_SIZE*8]);
__host__ __device__ void update_biases(double bias[MAXIMUM_STATE_SIZE_IN_BITS], uint32_t result[MAXIMUM_STATE_SIZE_IN_BITS], uint64_t N);
__host__ __device__ uint8_t check_parity_of_equation(uint32_t state[MAXIMUM_STATE_SIZE], uint32_t ODmask[MAXIMUM_STATE_SIZE], int size);
__host__ __device__ uint8_t check_parity_of_linear_relation(uint32_t inputMask[MAXIMUM_STATE_SIZE], 
    uint32_t inputState[MAXIMUM_STATE_SIZE], uint32_t outputMask[MAXIMUM_STATE_SIZE], uint32_t outputState[MAXIMUM_STATE_SIZE], int size);
__device__ int cudaCmp(uint8_t *v1, uint8_t *v2, int len);

__host__ __device__ void set_bit(uint32_t state[MAXIMUM_STATE_SIZE], uint32_t w, uint32_t bit);
__host__ __device__  void set_list_of_bits(uint32_t state[MAXIMUM_STATE_SIZE], uint32_t *w, uint32_t *bit, uint32_t numberOfBits);
__host__ __device__ int hamming_weight_state(uint32_t state[MAXIMUM_STATE_SIZE]);
__host__ __device__ uint32_t aop(uint32_t x);
__host__ __device__  uint32_t minHwMaxGammaDPA(uint32_t alpha, uint32_t beta, int *k);

int test_significance_of_correlation(double correlation, uint64_t number_of_trials);
void srand_by_rank();
uint64_t seed_by_rank(int iteration);
__host__ __device__  uint8_t get_bit_from_word_and_bit(uint32_t state[MAXIMUM_STATE_SIZE], uint32_t w, uint32_t bit);
int write_single_bit_differentials_to_file(
    const char *file_name, 
    differential_t diff[MAXIMUM_NUMBER_OF_POSSIBLE_SINGLE_BIT_DIFFERENTIALS],
    int alg_type
    );
int read_single_bit_differentials_from_file(
    const char *file_name, 
    differential_t diff[MAXIMUM_NUMBER_OF_POSSIBLE_SINGLE_BIT_DIFFERENTIALS],
    int alg_type
    );

int update_single_bit_differentials_from_file(
    const char *file_name, 
    differential_t diff[MAXIMUM_NUMBER_OF_POSSIBLE_SINGLE_BIT_DIFFERENTIALS],
    int alg_type
    );
void random_uint64(uint64_t *x);

void create_folder_if_doesnt_exist(const char *name);
#endif
