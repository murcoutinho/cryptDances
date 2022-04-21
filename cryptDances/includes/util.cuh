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

#define GENERATE_RANDOM_STATE(state)\
	state[0] = curand(&rng);\
	state[1] = curand(&rng);\
	state[2] = curand(&rng);\
	state[3] = curand(&rng);\
	state[4] = curand(&rng);\
	state[5] = curand(&rng);\
	state[6] = curand(&rng);\
	state[7] = curand(&rng);\
	state[8] = curand(&rng);\
	state[9] = curand(&rng);\
	state[10] = curand(&rng);\
	state[11] = curand(&rng);\
	state[12] = curand(&rng);\
	state[13] = curand(&rng);\
	state[14] = curand(&rng);\
	state[15] = curand(&rng);\


#define is_empty(...) ( sizeof( (char[]){#__VA_ARGS__} ) == 1 )

#define dprintf(fp, ...)\
	{\
			printf (__VA_ARGS__);\
			if(fp!=NULL) fprintf (fp, __VA_ARGS__);\
	}\

__host__ __device__ void get_difference_from_single_bit(uint32_t diff[STATE_SIZE], int word, int bit);
__host__ __device__ void print_state(uint32_t state[STATE_SIZE]);
__host__ __device__ void xor_array(uint32_t *z, uint32_t *x, uint32_t *y, int size);
__host__ __device__ void and_array(uint32_t *z, uint32_t *x, uint32_t *y, int size);
__host__ __device__ uint8_t xor_bits_of_state(uint32_t state[STATE_SIZE]);
__host__ __device__ void transform_state_to_bits(uint32_t state[STATE_SIZE], uint8_t bits[STATE_SIZE]);
__host__ __device__ void update_result(int result[STATE_SIZE_IN_BITS], uint8_t bits[STATE_SIZE*8]);
__host__ __device__ void update_biases(double bias[STATE_SIZE_IN_BITS], uint32_t result[STATE_SIZE_IN_BITS], uint64_t N);
__host__ __device__ uint8_t check_parity_of_equation(uint32_t state[STATE_SIZE], uint32_t ODmask[STATE_SIZE]);
__host__ __device__ uint8_t check_parity_of_linear_relation(uint32_t inputMask[STATE_SIZE], 
	uint32_t inputState[STATE_SIZE], uint32_t outputMask[STATE_SIZE], uint32_t outputState[STATE_SIZE]);
__device__ int cudaCmp(uint8_t *v1, uint8_t *v2, int len);

__host__ __device__ void set_bit(uint32_t state[STATE_SIZE], uint32_t w, uint32_t bit);
__host__ __device__  void set_list_of_bits(uint32_t state[STATE_SIZE], uint32_t *w, uint32_t *bit, uint32_t numberOfBits);
__host__ __device__ int hamming_weight_state(uint32_t state[STATE_SIZE]);
__host__ __device__ uint32_t aop(uint32_t x);
__host__ __device__  uint32_t minHwMaxGammaDPA(uint32_t alpha, uint32_t beta, int *k);

int test_significance_of_correlation(double correlation, uint64_t number_of_trials);
void srand_by_rank();
uint64_t seed_by_rank();
__host__ __device__  uint8_t get_bit_from_word_and_bit(uint32_t state[STATE_SIZE], uint32_t w, uint32_t bit);
int write_single_bit_differentials_to_file(
	const char *file_name, 
	differential_t diff[NUMBER_OF_POSSIBLE_SINGLE_BIT_DIFFERENTIALS]
	);
int read_single_bit_differentials_from_file(
	const char *file_name, 
	differential_t diff[NUMBER_OF_POSSIBLE_SINGLE_BIT_DIFFERENTIALS]
	);

int update_single_bit_differentials_from_file(
	const char *file_name, 
	differential_t diff[NUMBER_OF_POSSIBLE_SINGLE_BIT_DIFFERENTIALS]
	);
#endif