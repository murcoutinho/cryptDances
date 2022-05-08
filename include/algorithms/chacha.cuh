#ifndef CHACHA_CUH
#define CHACHA_CUH

#include "arx_cryptanalysis.cuh"

#define CHACHA_NUMBER_OF_ROUNDS 20
#define CHACHA_NUMBER_OF_SUBROUNDS_IN_EACH_ROUND 2
#define CHACHA_TOTAL_OF_DIFFERENT_SUBROUNDS 4

 __host__ __device__ void chacha_init(uint32_t state[STATE_SIZE], uint32_t k[KEY_SIZE], uint32_t nonce[NONCE_SIZE], uint32_t ctr[CTR_SIZE]);
 __host__ __device__ void chacha_odd_round(uint32_t x[STATE_SIZE]);
 __host__ __device__ void chacha_even_round(uint32_t x[STATE_SIZE]);
 __host__ __device__ void chacha_invert_odd_round(uint32_t x[STATE_SIZE]);
 __host__ __device__ void chacha_invert_even_round(uint32_t x[STATE_SIZE]);
 __host__ __device__ void chacha_rounds(uint32_t state[STATE_SIZE], uint32_t rounds, uint32_t last_round);
 __host__ __device__ void chacha_invert_rounds(uint32_t state[STATE_SIZE], uint32_t rounds, uint32_t last_round);
__host__ __device__ void chacha_subrounds(uint32_t state[STATE_SIZE], uint32_t subrounds, uint32_t last_subround);
__host__ __device__ void chacha_invert_subrounds(uint32_t state[STATE_SIZE], uint32_t subrounds, uint32_t last_subround);

 __host__ __device__ void chacha_encrypt_rounds(uint32_t output[STATE_SIZE], uint32_t input[STATE_SIZE], uint32_t rounds);

 __host__ __device__ void chacha_decrypt_rounds(uint32_t output[STATE_SIZE], uint32_t input[STATE_SIZE], 
	uint32_t intermediate[STATE_SIZE], uint32_t rounds, uint32_t last_round);

 __host__ __device__ void chacha_encrypt_subrounds(uint32_t output[STATE_SIZE], uint32_t input[STATE_SIZE], uint32_t subrounds);

 __host__ __device__ void chacha_decrypt_subrounds(uint32_t output[STATE_SIZE], uint32_t input[STATE_SIZE], 
	uint32_t intermediate[STATE_SIZE], uint32_t subrounds, uint32_t last_subround);

 __host__ __device__ void chacha_expand_bit(linear_approximation_t *L, int w, int bit, int expansion_type);

 __host__ __device__ int chacha_get_letter(int target_word, int subround);

__host__ __device__ void chacha_differential_update(
	uint32_t diff[STATE_SIZE], 
	int subrounds,
	int *correlation_exponent
	);
 #endif