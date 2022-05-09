#ifndef FORRO_CUH
#define FORRO_CUH

#include "arx_cryptanalysis.cuh"

#define FORRO_NUMBER_OF_ROUNDS 12
#define FORRO_NUMBER_OF_SUBROUNDS_IN_EACH_ROUND 4
#define FORRO_TOTAL_OF_DIFFERENT_SUBROUNDS 8

__host__ __device__ void forro_define_rotations(int r1, int r2, int r3);

__host__ __device__ void forro_init(uint32_t state[STATE_SIZE], uint32_t k[KEY_SIZE], uint32_t nonce[NONCE_SIZE], uint32_t ctr[CTR_SIZE]);

__host__ __device__ void forro_odd_round(uint32_t x[STATE_SIZE]);

__host__ __device__ void forro_even_round(uint32_t x[STATE_SIZE]);

__host__ __device__ void forro_invert_odd_round(uint32_t x[STATE_SIZE]);

__host__ __device__ void forro_invert_even_round(uint32_t x[STATE_SIZE]);

__host__ __device__ void forro_rounds(uint32_t state[STATE_SIZE], uint32_t rounds, uint32_t last_round);

__host__ __device__ void forro_subrounds(uint32_t state[STATE_SIZE], uint32_t subrounds, uint32_t last_subround);

__host__ __device__ void forro_invert_subrounds(uint32_t state[STATE_SIZE], uint32_t subrounds, uint32_t last_subround);

__host__ __device__ void forro_invert_rounds(uint32_t state[STATE_SIZE], uint32_t rounds, uint32_t last_round);

__host__ __device__ void forro_encrypt_rounds(uint32_t final_state[STATE_SIZE], uint32_t initial_state[STATE_SIZE], uint32_t rounds);

__host__ __device__ void forro_encrypt_subrounds(uint32_t final_state[STATE_SIZE], uint32_t initial_state[STATE_SIZE], uint32_t subrounds);

__host__ __device__ void forro_decrypt_rounds(uint32_t final_state[STATE_SIZE], uint32_t initial_state[STATE_SIZE], 
    uint32_t intermediate_state[STATE_SIZE], uint32_t rounds, uint32_t last_round);

__host__ __device__ void forro_decrypt_subrounds(uint32_t final_state[STATE_SIZE], uint32_t initial_state[STATE_SIZE], 
    uint32_t intermediate_state[STATE_SIZE], uint32_t subrounds, uint32_t last_subround);

__host__ __device__ void forro_expand_bit(linear_approximation_t *L, int w, int bit, int expansionType);

__host__ __device__  void forro_special_cases(linear_approximation_t *L, int bit, int expansionType);

__host__ __device__ void forro_words_of_subround(int *W, int subround);

__host__ __device__ int forro_get_letter(int w, int subround);

__host__ __device__ int forro_get_position_of_letter(int letter, int subround);
#endif