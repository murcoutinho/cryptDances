#ifndef SALSA_CUH
#define SALSA_CUH

#include "arx_cryptanalysis.cuh"

#define SALSA_NUMBER_OF_ROUNDS 20
#define SALSA_NUMBER_OF_SUBROUNDS_IN_EACH_ROUND 1
#define SALSA_TOTAL_OF_DIFFERENT_SUBROUNDS 2

__host__ __device__  void salsa_init(uint32_t state[16], uint32_t k[8], uint32_t nonce[2], uint32_t ctr[2]);

__host__ __device__  void salsa_rounds(uint32_t state[16], uint32_t rounds, uint32_t last_round);

__host__ __device__  void salsa_invert_rounds(uint32_t state[16], uint32_t rounds, uint32_t last_round);

__host__ __device__  void salsa_subrounds(uint32_t state[STATE_SIZE], uint32_t subrounds, uint32_t last_subround);

__host__ __device__  void salsa_invert_subrounds(uint32_t state[STATE_SIZE], uint32_t subrounds, uint32_t last_subround);

__host__ __device__  void salsa_encrypt_rounds(uint32_t final_state[STATE_SIZE], uint32_t initial_state[STATE_SIZE], uint32_t rounds);

__host__ __device__  void salsa_encrypt_subrounds(uint32_t final_state[STATE_SIZE], uint32_t initial_state[STATE_SIZE], uint32_t subrounds);

__host__ __device__  void salsa_decrypt_rounds(uint32_t final_state[STATE_SIZE], uint32_t initial_state[STATE_SIZE], 
    uint32_t intermediate_state[STATE_SIZE], uint32_t rounds, uint32_t last_round);

__host__ __device__  void salsa_decrypt_subrounds(uint32_t final_state[STATE_SIZE], uint32_t initial_state[STATE_SIZE], 
    uint32_t intermediate_state[STATE_SIZE], uint32_t subrounds, uint32_t last_subround);

__host__ __device__ void salsa_differential_update(
    uint32_t diff[STATE_SIZE], 
    int subrounds,
    int *correlation_exponent
    );

#endif