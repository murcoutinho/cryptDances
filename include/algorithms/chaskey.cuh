#ifndef CHASKEY_CUH
#define CHASKEY_CUH

#include "arx_cryptanalysis.cuh"

#define CHASKEY_STATE_SIZE 4
#define CHASKEY_KEY_SIZE 4
#define CHASKEY_NUMBER_OF_ROUNDS 8
#define CHASKEY_NUMBER_OF_SUBROUNDS_IN_EACH_ROUND 2

 __host__ __device__ void chaskey_init(uint32_t *state, uint32_t *k,uint32_t *nonce,uint32_t *ctr);
 __host__ __device__ void chaskey_round(uint32_t state[CHASKEY_STATE_SIZE]);
 __host__ __device__ void chaskey_subround(uint32_t state[CHASKEY_STATE_SIZE], uint32_t subround_number);
 __host__ __device__ void chaskey_rounds(uint32_t state[CHASKEY_STATE_SIZE], uint32_t rounds, uint32_t last_round);
 __host__ __device__ void chaskey_subrounds(uint32_t state[CHASKEY_STATE_SIZE], uint32_t subrounds, uint32_t last_subround);
 __host__ __device__ void chaskey_permute(uint32_t state[CHASKEY_STATE_SIZE]);

 #endif