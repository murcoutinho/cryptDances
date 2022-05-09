#include "arx_cryptanalysis.cuh"
#include "tests.cuh"

#define DETECT_ERROR(rv, d_rv, test_name, alg_name)\
    {\
        int de_flag = 0;\
        cudaMemcpy(rv, d_rv, sizeof(int)*TOTAL_THREADS_TESTS, cudaMemcpyDeviceToHost);\
        for (int de_i = 0; de_i < TOTAL_THREADS_TESTS; de_i++)\
            if (rv[de_i])\
                de_flag = 1;\
        if(de_flag){\
            printf("Results for test %s (alg = %s) at process %d: ", test_name, alg_name, my_rank);\
            printf("error!\n");\
        }\
    }\

__device__ void get_values_for_test_vector(uint32_t k[KEY_SIZE], uint32_t nonce[NONCE_SIZE], 
    uint32_t ctr[CTR_SIZE], uint8_t encrypt_input[STATE_SIZE_IN_BYTES], uint8_t test_vector[STATE_SIZE_IN_BYTES], int alg_type)
{
    uint8_t forro_test1_out[STATE_SIZE_IN_BYTES] =
    {
    0x5b,0x00,0xcd,0xfa,0xf0,0x51,0x25,0x54,0x12,0x9d,0xe2,0xe6,0xe7,0x85,0x81,0x7d,
    0x84,0xe6,0xb2,0x45,0x98,0xe0,0x37,0x66,0x99,0x5a,0x58,0xeb,0x3f,0x2d,0xd3,0x0f,
    0x05,0x1e,0xb1,0x88,0xbd,0x63,0x6b,0x7d,0x79,0x3b,0xda,0xe2,0x90,0x19,0xe6,0xf2,
    0x20,0xbd,0x5c,0x1a,0xec,0x7b,0x64,0x8e,0xa7,0x38,0x91,0xd3,0x82,0x3f,0xd4,0x5c
    };

    uint8_t chacha_test1_out[STATE_SIZE_IN_BYTES] = {0x76,0xb8,0xe0,0xad,0xa0,0xf1,0x3d,0x90,0x40,
    0x5d,0x6a,0xe5,0x53,0x86,0xbd,0x28,0xbd,0xd2,0x19,0xb8,0xa0,0x8d,0xed,0x1a,0xa8,0x36,0xef,
    0xcc,0x8b,0x77,0x0d,0xc7,0xda,0x41,0x59,0x7c,0x51,0x57,0x48,0x8d,0x77,0x24,0xe0,0x3f,0xb8,
    0xd8,0x4a,0x37,0x6a,0x43,0xb8,0xf4,0x15,0x18,0xa1,0x1c,0xc3,0x87,0xb6,0x69,0xb2,0xee,0x65,0x86};

    uint8_t salsa_test_in[STATE_SIZE_IN_BYTES] = {211,159, 13,115, 76, 55, 82,183, 3,117,222, 37,191,
    187,234,136, 49,237,179, 48, 1,106,178,219,175,199,166, 48, 86, 16,179,207, 31,240, 32, 63, 15, 
    83, 93,161,116,147, 48,113,238, 55,204, 36, 79,201,235, 79, 3, 81,156, 47,203, 26,244,243, 88,118,104, 54};

    uint8_t salsa_test_out[STATE_SIZE_IN_BYTES] = {109, 42,178,168,156,240,248,238,168,196,190,203, 
    26,110,170,154, 29, 29,150, 26,150, 30,235,249,190,163,251, 48, 69,144, 51, 57, 118, 40,152,157,
    180, 57, 27, 94,107, 42,236, 35, 27,111,114,114, 219,236,232,135,111,155,110, 18, 24,232, 95,158,179, 19, 48,202};

    const char *key8 = "minha vida e andar por este pais";
    const char *iv8 = "mostro a";
    const char *count8 = " alegria";

    switch (alg_type)
    {
    case ALG_TYPE_CHACHA:
        for(int i=0;i<KEY_SIZE;i++)
            k[i] = 0;
        for(int i=0;i<CTR_SIZE;i++)
            ctr[i] = 0;
        for(int i=0;i<NONCE_SIZE;i++)
            nonce[i] = 0;
        for(int i=0;i<STATE_SIZE_IN_BYTES;i++)
            test_vector[i] = chacha_test1_out[i];
        break;
    case ALG_TYPE_FORRO:
        for(int i=0;i<KEY_SIZE;i++)
            k[i] = *((uint32_t *)(key8 + 4*i));
        for(int i=0;i<CTR_SIZE;i++)
            ctr[i] = *((uint32_t *)(count8 + 4*i));
        for(int i=0;i<NONCE_SIZE;i++)
            nonce[i] = *((uint32_t *)(iv8 + 4*i));
        for(int i=0;i<STATE_SIZE_IN_BYTES;i++)
            test_vector[i] = forro_test1_out[i];
        break;
    case ALG_TYPE_SALSA:
        for(int i=0;i<STATE_SIZE_IN_BYTES;i++)
            encrypt_input[i] = salsa_test_in[i];
        for(int i=0;i<STATE_SIZE_IN_BYTES;i++)
            test_vector[i] = salsa_test_out[i];
        break;
    default:
        break;
    }
}

/*
* Tests the implementation of Forro, ChaCha and Salsa with test vectors.
*/
__global__ void ker_test_vectors(int *rv, int alg_type)
{
    uint32_t k[KEY_SIZE] = { 0 }, ctr[CTR_SIZE] = { 0 }, nonce[NONCE_SIZE] = { 0 }, 
        state[STATE_SIZE] = { 0 }, state_final[STATE_SIZE];
    uint8_t encrypt_input[STATE_SIZE_IN_BYTES], test_vector[STATE_SIZE_IN_BYTES];
    algorithm alg;
    int tx = threadIdx.x + blockIdx.x * blockDim.x;

    define_alg(&alg, alg_type);
    rv[tx] = RV_ERROR;

    get_values_for_test_vector(k, nonce, ctr, encrypt_input, test_vector, alg_type);

    if(alg_type == ALG_TYPE_SALSA)
        memcpy(state, encrypt_input, STATE_SIZE_IN_BYTES);
    else
        alg.init(state, k, nonce, ctr);

    alg.encrypt_rounds(state_final, state, alg.number_of_rounds);

    if (cudaCmp((uint8_t *)state_final, test_vector, STATE_SIZE_IN_BYTES) == 0)
        rv[tx] = RV_SUCESS;
}

__global__ void ker_test_round_vs_invert_round(int *rv, int alg_type, curandState *devStates)
{
    uint32_t k[KEY_SIZE] = { 0 }, ctr[CTR_SIZE] = { 0 }, nonce[NONCE_SIZE] = { 0 }, 
        state[STATE_SIZE] = { 0 }, state_alt[STATE_SIZE] ={0};
    algorithm alg;
    int tx = threadIdx.x + blockIdx.x * blockDim.x;
    int forward=0, backward=1;

    curandState rng = devStates[tx];

    GENERATE_KEY_NONCE_CTR(k,nonce,ctr,rng);
    define_alg(&alg, alg_type);

    rv[tx] = RV_SUCESS;

    //Randomize the number of rounds (forward and backward)
    while(forward<backward){
        forward = (curand(&rng))%20 + 1;
        backward = (curand(&rng))%20 + 1;
    }

    alg.init(state, k, nonce, ctr);
    alg.init(state_alt, k, nonce, ctr);

    alg.rounds(state, forward, 0);
    alg.rounds(state_alt, forward-backward, 0);

    alg.invert_rounds(state, backward, forward);

    if (cudaCmp((uint8_t *)state, (uint8_t *)state_alt, STATE_SIZE_IN_BYTES))
        rv[tx] = RV_ERROR;
}


__global__ void ker_test_subround_vs_invert_subround(int *rv, int alg_type, curandState *devStates)
{
   uint32_t k[KEY_SIZE] = { 0 }, ctr[CTR_SIZE] = { 0 }, nonce[NONCE_SIZE] = { 0 }, 
        state[STATE_SIZE] = { 0 }, state_alt[STATE_SIZE] ={0};
    algorithm alg;
    int tx = threadIdx.x + blockIdx.x * blockDim.x;
    int forward=0, backward=1;

    curandState rng = devStates[tx];

    GENERATE_KEY_NONCE_CTR(k,nonce,ctr,rng);
    define_alg(&alg, alg_type);

    rv[tx] = RV_SUCESS;

    //Randomize the number of rounds (forward and backward)
    while(forward<backward){
        forward = (curand(&rng))%40 + 1;
        backward = (curand(&rng))%40 + 1;
    }

    alg.init(state, k, nonce, ctr);
    alg.init(state_alt, k, nonce, ctr);

    alg.subrounds(state, forward, 0);
    alg.subrounds(state_alt, forward-backward, 0);

    alg.invert_subrounds(state, backward, forward);

    if (cudaCmp((uint8_t *)state, (uint8_t *)state_alt, STATE_SIZE_IN_BYTES))
        rv[tx] = RV_ERROR;
}


__global__ void ker_test_round_vs_subround(int *rv, int alg_type, curandState *devStates)
{
    uint32_t k[KEY_SIZE] = { 0 }, ctr[CTR_SIZE] = { 0 }, nonce[NONCE_SIZE] = { 0 }, 
        state[STATE_SIZE] = { 0 }, state_alt[STATE_SIZE] ={0};
    algorithm alg;
    int tx = threadIdx.x + blockIdx.x * blockDim.x;
    int rounds;

    curandState rng = devStates[tx];

    GENERATE_KEY_NONCE_CTR(k,nonce,ctr,rng);
    define_alg(&alg, alg_type);

    rv[tx] = RV_SUCESS;

    //Randomize the number of rounds
    rounds = (curand(&rng))%20 + 1;

    alg.init(state, k, nonce, ctr);
    alg.init(state_alt, k, nonce, ctr);

    alg.rounds(state, rounds, 0);
    alg.subrounds(state_alt, rounds * alg.number_of_subrounds_in_one_round, 0);

    if (cudaCmp((uint8_t *)state, (uint8_t *)state_alt, STATE_SIZE_IN_BYTES))
        rv[tx] = RV_ERROR;
}

__global__ void ker_test_encrypt_decrypt(int *rv, int alg_type, curandState *devStates)
{
    uint32_t k[KEY_SIZE] = { 0 }, ctr[CTR_SIZE] = { 0 }, nonce[NONCE_SIZE] = { 0 }, 
        state[STATE_SIZE] = { 0 }, state_final[STATE_SIZE], inverted[STATE_SIZE],
        state_final_subrounds[STATE_SIZE], inverted_subrounds[STATE_SIZE];
    algorithm alg;
    int tx = threadIdx.x + blockIdx.x * blockDim.x;

    curandState rng = devStates[tx];

    GENERATE_KEY_NONCE_CTR(k,nonce,ctr,rng);
    define_alg(&alg, alg_type);

    rv[tx] = RV_SUCESS;

    alg.init(state, k, nonce, ctr);
    alg.encrypt_rounds(state_final, state, alg.number_of_rounds);
    alg.encrypt_subrounds(state_final_subrounds, state, alg.number_of_rounds * alg.number_of_subrounds_in_one_round);
    if (cudaCmp((uint8_t *)state_final_subrounds, (uint8_t *)state_final, STATE_SIZE_IN_BYTES))
        rv[tx] = RV_ERROR;

    alg.decrypt_rounds(state_final, state, inverted, alg.number_of_rounds, alg.number_of_rounds);
    if (cudaCmp((uint8_t *)inverted, (uint8_t *)state, STATE_SIZE_IN_BYTES))
        rv[tx] = RV_ERROR;

    alg.decrypt_subrounds(state_final_subrounds, state, inverted_subrounds, 
        alg.number_of_rounds * alg.number_of_subrounds_in_one_round, alg.number_of_rounds * alg.number_of_subrounds_in_one_round);
    if (cudaCmp((uint8_t *)inverted_subrounds, (uint8_t *)state, STATE_SIZE_IN_BYTES))
        rv[tx] = RV_ERROR;
}


__global__ void ker_test_partial_decryption(int *rv, int alg_type, curandState *devStates)
{
    uint32_t k[KEY_SIZE] = { 0 }, ctr[CTR_SIZE] = { 0 }, nonce[NONCE_SIZE] = { 0 }, 
        state[STATE_SIZE] = { 0 }, intermediate_state[STATE_SIZE] = { 0 }, state_final[STATE_SIZE], inverted[STATE_SIZE];
    algorithm alg;
    int tx = threadIdx.x + blockIdx.x * blockDim.x, subrounds,subrounds_dec;

    curandState rng = devStates[tx];

    GENERATE_KEY_NONCE_CTR(k,nonce,ctr,rng);
    define_alg(&alg, alg_type);

    rv[tx] = RV_SUCESS;

    //Randomize the number of rounds
    subrounds = (curand(&rng))%20 + 1;

    subrounds_dec = (curand(&rng))%20 + 1;
    while(subrounds_dec>subrounds)
        subrounds_dec--;

    alg.init(state, k, nonce, ctr);
    alg.encrypt_subrounds(state_final, state, subrounds);

    alg.init(intermediate_state, k, nonce, ctr);
    alg.subrounds(intermediate_state, subrounds - subrounds_dec, 0);

    alg.decrypt_subrounds(state_final, state, inverted, subrounds_dec, subrounds);

    if (cudaCmp((uint8_t *)inverted, (uint8_t *)intermediate_state, STATE_SIZE_IN_BYTES))
        rv[tx] = RV_ERROR;
}

int test_algorithms_on_gpu(curandState *dev_states, int alg_type)
{
    int *d_rvs;
    int results[TOTAL_THREADS_TESTS];
    char name[10];
    
    get_alg_name(name ,alg_type);
    
    if(my_rank == 0)
        return RV_SUCESS;

    cudaSetDevice((my_rank-1)%8);
    cudaMalloc((void **)&d_rvs, sizeof(int) * TOTAL_THREADS_TESTS);

    ker_test_vectors<<<NUMBER_OF_THREADS_TESTS, NUMBER_OF_BLOCKS_TESTS>>>(d_rvs, alg_type);
    DETECT_ERROR(results, d_rvs, "ker_test_vectors", name);

    ker_test_round_vs_invert_round<<<NUMBER_OF_THREADS_TESTS, NUMBER_OF_BLOCKS_TESTS>>>(d_rvs, alg_type, dev_states);
    DETECT_ERROR(results, d_rvs, "ker_test_round_vs_invert_round", name);

    ker_test_subround_vs_invert_subround<<<NUMBER_OF_THREADS_TESTS, NUMBER_OF_BLOCKS_TESTS>>>(d_rvs, alg_type, dev_states);
    DETECT_ERROR(results, d_rvs, "ker_test_subround_vs_invert_subround", name);

    ker_test_round_vs_subround<<<NUMBER_OF_THREADS_TESTS, NUMBER_OF_BLOCKS_TESTS>>>(d_rvs, alg_type, dev_states);
    DETECT_ERROR(results, d_rvs, "ker_test_round_vs_subround", name);

    ker_test_encrypt_decrypt<<<NUMBER_OF_THREADS_TESTS, NUMBER_OF_BLOCKS_TESTS>>>(d_rvs, alg_type, dev_states);
    DETECT_ERROR(results, d_rvs, "ker_test_encrypt_decrypt", name);

    ker_test_partial_decryption<<<NUMBER_OF_THREADS_TESTS, NUMBER_OF_BLOCKS_TESTS>>>(d_rvs, alg_type, dev_states);
    DETECT_ERROR(results, d_rvs, "ker_test_partial_decryption", name);

    cudaFree(d_rvs);

    return RV_SUCESS;
}