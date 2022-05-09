#include "arx_cryptanalysis.cuh"
#include "salsa.cuh"

__host__ __device__ void salsa_odd_round(uint32_t x[16])
{
    x[4] = XOR(x[4], ROTATE(PLUS(x[0], x[12]), 7));
    x[8] = XOR(x[8], ROTATE(PLUS(x[4], x[0]), 9));
    x[12] = XOR(x[12], ROTATE(PLUS(x[8], x[4]), 13));
    x[0] = XOR(x[0], ROTATE(PLUS(x[12], x[8]), 18));
    x[9] = XOR(x[9], ROTATE(PLUS(x[5], x[1]), 7));
    x[13] = XOR(x[13], ROTATE(PLUS(x[9], x[5]), 9));
    x[1] = XOR(x[1], ROTATE(PLUS(x[13], x[9]), 13));
    x[5] = XOR(x[5], ROTATE(PLUS(x[1], x[13]), 18));
    x[14] = XOR(x[14], ROTATE(PLUS(x[10], x[6]), 7));
    x[2] = XOR(x[2], ROTATE(PLUS(x[14], x[10]), 9));
    x[6] = XOR(x[6], ROTATE(PLUS(x[2], x[14]), 13));
    x[10] = XOR(x[10], ROTATE(PLUS(x[6], x[2]), 18));
    x[3] = XOR(x[3], ROTATE(PLUS(x[15], x[11]), 7));
    x[7] = XOR(x[7], ROTATE(PLUS(x[3], x[15]), 9));
    x[11] = XOR(x[11], ROTATE(PLUS(x[7], x[3]), 13));
    x[15] = XOR(x[15], ROTATE(PLUS(x[11], x[7]), 18));
}

__host__ __device__  void salsa_even_round(uint32_t x[16])
{
    x[1] = XOR(x[1], ROTATE(PLUS(x[0], x[3]), 7));
    x[2] = XOR(x[2], ROTATE(PLUS(x[1], x[0]), 9));
    x[3] = XOR(x[3], ROTATE(PLUS(x[2], x[1]), 13));
    x[0] = XOR(x[0], ROTATE(PLUS(x[3], x[2]), 18));
    x[6] = XOR(x[6], ROTATE(PLUS(x[5], x[4]), 7));
    x[7] = XOR(x[7], ROTATE(PLUS(x[6], x[5]), 9));
    x[4] = XOR(x[4], ROTATE(PLUS(x[7], x[6]), 13));
    x[5] = XOR(x[5], ROTATE(PLUS(x[4], x[7]), 18));
    x[11] = XOR(x[11], ROTATE(PLUS(x[10], x[9]), 7));
    x[8] = XOR(x[8], ROTATE(PLUS(x[11], x[10]), 9));
    x[9] = XOR(x[9], ROTATE(PLUS(x[8], x[11]), 13));
    x[10] = XOR(x[10], ROTATE(PLUS(x[9], x[8]), 18));
    x[12] = XOR(x[12], ROTATE(PLUS(x[15], x[14]), 7));
    x[13] = XOR(x[13], ROTATE(PLUS(x[12], x[15]), 9));
    x[14] = XOR(x[14], ROTATE(PLUS(x[13], x[12]), 13));
    x[15] = XOR(x[15], ROTATE(PLUS(x[14], x[13]), 18));
}

__host__ __device__  void salsa_invert_odd_round(uint32_t x[16])
{
    x[15] = XOR(x[15], ROTATE(PLUS(x[11], x[7]), 18));
    x[11] = XOR(x[11], ROTATE(PLUS(x[7], x[3]), 13));
    x[7] = XOR(x[7], ROTATE(PLUS(x[3], x[15]), 9));
    x[3] = XOR(x[3], ROTATE(PLUS(x[15], x[11]), 7));
    x[10] = XOR(x[10], ROTATE(PLUS(x[6], x[2]), 18));
    x[6] = XOR(x[6], ROTATE(PLUS(x[2], x[14]), 13));
    x[2] = XOR(x[2], ROTATE(PLUS(x[14], x[10]), 9));
    x[14] = XOR(x[14], ROTATE(PLUS(x[10], x[6]), 7));
    x[5] = XOR(x[5], ROTATE(PLUS(x[1], x[13]), 18));
    x[1] = XOR(x[1], ROTATE(PLUS(x[13], x[9]), 13));
    x[13] = XOR(x[13], ROTATE(PLUS(x[9], x[5]), 9));
    x[9] = XOR(x[9], ROTATE(PLUS(x[5], x[1]), 7));
    x[0] = XOR(x[0], ROTATE(PLUS(x[12], x[8]), 18));
    x[12] = XOR(x[12], ROTATE(PLUS(x[8], x[4]), 13));
    x[8] = XOR(x[8], ROTATE(PLUS(x[4], x[0]), 9));
    x[4] = XOR(x[4], ROTATE(PLUS(x[0], x[12]), 7));
}

__host__ __device__  void salsa_invert_even_round(uint32_t x[16])
{
    x[15] = XOR(x[15], ROTATE(PLUS(x[14], x[13]), 18));
    x[14] = XOR(x[14], ROTATE(PLUS(x[13], x[12]), 13));
    x[13] = XOR(x[13], ROTATE(PLUS(x[12], x[15]), 9));
    x[12] = XOR(x[12], ROTATE(PLUS(x[15], x[14]), 7));
    x[10] = XOR(x[10], ROTATE(PLUS(x[9], x[8]), 18));
    x[9] = XOR(x[9], ROTATE(PLUS(x[8], x[11]), 13));
    x[8] = XOR(x[8], ROTATE(PLUS(x[11], x[10]), 9));
    x[11] = XOR(x[11], ROTATE(PLUS(x[10], x[9]), 7));
    x[5] = XOR(x[5], ROTATE(PLUS(x[4], x[7]), 18));
    x[4] = XOR(x[4], ROTATE(PLUS(x[7], x[6]), 13));
    x[7] = XOR(x[7], ROTATE(PLUS(x[6], x[5]), 9));
    x[6] = XOR(x[6], ROTATE(PLUS(x[5], x[4]), 7));
    x[0] = XOR(x[0], ROTATE(PLUS(x[3], x[2]), 18));
    x[3] = XOR(x[3], ROTATE(PLUS(x[2], x[1]), 13));
    x[2] = XOR(x[2], ROTATE(PLUS(x[1], x[0]), 9));
    x[1] = XOR(x[1], ROTATE(PLUS(x[0], x[3]), 7));
}

__host__ __device__  void salsa_init(uint32_t state[16], uint32_t k[8], uint32_t nonce[2], uint32_t ctr[2])
{
    const char sigma[17] = "expand 32-byte k";
    const char *constants;
    constants = sigma;

    state[1] = k[0];
    state[2] = k[1];
    state[3] = k[2];
    state[4] = k[3];
    state[11] = k[4];
    state[12] = k[5];
    state[13] = k[6];
    state[14] = k[7];
#if 1
    state[0] = U8TO32_LITTLE(constants + 0);
    state[5] = U8TO32_LITTLE(constants + 4);
    state[10] = U8TO32_LITTLE(constants + 8);
    state[15] = U8TO32_LITTLE(constants + 12);
#endif
    state[6] = nonce[0];
    state[7] = nonce[1];
    state[8] = ctr[0];
    state[9] = ctr[1];
}


__host__ __device__  void salsa_rounds(uint32_t state[16], uint32_t rounds, uint32_t last_round)
{
    uint32_t i;

    for (i = 1; i <= rounds; i++) {
        if ((i + last_round) % 2)
            salsa_odd_round(state);
        else
            salsa_even_round(state);
    }
}

__host__ __device__  void salsa_invert_rounds(uint32_t state[16], uint32_t rounds, uint32_t last_round)
{
    uint32_t i;

    for (i = 1; i <= rounds; i++) {
        if ((last_round + i) % 2)
            salsa_invert_even_round(state);
        else
            salsa_invert_odd_round(state);
    }
}

__host__ __device__  void salsa_subrounds(uint32_t state[STATE_SIZE], uint32_t subrounds, uint32_t last_subround)
{
    salsa_rounds(state, subrounds, last_subround);
}

__host__ __device__  void salsa_invert_subrounds(uint32_t state[STATE_SIZE], uint32_t subrounds, uint32_t last_subround)
{
    salsa_invert_rounds(state, subrounds, last_subround);
}

__host__ __device__  void salsa_encrypt_rounds(uint32_t final_state[STATE_SIZE], uint32_t initial_state[STATE_SIZE], uint32_t rounds)
{
    uint32_t x[16];
    uint32_t i;
    for (i = 0; i < 16; ++i) x[i] = initial_state[i];
    salsa_rounds(x, rounds, 0);
    for (i = 0; i < 16; ++i) final_state[i] = PLUS(x[i], initial_state[i]);
}

__host__ __device__  void salsa_encrypt_subrounds(uint32_t final_state[STATE_SIZE], uint32_t initial_state[STATE_SIZE], uint32_t subrounds)
{
    salsa_encrypt_rounds(final_state, initial_state, subrounds);
}

__host__ __device__  void salsa_decrypt_rounds(uint32_t final_state[STATE_SIZE], uint32_t initial_state[STATE_SIZE], 
    uint32_t intermediate_state[STATE_SIZE], uint32_t rounds, uint32_t last_round)
{
    for (int i = 0; i < 16; ++i) intermediate_state[i] = MINUS(final_state[i], initial_state[i]);
    salsa_invert_rounds(intermediate_state, rounds, last_round);
}

__host__ __device__  void salsa_decrypt_subrounds(uint32_t final_state[STATE_SIZE], uint32_t initial_state[STATE_SIZE], 
    uint32_t intermediate_state[STATE_SIZE], uint32_t subrounds, uint32_t last_subround)
{
    salsa_decrypt_rounds(final_state, initial_state, intermediate_state, subrounds, last_subround);
}

__host__ __device__ void salsa_differential_update_odd_round(
    uint32_t diff[STATE_SIZE], 
    int *correlation_exponent
    )
{
    uint32_t aux;

    aux = minHwMaxGammaDPA(diff[0], diff[12], correlation_exponent);
    diff[4] = XOR(diff[4], ROTATE(aux, 7));
    aux = minHwMaxGammaDPA(diff[4], diff[0], correlation_exponent);
    diff[8] = XOR(diff[8], ROTATE(aux, 9));
    aux = minHwMaxGammaDPA(diff[8], diff[4], correlation_exponent);
    diff[12] = XOR(diff[12], ROTATE(aux, 13));
    aux = minHwMaxGammaDPA(diff[12], diff[8], correlation_exponent);
    diff[0] = XOR(diff[0], ROTATE(aux, 18));
    aux = minHwMaxGammaDPA(diff[5], diff[1], correlation_exponent);
    diff[9] = XOR(diff[9], ROTATE(aux, 7));
    aux = minHwMaxGammaDPA(diff[9], diff[5],correlation_exponent);
    diff[13] = XOR(diff[13], ROTATE(aux, 9));
    aux = minHwMaxGammaDPA(diff[13], diff[9], correlation_exponent);
    diff[1] = XOR(diff[1], ROTATE(aux, 13));
    aux = minHwMaxGammaDPA(diff[1], diff[13], correlation_exponent);
    diff[5] = XOR(diff[5], ROTATE(aux, 18));
    aux = minHwMaxGammaDPA(diff[10], diff[6], correlation_exponent);
    diff[14] = XOR(diff[14], ROTATE(aux, 7));
    aux = minHwMaxGammaDPA(diff[14], diff[10], correlation_exponent);
    diff[2] = XOR(diff[2], ROTATE(aux, 9));
    aux = minHwMaxGammaDPA(diff[2], diff[14], correlation_exponent);
    diff[6] = XOR(diff[6], ROTATE(aux, 13));
    aux = minHwMaxGammaDPA(diff[6], diff[2], correlation_exponent);
    diff[10] = XOR(diff[10], ROTATE(aux, 18));
    aux = minHwMaxGammaDPA(diff[15], diff[11], correlation_exponent);
    diff[3] = XOR(diff[3], ROTATE(aux, 7));
    aux = minHwMaxGammaDPA(diff[3], diff[15], correlation_exponent);
    diff[7] = XOR(diff[7], ROTATE(aux, 9));
    aux = minHwMaxGammaDPA(diff[7], diff[3], correlation_exponent);
    diff[11] = XOR(diff[11], ROTATE(aux, 13));
    aux = minHwMaxGammaDPA(diff[11], diff[7], correlation_exponent);
    diff[15] = XOR(diff[15], ROTATE(aux, 18));
}


__host__ __device__ void salsa_differential_update_even_round(
    uint32_t diff[STATE_SIZE], 
    int *correlation_exponent
    )
{
    uint32_t aux;
    
    aux = minHwMaxGammaDPA(diff[0], diff[3], correlation_exponent);
    diff[1] = XOR(diff[1], ROTATE(aux, 7));
    aux = minHwMaxGammaDPA(diff[1], diff[0], correlation_exponent);
    diff[2] = XOR(diff[2], ROTATE(aux, 9));
    aux = minHwMaxGammaDPA(diff[2], diff[1], correlation_exponent);
    diff[3] = XOR(diff[3], ROTATE(aux, 13));
    aux = minHwMaxGammaDPA(diff[3], diff[2], correlation_exponent);
    diff[0] = XOR(diff[0], ROTATE(aux, 18));
    aux = minHwMaxGammaDPA(diff[5], diff[4], correlation_exponent);
    diff[6] = XOR(diff[6], ROTATE(aux, 7));
    aux = minHwMaxGammaDPA(diff[6], diff[5], correlation_exponent);
    diff[7] = XOR(diff[7], ROTATE(aux, 9));
    aux = minHwMaxGammaDPA(diff[7], diff[6], correlation_exponent);
    diff[4] = XOR(diff[4], ROTATE(aux, 13));
    aux = minHwMaxGammaDPA(diff[4], diff[7], correlation_exponent);
    diff[5] = XOR(diff[5], ROTATE(aux, 18));
    aux = minHwMaxGammaDPA(diff[10], diff[9],correlation_exponent);
    diff[11] = XOR(diff[11], ROTATE(aux, 7));
    aux = minHwMaxGammaDPA(diff[11], diff[10],correlation_exponent);
    diff[8] = XOR(diff[8], ROTATE(aux, 9));
    aux = minHwMaxGammaDPA(diff[8], diff[11], correlation_exponent);
    diff[9] = XOR(diff[9], ROTATE(aux, 13));
    aux = minHwMaxGammaDPA(diff[9], diff[8], correlation_exponent);
    diff[10] = XOR(diff[10], ROTATE(aux, 18));
    aux = minHwMaxGammaDPA(diff[15], diff[14], correlation_exponent);
    diff[12] = XOR(diff[12], ROTATE(aux, 7));
    aux = minHwMaxGammaDPA(diff[12], diff[15], correlation_exponent);
    diff[13] = XOR(diff[13], ROTATE(aux, 9));
    aux = minHwMaxGammaDPA(diff[13], diff[12], correlation_exponent);
    diff[14] = XOR(diff[14], ROTATE(aux, 13));
    aux = minHwMaxGammaDPA(diff[14], diff[13], correlation_exponent);
    diff[15] = XOR(diff[15], ROTATE(aux, 18));	
}

__host__ __device__ void salsa_differential_update(
    uint32_t diff[STATE_SIZE], 
    int subrounds,
    int *correlation_exponent
    )
{
    *correlation_exponent = 0;

    for (int i = 1; i <= subrounds; i++) {
        if (i % 2)
            salsa_differential_update_odd_round(diff, correlation_exponent);
        else
            salsa_differential_update_even_round(diff, correlation_exponent);
    }
}