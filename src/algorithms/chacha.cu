#include "arx_cryptanalysis.cuh"
#include "chacha.cuh"

#define HALF_1_CHACHAQUARTERROUND(a, b, c, d) \
a = PLUS(a, b);              \
d = ROTATE(XOR(d, a), 16);   \
c = PLUS(c, d);              \
b = ROTATE(XOR(b, c), 12);   \

#define HALF_2_CHACHAQUARTERROUND(a, b, c, d) \
a = PLUS(a, b);              \
d = ROTATE(XOR(d, a), 8);    \
c = PLUS(c, d);              \
b = ROTATE(XOR(b, c), 7);

#define INVERT_HALF_1_CHACHAQUARTERROUND(a, b, c, d) \
b = XOR(ROTATE(b,20), c); \
c = MINUS(c, d);              \
d = XOR(ROTATE(d,16), a); \
a = MINUS(a, b);

#define INVERT_HALF_2_CHACHAQUARTERROUND(a, b, c, d) \
b = XOR(ROTATE(b,25), c); \
c = MINUS(c, d);              \
d = XOR(ROTATE(d,24), a); \
a = MINUS(a, b);              \

#define CHACHAQUARTERROUND(a, b, c, d) \
HALF_1_CHACHAQUARTERROUND(a,b,c,d)\
HALF_2_CHACHAQUARTERROUND(a,b,c,d)\

#define INVERT_CHACHAQUARTERROUND(a,b,c,d)\
INVERT_HALF_2_CHACHAQUARTERROUND(a,b,c,d)\
INVERT_HALF_1_CHACHAQUARTERROUND(a,b,c,d)\

  __host__ __device__ void chacha_init(uint32_t state[MAXIMUM_STATE_SIZE], uint32_t k[KEY_SIZE], uint32_t nonce[NONCE_SIZE], uint32_t ctr[CTR_SIZE])
{
    state[0] = U32C(0x61707865);
    state[1] = U32C(0x3320646e);
    state[2] = U32C(0x79622d32);
    state[3] = U32C(0x6b206574);
    state[4] = k[0];
    state[5] = k[1];
    state[6] = k[2];
    state[7] = k[3];
    state[8] = k[4];
    state[9] = k[5];
    state[10] = k[6];
    state[11] = k[7];
    state[12] = ctr[0];
    state[13] = ctr[1];
    state[14] = nonce[0];
    state[15] = nonce[1];
}

  __host__ __device__ void chacha_odd_round(uint32_t x[MAXIMUM_STATE_SIZE])
{
    CHACHAQUARTERROUND(x[0], x[4], x[8], x[12])
    CHACHAQUARTERROUND(x[1], x[5], x[9], x[13])
    CHACHAQUARTERROUND(x[2], x[6], x[10], x[14])
    CHACHAQUARTERROUND(x[3], x[7], x[11], x[15])
}

  __host__ __device__ void chacha_even_round(uint32_t x[MAXIMUM_STATE_SIZE])
{
    CHACHAQUARTERROUND(x[0], x[5], x[10], x[15])
    CHACHAQUARTERROUND(x[1], x[6], x[11], x[12])
    CHACHAQUARTERROUND(x[2], x[7], x[8], x[13])
    CHACHAQUARTERROUND(x[3], x[4], x[9], x[14])
}

  __host__ __device__ void chacha_invert_odd_round(uint32_t x[MAXIMUM_STATE_SIZE])
{
    INVERT_CHACHAQUARTERROUND(x[3], x[7], x[11], x[15])
    INVERT_CHACHAQUARTERROUND(x[2], x[6], x[10], x[14])
    INVERT_CHACHAQUARTERROUND(x[1], x[5], x[9], x[13])
    INVERT_CHACHAQUARTERROUND(x[0], x[4], x[8], x[12])
}

  __host__ __device__ void chacha_invert_even_round(uint32_t x[MAXIMUM_STATE_SIZE])
{
    INVERT_CHACHAQUARTERROUND(x[3], x[4], x[9], x[14])
    INVERT_CHACHAQUARTERROUND(x[2], x[7], x[8], x[13])
    INVERT_CHACHAQUARTERROUND(x[1], x[6], x[11], x[12])
    INVERT_CHACHAQUARTERROUND(x[0], x[5], x[10], x[15])
}

  __host__ __device__ void chacha_rounds(uint32_t state[MAXIMUM_STATE_SIZE], uint32_t rounds, uint32_t last_round)
{
    uint32_t i;

    for (i = 1; i <= rounds; i++) {
        if ((i+last_round) % 2)
            chacha_odd_round(state);
        else
            chacha_even_round(state);
    }
}

  __host__ __device__ void chacha_invert_rounds(uint32_t state[MAXIMUM_STATE_SIZE], uint32_t rounds, uint32_t last_round)
{
    uint32_t i;

    for (i = 1; i <= rounds; i++) {
        if ((last_round + i) % 2)
            chacha_invert_even_round(state);
        else
            chacha_invert_odd_round(state);
    }
}



 __host__ __device__ void chacha_subrounds(uint32_t state[MAXIMUM_STATE_SIZE], uint32_t subrounds, uint32_t last_subround)
{
    uint32_t i, rounds;
    
    if((last_subround%CHACHA_TOTAL_OF_DIFFERENT_SUBROUNDS) != 0)
    {
        last_subround %= CHACHA_TOTAL_OF_DIFFERENT_SUBROUNDS;
        switch(last_subround)
        {
            case 1:    
                HALF_2_CHACHAQUARTERROUND(state[0], state[4], state[8], state[12])
                HALF_2_CHACHAQUARTERROUND(state[1], state[5], state[9], state[13])
                HALF_2_CHACHAQUARTERROUND(state[2], state[6], state[10], state[14])
                HALF_2_CHACHAQUARTERROUND(state[3], state[7], state[11], state[15])
                subrounds--;
                if(subrounds == 0) return;
            case 2:
                HALF_1_CHACHAQUARTERROUND(state[0], state[5], state[10], state[15])
                HALF_1_CHACHAQUARTERROUND(state[1], state[6], state[11], state[12])
                HALF_1_CHACHAQUARTERROUND(state[2], state[7], state[8], state[13])
                HALF_1_CHACHAQUARTERROUND(state[3], state[4], state[9], state[14])
                subrounds--;
                if(subrounds == 0) return;
            case 3:
                HALF_2_CHACHAQUARTERROUND(state[0], state[5], state[10], state[15])
                HALF_2_CHACHAQUARTERROUND(state[1], state[6], state[11], state[12])
                HALF_2_CHACHAQUARTERROUND(state[2], state[7], state[8], state[13])
                HALF_2_CHACHAQUARTERROUND(state[3], state[4], state[9], state[14])
                subrounds--;
                if(subrounds == 0) return;
        }
    }
    
    rounds = subrounds/CHACHA_TOTAL_OF_DIFFERENT_SUBROUNDS;
    subrounds %=CHACHA_TOTAL_OF_DIFFERENT_SUBROUNDS;

    for (i = 0; i < rounds; i++) {
        CHACHAQUARTERROUND(state[0], state[4], state[8], state[12])
        CHACHAQUARTERROUND(state[1], state[5], state[9], state[13])
        CHACHAQUARTERROUND(state[2], state[6], state[10], state[14])
        CHACHAQUARTERROUND(state[3], state[7], state[11], state[15])

        CHACHAQUARTERROUND(state[0], state[5], state[10], state[15])
        CHACHAQUARTERROUND(state[1], state[6], state[11], state[12])
        CHACHAQUARTERROUND(state[2], state[7], state[8], state[13])
        CHACHAQUARTERROUND(state[3], state[4], state[9], state[14])
    }

    if(subrounds==0)return;
    HALF_1_CHACHAQUARTERROUND(state[0], state[4], state[8], state[12])
    HALF_1_CHACHAQUARTERROUND(state[1], state[5], state[9], state[13])
    HALF_1_CHACHAQUARTERROUND(state[2], state[6], state[10], state[14])
    HALF_1_CHACHAQUARTERROUND(state[3], state[7], state[11], state[15])
    if(subrounds==1)return;
    HALF_2_CHACHAQUARTERROUND(state[0], state[4], state[8], state[12])
    HALF_2_CHACHAQUARTERROUND(state[1], state[5], state[9], state[13])
    HALF_2_CHACHAQUARTERROUND(state[2], state[6], state[10], state[14])
    HALF_2_CHACHAQUARTERROUND(state[3], state[7], state[11], state[15])
    if(subrounds==2)return;
    HALF_1_CHACHAQUARTERROUND(state[0], state[5], state[10], state[15])
    HALF_1_CHACHAQUARTERROUND(state[1], state[6], state[11], state[12])
    HALF_1_CHACHAQUARTERROUND(state[2], state[7], state[8], state[13])
    HALF_1_CHACHAQUARTERROUND(state[3], state[4], state[9], state[14])
}

 __host__ __device__ void chacha_invert_subrounds(uint32_t state[MAXIMUM_STATE_SIZE], uint32_t subrounds, uint32_t last_subround)
{
    uint32_t i, rounds;
    
    if((last_subround%CHACHA_TOTAL_OF_DIFFERENT_SUBROUNDS) != 0)
    {
        last_subround %= CHACHA_TOTAL_OF_DIFFERENT_SUBROUNDS;
        switch(last_subround)
        {
            case 3:    
                INVERT_HALF_1_CHACHAQUARTERROUND(state[0], state[5], state[10], state[15])
                INVERT_HALF_1_CHACHAQUARTERROUND(state[1], state[6], state[11], state[12])
                INVERT_HALF_1_CHACHAQUARTERROUND(state[2], state[7], state[8], state[13])
                INVERT_HALF_1_CHACHAQUARTERROUND(state[3], state[4], state[9], state[14])
                subrounds--;
                if(subrounds == 0) return;
            case 2:
                INVERT_HALF_2_CHACHAQUARTERROUND(state[0], state[4], state[8], state[12])
                INVERT_HALF_2_CHACHAQUARTERROUND(state[1], state[5], state[9], state[13])
                INVERT_HALF_2_CHACHAQUARTERROUND(state[2], state[6], state[10], state[14])
                INVERT_HALF_2_CHACHAQUARTERROUND(state[3], state[7], state[11], state[15])
                subrounds--;
                if(subrounds == 0) return;
            case 1:
                INVERT_HALF_1_CHACHAQUARTERROUND(state[0], state[4], state[8], state[12])
                INVERT_HALF_1_CHACHAQUARTERROUND(state[1], state[5], state[9], state[13])
                INVERT_HALF_1_CHACHAQUARTERROUND(state[2], state[6], state[10], state[14])
                INVERT_HALF_1_CHACHAQUARTERROUND(state[3], state[7], state[11], state[15])
                subrounds--;
                if(subrounds == 0) return;
        }
    }
    
    rounds = subrounds/CHACHA_TOTAL_OF_DIFFERENT_SUBROUNDS;
    subrounds %=CHACHA_TOTAL_OF_DIFFERENT_SUBROUNDS;

    for (i = 0; i < rounds; i++) {   
        INVERT_CHACHAQUARTERROUND(state[3], state[4], state[9], state[14])
        INVERT_CHACHAQUARTERROUND(state[2], state[7], state[8], state[13])
        INVERT_CHACHAQUARTERROUND(state[1], state[6], state[11], state[12])
        INVERT_CHACHAQUARTERROUND(state[0], state[5], state[10], state[15])
        INVERT_CHACHAQUARTERROUND(state[3], state[7], state[11], state[15])
        INVERT_CHACHAQUARTERROUND(state[2], state[6], state[10], state[14])
        INVERT_CHACHAQUARTERROUND(state[1], state[5], state[9], state[13])
        INVERT_CHACHAQUARTERROUND(state[0], state[4], state[8], state[12])
    }

    if(subrounds==0)return;
    INVERT_HALF_2_CHACHAQUARTERROUND(state[0], state[5], state[10], state[15])
    INVERT_HALF_2_CHACHAQUARTERROUND(state[1], state[6], state[11], state[12])
    INVERT_HALF_2_CHACHAQUARTERROUND(state[2], state[7], state[8], state[13])
    INVERT_HALF_2_CHACHAQUARTERROUND(state[3], state[4], state[9], state[14])
    if(subrounds==1)return;
    INVERT_HALF_1_CHACHAQUARTERROUND(state[0], state[5], state[10], state[15])
    INVERT_HALF_1_CHACHAQUARTERROUND(state[1], state[6], state[11], state[12])
    INVERT_HALF_1_CHACHAQUARTERROUND(state[2], state[7], state[8], state[13])
    INVERT_HALF_1_CHACHAQUARTERROUND(state[3], state[4], state[9], state[14])
    if(subrounds==2)return;
    INVERT_HALF_2_CHACHAQUARTERROUND(state[0], state[4], state[8], state[12])
    INVERT_HALF_2_CHACHAQUARTERROUND(state[1], state[5], state[9], state[13])
    INVERT_HALF_2_CHACHAQUARTERROUND(state[2], state[6], state[10], state[14])
    INVERT_HALF_2_CHACHAQUARTERROUND(state[3], state[7], state[11], state[15])
}

  __host__ __device__ void chacha_encrypt_rounds(uint32_t output[MAXIMUM_STATE_SIZE], uint32_t input[MAXIMUM_STATE_SIZE], uint32_t rounds)
{
    uint32_t x[MAXIMUM_STATE_SIZE];
    uint32_t i;

    for (i = 0; i < MAXIMUM_STATE_SIZE; ++i) x[i] = input[i];
    chacha_rounds(x, rounds, 0);
    for (i = 0; i < MAXIMUM_STATE_SIZE; ++i) x[i] = PLUS(x[i], input[i]);

    memcpy(output, x, MAXIMUM_STATE_SIZE_IN_BYTES);
}

  __host__ __device__ void chacha_decrypt_rounds(uint32_t output[MAXIMUM_STATE_SIZE], uint32_t input[MAXIMUM_STATE_SIZE], 
    uint32_t intermediate[MAXIMUM_STATE_SIZE], uint32_t rounds, uint32_t last_round)
{
    for (int i = 0; i < MAXIMUM_STATE_SIZE; ++i) intermediate[i] = MINUS(output[i], input[i]);
    chacha_invert_rounds(intermediate, rounds, last_round);
}

  __host__ __device__ void chacha_encrypt_subrounds(uint32_t output[MAXIMUM_STATE_SIZE], uint32_t input[MAXIMUM_STATE_SIZE], uint32_t subrounds)
{
    uint32_t x[MAXIMUM_STATE_SIZE];
    uint32_t i;

    for (i = 0; i < MAXIMUM_STATE_SIZE; ++i) x[i] = input[i];
    chacha_subrounds(x, subrounds, 0);
    for (i = 0; i < MAXIMUM_STATE_SIZE; ++i) x[i] = PLUS(x[i], input[i]);

    memcpy(output, x, MAXIMUM_STATE_SIZE_IN_BYTES);
}

  __host__ __device__ void chacha_decrypt_subrounds(uint32_t output[MAXIMUM_STATE_SIZE], uint32_t input[MAXIMUM_STATE_SIZE], 
    uint32_t intermediate[MAXIMUM_STATE_SIZE], uint32_t subrounds, uint32_t last_subround)
{
    for (int i = 0; i < MAXIMUM_STATE_SIZE; ++i) intermediate[i] = MINUS(output[i], input[i]);
    chacha_invert_subrounds(intermediate, subrounds, last_subround);
}


 __host__ __device__ void chacha_get_words_of_subround_for_target_word(int words[4], int subround, int target_word)
{
    int pattern[4][4][4] = { //[subrounds][number of words used] 
        {
            {0,4,8,12},
            {1,5,9,13},
            {2,6,10,14},
            {3,7,11,15},
        },
        {
            {0,4,8,12},
            {1,5,9,13},
            {2,6,10,14},
            {3,7,11,15},
        },
        {
            {0,5,10,15},
            {1,6,11,12},
            {2,7,8,13},
            {3,4,9,14}
        },
        {
            {0,5,10,15},
            {1,6,11,12},
            {2,7,8,13},
            {3,4,9,14}
        }
    };
    
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            if(pattern[(subround-1)%CHACHA_TOTAL_OF_DIFFERENT_SUBROUNDS][i][j] == target_word)
            {
                memcpy(words, pattern[(subround-1)%CHACHA_TOTAL_OF_DIFFERENT_SUBROUNDS][i], 4*sizeof(int));
                return;
            }
        }
    }
}

 __host__ __device__ void chacha_get_rotations_of_subround(int r[2], int subround)
{
    if(subround%2)
    {
        r[0] = 16;
        r[1] = 12;
    }
    else
    {
        r[0] = 8;
        r[1] = 7;
    }
}

 __host__ __device__ int chacha_get_letter(int target_word, int subround)
{
    int words[4];
    
    chacha_get_words_of_subround_for_target_word(words, subround, target_word);
    
    for(int i=0;i<4;i++)
        if(target_word == words[i])
            return i;
    
    return -1;
}

 __host__ __device__ void chacha_get_expansion(expansion e[4], int subround, int target_word)
{
    int words[4];
    int rotations[2];
    
    chacha_get_words_of_subround_for_target_word(words, subround, target_word);
    chacha_get_rotations_of_subround(rotations, subround);

    expansion exp[4] = {
        {
            {(uint32_t) words[LetterA],(uint32_t)  words[LetterB],(uint32_t)  words[LetterC],(uint32_t)  words[LetterB],(uint32_t)  words[LetterC] },
            {0,(uint32_t) rotations[1], 0,(uint32_t)  (rotations[1]-1), 31}, 5, 3, {0}
        },
        {
            {(uint32_t) words[LetterB], (uint32_t) words[LetterC]},
            {(uint32_t) rotations[1],0}, 2, 2, {0}
        },
        {
            {(uint32_t) words[LetterC],(uint32_t)  words[LetterD],(uint32_t)  words[LetterD]},
            {0,0,31}, 3, 2, {0}
        },
        {
            {(uint32_t) words[LetterA],(uint32_t)  words[LetterD]},
            {0,(uint32_t) rotations[0]}, 2, 2, {0}
        }
    };
    
    memcpy(e, exp, sizeof(expansion)*4);
}

 __host__ __device__ void chacha_expand_bit(linear_approximation_t *L, int w, int bit, int expansion_type)
{
    expansion exp[4];
    int subround = L->input.subround + 1, size_linear, letter;

    if(bit==0)
        expansion_type = EXPANSION_LINEAR;

    chacha_get_expansion(exp, subround, w);
    letter = chacha_get_letter(w, subround);

    if(expansion_type == EXPANSION_LINEAR)
        size_linear = exp[letter].number_of_linear_terms;
    else
        size_linear = exp[letter].number_of_terms;
    
    for (int j = 0; j < size_linear; j++)
        exp[letter].list_of_bits[j] = (exp[letter].list_of_bits[j] + 32 + bit) % 32;
    set_list_of_bits(L->output.mask, exp[letter].list_of_words, exp[letter].list_of_bits, size_linear);
}


#define HALF_1_CHACHAQUARTERROUND_DIFFERENTIAL_UPDATE(a, b, c, d, correlation_exponent) \
a = minHwMaxGammaDPA(a, b, correlation_exponent);\
d = ROTATE(XOR(d, a), 16);   \
c = minHwMaxGammaDPA(c, d, correlation_exponent);              \
b = ROTATE(XOR(b, c), 12);   \

#define HALF_2_CHACHAQUARTERROUND_DIFFERENTIAL_UPDATE(a, b, c, d, correlation_exponent) \
a = minHwMaxGammaDPA(a, b,correlation_exponent);              \
d = ROTATE(XOR(d, a), 8);    \
c = minHwMaxGammaDPA(c, d,correlation_exponent);              \
b = ROTATE(XOR(b, c), 7);

__host__ __device__ void chacha_differential_update(
    uint32_t diff[MAXIMUM_STATE_SIZE], 
    int subrounds,
    int *correlation_exponent
    )
{
    *correlation_exponent = 0;

    for(int r=0;r<subrounds;r++)
    {
        switch (r%4)
        {
        case 0:
            HALF_1_CHACHAQUARTERROUND_DIFFERENTIAL_UPDATE(diff[0], diff[4], diff[8], diff[12],correlation_exponent)
            HALF_1_CHACHAQUARTERROUND_DIFFERENTIAL_UPDATE(diff[1], diff[5], diff[9], diff[13],correlation_exponent)
            HALF_1_CHACHAQUARTERROUND_DIFFERENTIAL_UPDATE(diff[2], diff[6], diff[10], diff[14],correlation_exponent)
            HALF_1_CHACHAQUARTERROUND_DIFFERENTIAL_UPDATE(diff[3], diff[7], diff[11], diff[15],correlation_exponent)
            break;
        case 1:
            HALF_2_CHACHAQUARTERROUND_DIFFERENTIAL_UPDATE(diff[0], diff[4], diff[8], diff[12],correlation_exponent)
            HALF_2_CHACHAQUARTERROUND_DIFFERENTIAL_UPDATE(diff[1], diff[5], diff[9], diff[13],correlation_exponent)
            HALF_2_CHACHAQUARTERROUND_DIFFERENTIAL_UPDATE(diff[2], diff[6], diff[10], diff[14],correlation_exponent)
            HALF_2_CHACHAQUARTERROUND_DIFFERENTIAL_UPDATE(diff[3], diff[7], diff[11], diff[15],correlation_exponent)
            break;
        case 2:
            HALF_1_CHACHAQUARTERROUND_DIFFERENTIAL_UPDATE(diff[0], diff[5], diff[10], diff[15],correlation_exponent)
            HALF_1_CHACHAQUARTERROUND_DIFFERENTIAL_UPDATE(diff[1], diff[6], diff[11], diff[12],correlation_exponent)
            HALF_1_CHACHAQUARTERROUND_DIFFERENTIAL_UPDATE(diff[2], diff[7], diff[8], diff[13],correlation_exponent)
            HALF_1_CHACHAQUARTERROUND_DIFFERENTIAL_UPDATE(diff[3], diff[4], diff[9], diff[14],correlation_exponent)
            break;
        case 3:
            HALF_2_CHACHAQUARTERROUND_DIFFERENTIAL_UPDATE(diff[0], diff[5], diff[10], diff[15],correlation_exponent)
            HALF_2_CHACHAQUARTERROUND_DIFFERENTIAL_UPDATE(diff[1], diff[6], diff[11], diff[12],correlation_exponent)
            HALF_2_CHACHAQUARTERROUND_DIFFERENTIAL_UPDATE(diff[2], diff[7], diff[8], diff[13],correlation_exponent)
            HALF_2_CHACHAQUARTERROUND_DIFFERENTIAL_UPDATE(diff[3], diff[4], diff[9], diff[14],correlation_exponent)
            break;
        default:
            break;
        }
    }
}