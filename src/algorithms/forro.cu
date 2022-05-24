#include "arx_cryptanalysis.cuh"
#include "forro.cuh"

#define QUARTERROUND(a, b, c, d, e) \
d = PLUS(d, e); \
c = XOR(c, d); \
b = PLUS(b, c); \
b = ROTATE(b,R1);\
a = PLUS(a, b); \
e = XOR(e, a); \
d = PLUS(d, e); \
d = ROTATE(d,R2);\
c = PLUS(c, d); \
b = XOR(b, c); \
a = PLUS(a, b); \
a = ROTATE(a,R3);

#define INVERT_QUARTERROUND(a,b,c,d,e)\
a = ROTATE(a,(32-R3));\
a = MINUS(a, b); \
b = XOR(b, c); \
c = MINUS(c, d); \
d = ROTATE(d,(32-R2));\
d = MINUS(d, e); \
e = XOR(e, a); \
a = MINUS(a, b); \
b = ROTATE(b,(32-R1));\
b = MINUS(b, c); \
c = XOR(c, d); \
d = MINUS(d, e);

#define LOAD32_LE(v) (*((uint32_t *) (v)))
#define STORE32_LE(c,x) (memcpy(c,&x,4))

__host__ __device__ void forro_init(uint32_t state[STATE_SIZE], uint32_t k[KEY_SIZE], uint32_t nonce[NONCE_SIZE], uint32_t ctr[CTR_SIZE])
{
    const char *constants = "voltadaasabranca";
    state[0] = k[0];
    state[1] = k[1];
    state[2] = k[2];
    state[3] = k[3];
    state[4] = ctr[0];
    state[5] = ctr[1];
    state[6] = U8TO32_LITTLE(constants + 0);
    state[7] = U8TO32_LITTLE(constants + 4);
    state[8] = k[4];
    state[9] = k[5];
    state[10] = k[6];
    state[11] = k[7];
    state[12] = nonce[0];
    state[13] = nonce[1];
    state[14] = U8TO32_LITTLE(constants + 8);
    state[15] = U8TO32_LITTLE(constants + 12);
}

__host__ __device__ void forro_odd_round(uint32_t x[STATE_SIZE])
{
    QUARTERROUND(x[0], x[4], x[8], x[12], x[3])
    QUARTERROUND(x[1], x[5], x[9], x[13], x[0])
    QUARTERROUND(x[2], x[6], x[10], x[14], x[1])
    QUARTERROUND(x[3], x[7], x[11], x[15], x[2])
}

__host__ __device__ void forro_even_round(uint32_t x[STATE_SIZE])
{
    QUARTERROUND(x[0], x[5], x[10], x[15], x[3])
    QUARTERROUND(x[1], x[6], x[11], x[12], x[0])
    QUARTERROUND(x[2], x[7], x[8], x[13], x[1])
    QUARTERROUND(x[3], x[4], x[9], x[14], x[2])
}

__host__ __device__ void forro_invert_odd_round(uint32_t x[STATE_SIZE])
{
    INVERT_QUARTERROUND(x[3], x[7], x[11], x[15], x[2])
    INVERT_QUARTERROUND(x[2], x[6], x[10], x[14], x[1])
    INVERT_QUARTERROUND(x[1], x[5], x[9], x[13], x[0])
    INVERT_QUARTERROUND(x[0], x[4], x[8], x[12], x[3])
}

__host__ __device__ void forro_invert_even_round(uint32_t x[STATE_SIZE])
{
    INVERT_QUARTERROUND(x[3], x[4], x[9], x[14], x[2])
    INVERT_QUARTERROUND(x[2], x[7], x[8], x[13], x[1])
    INVERT_QUARTERROUND(x[1], x[6], x[11], x[12], x[0])
    INVERT_QUARTERROUND(x[0], x[5], x[10], x[15], x[3])
}

__host__ __device__ void forro_rounds(uint32_t state[STATE_SIZE], uint32_t rounds, uint32_t last_round)
{
    uint32_t i;

    for (i = 1; i <= rounds; i++) {
        if ((last_round + i) % 2)
            forro_odd_round(state);
        else
            forro_even_round(state);
    }
}

__host__ __device__ void forro_invert_rounds(uint32_t state[STATE_SIZE], uint32_t rounds, uint32_t last_round)
{
    uint32_t i;

    for (i = 1; i <= rounds; i++) {
        if ((last_round + i) % 2)
            forro_invert_even_round(state);
        else
            forro_invert_odd_round(state);
    }
}

__host__ __device__ void forro_subrounds(uint32_t state[STATE_SIZE], uint32_t subrounds, uint32_t last_subround)
{
    uint32_t i, rounds;
    
    if((last_subround%FORRO_TOTAL_OF_DIFFERENT_SUBROUNDS) != 0)
    {
        last_subround %= FORRO_TOTAL_OF_DIFFERENT_SUBROUNDS;
        switch(last_subround)
        {
            case 1:
                QUARTERROUND(state[1], state[5], state[9], state[13], state[0]);
                subrounds--;
                if(subrounds == 0) return;
            case 2:
                QUARTERROUND(state[2], state[6], state[10], state[14], state[1]);
                subrounds--;
                if(subrounds == 0) return;
            case 3:
                QUARTERROUND(state[3], state[7], state[11], state[15], state[2]);
                subrounds--;
                if(subrounds == 0) return;
            case 4:
                QUARTERROUND(state[0], state[5], state[10], state[15], state[3]);
                subrounds--;
                if(subrounds == 0) return;
            case 5:
                QUARTERROUND(state[1], state[6], state[11], state[12], state[0]);
                subrounds--;
                if(subrounds == 0) return;
            case 6:
                QUARTERROUND(state[2], state[7], state[8], state[13], state[1]);
                subrounds--;
                if(subrounds == 0) return;
            case 7:
                QUARTERROUND(state[3], state[4], state[9], state[14], state[2]);
                subrounds--;
                if(subrounds == 0) return;
        }
    }
    
    rounds = subrounds/FORRO_TOTAL_OF_DIFFERENT_SUBROUNDS;
    subrounds %=FORRO_TOTAL_OF_DIFFERENT_SUBROUNDS;

    for (i = 0; i < rounds; i++) {
        QUARTERROUND(state[0], state[4], state[8], state[12], state[3])
        QUARTERROUND(state[1], state[5], state[9], state[13], state[0])
        QUARTERROUND(state[2], state[6], state[10], state[14], state[1])
        QUARTERROUND(state[3], state[7], state[11], state[15], state[2])
        QUARTERROUND(state[0], state[5], state[10], state[15], state[3])
        QUARTERROUND(state[1], state[6], state[11], state[12], state[0])
        QUARTERROUND(state[2], state[7], state[8], state[13], state[1])
        QUARTERROUND(state[3], state[4], state[9], state[14], state[2])
    }

    if(subrounds==0)return;
    QUARTERROUND(state[0], state[4], state[8], state[12], state[3])
    if(subrounds==1)return;
    QUARTERROUND(state[1], state[5], state[9], state[13], state[0])
    if(subrounds==2)return;
    QUARTERROUND(state[2], state[6], state[10], state[14], state[1])
    if(subrounds==3)return;
    QUARTERROUND(state[3], state[7], state[11], state[15], state[2])
    if(subrounds==4)return;
    QUARTERROUND(state[0], state[5], state[10], state[15], state[3])
    if(subrounds==5)return;
    QUARTERROUND(state[1], state[6], state[11], state[12], state[0])
    if(subrounds==6)return;
    QUARTERROUND(state[2], state[7], state[8], state[13], state[1])
}




__host__ __device__ void forro_invert_subrounds(uint32_t state[STATE_SIZE], uint32_t subrounds, uint32_t last_subround)
{
    uint32_t i, rounds;
    
    if((last_subround%FORRO_TOTAL_OF_DIFFERENT_SUBROUNDS) != 0)
    {
        last_subround %= FORRO_TOTAL_OF_DIFFERENT_SUBROUNDS;
        switch(last_subround)
        {
            case 7:
                INVERT_QUARTERROUND(state[2], state[7], state[8], state[13], state[1]);
                subrounds--;
                if(subrounds == 0) return;
            case 6:
                INVERT_QUARTERROUND(state[1], state[6], state[11], state[12], state[0]);
                subrounds--;
                if(subrounds == 0) return;
            case 5:
                INVERT_QUARTERROUND(state[0], state[5], state[10], state[15], state[3]);
                subrounds--;
                if(subrounds == 0) return;
            case 4:
                INVERT_QUARTERROUND(state[3], state[7], state[11], state[15], state[2]);
                subrounds--;
                if(subrounds == 0) return;
            case 3:
                INVERT_QUARTERROUND(state[2], state[6], state[10], state[14], state[1]);
                subrounds--;
                if(subrounds == 0) return;
            case 2:
                INVERT_QUARTERROUND(state[1], state[5], state[9], state[13], state[0]);
                subrounds--;
                if(subrounds == 0) return;
            case 1:
                INVERT_QUARTERROUND(state[0], state[4], state[8], state[12], state[3])
                subrounds--;
                if(subrounds == 0) return;
        }
    }
    
    rounds = subrounds/FORRO_TOTAL_OF_DIFFERENT_SUBROUNDS;
    subrounds %=FORRO_TOTAL_OF_DIFFERENT_SUBROUNDS;

    for (i = 0; i < rounds; i++) {
        INVERT_QUARTERROUND(state[3], state[4], state[9], state[14], state[2])
        INVERT_QUARTERROUND(state[2], state[7], state[8], state[13], state[1])
        INVERT_QUARTERROUND(state[1], state[6], state[11], state[12], state[0])
        INVERT_QUARTERROUND(state[0], state[5], state[10], state[15], state[3])
        INVERT_QUARTERROUND(state[3], state[7], state[11], state[15], state[2])
        INVERT_QUARTERROUND(state[2], state[6], state[10], state[14], state[1])
        INVERT_QUARTERROUND(state[1], state[5], state[9], state[13], state[0])
        INVERT_QUARTERROUND(state[0], state[4], state[8], state[12], state[3])
    }

    if(subrounds==0)return;
    INVERT_QUARTERROUND(state[3], state[4], state[9], state[14], state[2])
    if(subrounds==1)return;
    INVERT_QUARTERROUND(state[2], state[7], state[8], state[13], state[1])
    if(subrounds==2)return;
    INVERT_QUARTERROUND(state[1], state[6], state[11], state[12], state[0])
    if(subrounds==3)return;
    INVERT_QUARTERROUND(state[0], state[5], state[10], state[15], state[3])
    if(subrounds==4)return;
    INVERT_QUARTERROUND(state[3], state[7], state[11], state[15], state[2])
    if(subrounds==5)return;
    INVERT_QUARTERROUND(state[2], state[6], state[10], state[14], state[1])
    if(subrounds==6)return;
    INVERT_QUARTERROUND(state[1], state[5], state[9], state[13], state[0])
}


__host__ __device__ void forro_encrypt_rounds(uint32_t final_state[STATE_SIZE], uint32_t initial_state[STATE_SIZE], uint32_t rounds)
{
    uint32_t x[STATE_SIZE];
    uint32_t i;

    for (i = 0; i < STATE_SIZE; ++i) x[i] = initial_state[i];
    forro_rounds(x, rounds,0);
    for (i = 0; i < STATE_SIZE; ++i) x[i] = PLUS(x[i], initial_state[i]);

    for(i=0; i<STATE_SIZE;i++)
        U32TO8_LITTLE(final_state+i, x[i]);
}

__host__ __device__ void forro_encrypt_subrounds(uint32_t final_state[STATE_SIZE], uint32_t initial_state[STATE_SIZE], uint32_t subrounds)
{
    uint32_t x[STATE_SIZE];
    uint32_t i;

    for (i = 0; i < STATE_SIZE; ++i) x[i] = initial_state[i];
    forro_subrounds(x, subrounds, 0);
    for (i = 0; i < STATE_SIZE; ++i) x[i] = PLUS(x[i], initial_state[i]);

    for(i=0; i<STATE_SIZE;i++)
        U32TO8_LITTLE(final_state+i, x[i]);
}

__host__ __device__ void forro_decrypt_rounds(uint32_t final_state[STATE_SIZE], uint32_t initial_state[STATE_SIZE], 
    uint32_t intermediate_state[STATE_SIZE], uint32_t rounds, uint32_t last_round)
{
    for (int i = 0; i < STATE_SIZE; ++i) intermediate_state[i] = MINUS(final_state[i], initial_state[i]);
    forro_invert_rounds(intermediate_state, rounds, last_round);
}

__host__ __device__ void forro_decrypt_subrounds(uint32_t final_state[STATE_SIZE], uint32_t initial_state[STATE_SIZE], 
    uint32_t intermediate_state[STATE_SIZE], uint32_t subrounds, uint32_t last_subround)
{
    for (int i = 0; i < STATE_SIZE; ++i) intermediate_state[i] = MINUS(final_state[i], initial_state[i]);
    forro_invert_subrounds(intermediate_state, subrounds, last_subround);
}

__host__ __device__ void forro_words_of_subround(int W[5], int subround)
{
    uint32_t pattern[8][5] = {
        {0,4,8,12,3},
        {1,5,9,13,0},
        {2,6,10,14,1},
        {3,7,11,15,2},
        {0,5,10,15,3},
        {1,6,11,12,0},
        {2,7,8,13,1},
        {3,4,9,14,2}
    };
    
    memcpy(W, pattern[(subround-1)%8], 5*sizeof(int));
}

__host__ __device__ int forro_get_letter(int w, int subround)
{
    int W[5];
    
    forro_words_of_subround(W, subround);
    
    for(int i=0;i<5;i++)
        if(w == W[i])
            return i;
    
    return -1;
}

__host__ __device__ void forro_get_expansion(expansion e[5], int subround)
{
    int W[5];
    
    forro_words_of_subround(W, subround);
    expansion exp[5] = {
        {
            {(uint32_t) W[LetterA], (uint32_t) W[LetterC], (uint32_t) W[LetterC]},
            {8,0,31}, 3, 2, {0}
        },
        {
            {(uint32_t) W[LetterB], (uint32_t) W[LetterC], (uint32_t) W[LetterC], (uint32_t) W[LetterD], (uint32_t) W[LetterC], (uint32_t) W[LetterD]},
            {10, 10,0, 0, 31, 31}, 6, 4, {0}
        },
        {
            {(uint32_t) W[LetterC], (uint32_t) W[LetterD],(uint32_t)  W[LetterD],(uint32_t)  W[LetterE],(uint32_t)  W[LetterD],(uint32_t)  W[LetterE]},
            {0,0,27,0,31,31}, 6, 4, {0}
        },
        {
            {(uint32_t) W[LetterD],(uint32_t)  W[LetterA],(uint32_t)  W[LetterB],(uint32_t)  W[LetterA],(uint32_t)  W[LetterB]},
            {27,8,0,7,31}, 5, 3, {0}
        },
        {
            {(uint32_t) W[LetterE],(uint32_t)  W[LetterA],(uint32_t)  W[LetterB],(uint32_t)  W[LetterB]},
            {0,8,0,31}, 4, 3, { {SPECIAL_EXPANSION_EDD, 3} }
        }
    };
    
    memcpy(e, exp, sizeof(expansion)*5);
}

__host__ __device__ void forro_expand_bit(linear_approximation_t *L, int w, int bit, int expansionType)
{
    expansion exp[5];
    int subround = L->input.subround + 1, size, letter;

    forro_get_expansion(exp, subround);
    letter = forro_get_letter(w, subround);
    if(expansionType == EXPANSION_LINEAR)
        size = exp[letter].number_of_linear_terms;
    else
        size = exp[letter].number_of_terms;
    
    for (int j = 0; j < size; j++)
        exp[letter].list_of_bits[j] = (exp[letter].list_of_bits[j] + 32 + bit) % 32;
    set_list_of_bits(L->output.mask, exp[letter].list_of_words, exp[letter].list_of_bits, size);
}


__host__ __device__  void forro_special_cases(linear_approximation_t *L, int bit, int expansionType)
{
    expansion exp[5];
    int subround = L->input.subround + 1;
    forro_get_expansion(exp, subround);
    
    for(int i=0; i<5; i++)
    {
        for(int s=0; s<NUMBER_OF_SPECIAL_CASES; s++)
        {
            if(exp[i].special_cases_bits[s][0] == expansionType)
            {
                int position = exp[i].special_cases_bits[s][1];
                set_bit(L->output.mask, exp[i].list_of_words[position], (exp[i].list_of_bits[position] + 32 + bit) % 32);
            }
        }
    }
}