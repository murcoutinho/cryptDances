#include "arx_cryptanalysis.cuh"
#include <inttypes.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

__host__ __device__ uint32_t aop(uint32_t x)
{
    uint32_t X[6], Y[6];
    uint32_t L = 4;
    uint32_t one = 1;

    X[0] = x & (x >> 1);
    for (int i = 1; i < L; i++)
        X[i] = X[i - 1] & (X[i - 1] >> (one << (i)));

    Y[0] = x & (~X[0]);

    for (int i = 1; i < (L + 1); i++)
        Y[i] = Y[i - 1] | ((Y[i - 1] >> (one << (i))) & X[i - 1]);

    return Y[L];
}

//minimizes the hamming weight of the result of the sum
__host__ __device__  uint32_t minHwMaxGammaDPA(uint32_t alpha, uint32_t beta, int *k)
{
    uint32_t L = 32;
    uint32_t p = 0, a0, b0, a1, b1, g0, g1, gamma = 0, p1;

    a0 = alpha & 1;
    b0 = beta & 1;
    g0 = a0 ^ b0;

    p = aop(~(alpha^beta)) & ( (~(alpha^beta))>>1 & (alpha ^ (alpha>>1)) );

    for (int i = 1; i < L; i++)
    {
        a1 = (alpha >> i) & 1;
        b1 = (beta >> i) & 1;
        p1 = (p >> i) & 1;

        if (EQ(a0, b0, g0) & 1)
            g1 = a1 ^ b1 ^ a0;
        else if ((i == (L - 1)) || (a1 != b1) || (p1 == 1))
        {
            (*k)++;
            g1 = 0;
        }
        else
        {
            (*k)++;
            g1 = a1;
        }

        gamma |= (g1 << i);
        a0 = a1;
        b0 = b1;
        g0 = g1;
    }

    return gamma;
}

int test_significance_of_correlation(double correlation, uint64_t number_of_trials)
{
    double comp = (2 * (0.5 + 4 * sqrt(0.25 / number_of_trials)) - 1);
    if ((fabs(correlation) > comp))
        return TRUE;
    else
        return FALSE;
}

void srand_by_rank()
{
    srand((unsigned int) (my_rank*time(NULL)));
}

uint64_t seed_by_rank()
{
    //Guarantees that each MPI process has a different seed
    int shift = 1 + (int) log2((double)num_procs);
    return((rand()<<shift) ^ my_rank);
}

void random_uint64(uint64_t *x)
{
    *x = 0;
    *x |= (rand() ^ my_rank & 0xFFFF);
    *x |= (((uint64_t)rand()  ^ my_rank& 0xFFFF)) << 16;
    *x |= ((uint64_t)(rand()  ^ my_rank& 0xFFFF)) << 32;
    *x |= ((uint64_t)(rand()  ^ my_rank & 0xFFFF)) << 48;
}

__host__ __device__ void get_difference_from_single_bit(uint32_t diff[MAXIMUM_STATE_SIZE], int word, int bit)
{
    memset(diff, 0x00, sizeof(uint32_t) * MAXIMUM_STATE_SIZE);
    diff[word] = 1<<bit;
}


__host__ __device__ void print_state(uint32_t state[MAXIMUM_STATE_SIZE])
{
    for(int i=0;i<MAXIMUM_STATE_SIZE;i++)
    {
        printf("%02X ", state[i]);
        if((i%4) == 3)
        printf("\n");
    }
}



__host__ __device__ void xor_array(uint32_t *z, uint32_t *x, uint32_t *y, int size)
{
    int i;
    for (i = 0; i < size; i++)
        z[i] = x[i] ^ y[i];
}


__host__ __device__ void and_array(uint32_t *z, uint32_t *x, uint32_t *y, int size)
{
    for (int i = 0; i < size; i++)
        z[i] = x[i] & y[i];
}


__host__ __device__ uint8_t xor_bits_of_state(uint32_t state[MAXIMUM_STATE_SIZE], int size)
{
    uint32_t x = state[0];
    for (int i = 1; i < size; i++)
        x ^= state[i];

    x = x ^ (x >> 16);
    x = x ^ (x >> 8);
    x = x ^ (x >> 4);
    x = x ^ (x >> 2);
    return ((x ^ (x >> 1)) & 1);
}

__host__ __device__ void transform_state_to_bits(uint32_t state[MAXIMUM_STATE_SIZE], uint8_t bits[STATE_SIZE_IN_BITS])
{
    int i, b;
    int count = 0;
    for (i = 0; i < MAXIMUM_STATE_SIZE; i++)
    {
        for (b = 0; b < 32; b++)
        {
            bits[count] = (state[i] >> b) & 1;
            count++;
        }
    }
}

__host__ __device__ void update_result(int result[STATE_SIZE_IN_BITS], uint8_t bits[STATE_SIZE_IN_BITS])
{
    int i;
    for (i = 0; i < 512; i++)
        result[i] += (int)bits[i];
}

__host__ __device__ void update_biases(double bias[STATE_SIZE_IN_BITS], uint32_t result[STATE_SIZE_IN_BITS], uint64_t N)
{
    uint32_t i;
    for (i = 0; i < 512; i++)
        bias[i] += ((double)result[i]) / N;
}

__host__ __device__ uint8_t check_parity_of_equation(uint32_t state[MAXIMUM_STATE_SIZE], uint32_t ODmask[MAXIMUM_STATE_SIZE], int size)
{
    uint32_t aux[MAXIMUM_STATE_SIZE];

    and_array(aux, state, ODmask, size);
    return(xor_bits_of_state(aux, size));
}

__host__ __device__ uint8_t check_parity_of_linear_relation(uint32_t input_mask[MAXIMUM_STATE_SIZE], 
    uint32_t input_state[MAXIMUM_STATE_SIZE], uint32_t output_mask[MAXIMUM_STATE_SIZE], uint32_t output_state[MAXIMUM_STATE_SIZE], int size)
{
    uint32_t aux[MAXIMUM_STATE_SIZE];

    for (int i = 0; i < size; i++)
        aux[i] = (input_state[i] & input_mask[i]) ^ (output_state[i] & output_mask[i]);

    return(xor_bits_of_state(aux,size));
}


__device__ int cudaCmp(uint8_t *v1, uint8_t *v2, int len)
{
    for (int i = 0; i < len; i++)
        if (v1[i] != v2[i])
            return 1;

    return 0;
}

__host__ __device__ void set_bit(uint32_t state[MAXIMUM_STATE_SIZE], uint32_t w, uint32_t bit)
{
    state[w] ^= (1 << bit);
}

__host__ __device__ void set_list_of_bits(uint32_t state[MAXIMUM_STATE_SIZE], uint32_t *w, uint32_t *bit, uint32_t numberOfBits)
{
    for (uint32_t i = 0; i < numberOfBits; i++)
        set_bit(state, w[i], bit[i]);
}

__host__ __device__  uint8_t get_bit_from_word_and_bit(uint32_t state[MAXIMUM_STATE_SIZE], uint32_t w, uint32_t bit)
{
    if(bit>31) return 0;
    return((state[w] >> bit) & 1);
}

__host__ __device__ int hamming_weight_state(uint32_t state[MAXIMUM_STATE_SIZE], int size)
{
    int hw = 0;
    
    for (int w = 0; w < size; w++)
    {
        for (int b = 0; b < 32; b++)
        {
            hw += get_bit_from_word_and_bit(state, w,b);
        }
    }

    return hw;
}


int write_single_bit_differentials_to_file(
    const char *file_name, 
    differential_t diff[MAXIMUM_NUMBER_OF_POSSIBLE_SINGLE_BIT_DIFFERENTIALS],
    int alg_type
    )
{
    FILE *p = NULL;
    uint32_t state_size_in_bits = sizeof(uint32_t)*get_state_size(alg_type);
    uint32_t iv_size_in_bits = sizeof(uint32_t)*get_state_size(alg_type);
    uint32_t number_of_possible_single_bit_differentials = state_size_in_bits * iv_size_in_bits;
    
    //First check if input is a list of single bit differentials
    for(int i =0;i<number_of_possible_single_bit_differentials; i++)
    {
        if(diff[i].input.number_of_bits != 1)
            return RV_ERROR;
        if(diff[i].output.number_of_bits != 1)
            return RV_ERROR;
    }

    p = fopen(file_name, "w");
    if(p == NULL)
        return RV_ERROR;

    fprintf(p, "IDW IDB ODW ODB number_of_trials correlation_count correlation is_significant\n");

    for(int i =0;i<number_of_possible_single_bit_differentials; i++)
    {
        fprintf(p, "%d %d %d %d ", diff[i].input.words[0], diff[i].input.bits[0],diff[i].output.words[0], diff[i].output.bits[0]);
        fprintf(p, "%" PRIu64 " %" PRIu64 " ", diff[i].correlation.number_of_trials, diff[i].correlation.correlation_count);
        fprintf(p, "%.15f %d\n", diff[i].correlation, diff[i].correlation.is_significant);
    }

    fclose(p);
    return RV_SUCESS;
}


int read_single_bit_differentials_from_file(
    const char *file_name, 
    differential_t diff[MAXIMUM_NUMBER_OF_POSSIBLE_SINGLE_BIT_DIFFERENTIALS]
    )
{
    FILE *p = NULL;
    char header[200], *pheader = (char *)header;
    size_t header_size = 200;
    uint32_t state_size_in_bits = sizeof(uint32_t)*get_state_size(alg_type);
    uint32_t iv_size_in_bits = sizeof(uint32_t)*get_state_size(alg_type);
    uint32_t number_of_possible_single_bit_differentials = state_size_in_bits * iv_size_in_bits;

    p = fopen(file_name, "r");
    if(p == NULL)
        return RV_ERROR;

    getline(&pheader, &header_size, p);

    for(int i =0;i<number_of_possible_single_bit_differentials; i++)
    {
        fscanf(p, "%d %d %d %d %" SCNu64 " %" SCNu64 " %lf %d\n", diff[i].input.words, diff[i].input.bits,
            diff[i].output.words, diff[i].output.bits, &(diff[i].correlation.number_of_trials), &(diff[i].correlation.correlation_count), 
            &(diff[i].correlation), &(diff[i].correlation.is_significant));

        diff[i].input.number_of_bits = 1;
        diff[i].output.number_of_bits = 1;
    }

    fclose(p);
    return RV_SUCESS;
}


int update_single_bit_differentials_from_file(
    const char *file_name, 
    differential_t diff[MAXIMUM_NUMBER_OF_POSSIBLE_SINGLE_BIT_DIFFERENTIALS]
    )
{
    FILE *p = NULL;
    differential_t *old_diff = NULL;
    uint32_t state_size_in_bits = sizeof(uint32_t)*get_state_size(alg_type);
    uint32_t iv_size_in_bits = sizeof(uint32_t)*get_state_size(alg_type);
    uint32_t number_of_possible_single_bit_differentials = state_size_in_bits * iv_size_in_bits;

    p = fopen(file_name, "r");
    if(p == NULL) return(write_single_bit_differentials_to_file(file_name, diff));

    old_diff = (differential_t * ) malloc(sizeof(differential_t) * number_of_possible_single_bit_differentials);
    if(old_diff == NULL)
    {
        fclose(p);
        printf("Not enough memory\n");
        return RV_ERROR;
    }

    read_single_bit_differentials_from_file(file_name, old_diff, alg_type);

    for(int i=0;i<number_of_possible_single_bit_differentials;i++)
    {
        diff[i].correlation.correlation_count += old_diff[i].correlation.correlation_count;
        diff[i].correlation.number_of_trials += old_diff[i].correlation.number_of_trials;
        ct_compute_and_test_correlation(&(diff[i].correlation));
    }

    write_single_bit_differentials_to_file(file_name, diff, alg_type);
    
    fclose(p);
    free(old_diff);
    return RV_SUCESS;
}


void create_folder_if_doesnt_exist(const char *name) {
    struct stat sb;
    int e = stat(name, &sb);
    if (e != 0)
    {
        if (errno == ENOENT)
            {
            // fprintf(stderr, "The directory does not exist. Creating new directory...\n");
            // Add more flags to the mode if necessary.
            e = mkdir(name, S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH | S_IRWXU);
            if (e != 0)
                {
                fprintf(stderr, "mkdir failed; errno=%d\n",errno);
                exit(1);
                }
            }
    } else {
        if (sb.st_mode & S_IFREG)
            fprintf(stderr, "%s is a regular file.\n",name);
            exit(1);  
    }
}
