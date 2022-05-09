#include "pnb_kernels.cuh"
#include <inttypes.h>
#include "arx_cryptanalysis.cuh"

#define GENERATE_KTV()\
    K[0] = curand(&rng);\
    K[1] = curand(&rng);\
    K[2] = curand(&rng);\
    K[3] = curand(&rng);\
    K[4] = curand(&rng);\
    K[5] = curand(&rng);\
    K[6] = curand(&rng);\
    K[7] = curand(&rng);\
    nonce[0] = curand(&rng);\
    nonce[1] = curand(&rng);\
    ctr[0] = curand(&rng);\
    ctr[1] = curand(&rng);\

//---------------------------------------------------------------------------------------
//----------------------  PARTE 2
//----------------------  Kernels
//---------------------------------------------------------------------------------------

/*
Alg of sec 3.3 of New features of latin dances.
*/
__global__ void compute_neutrality_kernel(unsigned long long seed, uint32_t enc_subrounds, uint32_t dec_subrounds,
uint32_t *id, uint32_t *od, int ntest, unsigned long long *d_result, int alg_type)
{
    int word = blockIdx.y, neutral_bit = blockIdx.z, neutral_word;

    uint32_t K[8] = { 0 }, state[STATE_SIZE] = { 0 }, alt_state[STATE_SIZE] = { 0 };
    uint32_t final_state[STATE_SIZE] = { 0 }, alt_final_state[STATE_SIZE] = { 0 };
    uint32_t intermediate_state[STATE_SIZE] = { 0 }, alt_intermediate_state[STATE_SIZE] = { 0 }, aux[STATE_SIZE];
    uint32_t nonce[2] = { 0 }, ctr[2] = { 0 };
    uint32_t diff, new_diff, neutral_diff;
    int count = 0;
    curandState_t rng;
    algorithm alg;

    const unsigned long long blockId = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;
    const unsigned long long tid = blockId * blockDim.x + threadIdx.x;

    define_alg(&alg, alg_type);
    neutral_word = alg.key_positions[word];

    curand_init(seed, tid, 0, &rng);
    neutral_diff = 1 << neutral_bit;

    for (uint64_t i = 0; i < ntest; i++)
    {
        GENERATE_KTV();
        alg.init(state, K, nonce, ctr);
        xor_array(alt_state, state, id, STATE_SIZE);
        alg.encrypt_subrounds(final_state, state, enc_subrounds);
        alg.encrypt_subrounds(alt_final_state, alt_state, enc_subrounds);
        
        alg.decrypt_subrounds(final_state, state, intermediate_state, dec_subrounds, enc_subrounds);
        alg.decrypt_subrounds(alt_final_state, alt_state, alt_intermediate_state, dec_subrounds, enc_subrounds);
        xor_array(aux, intermediate_state, alt_intermediate_state, STATE_SIZE);
        diff = check_parity_of_equation(aux, od);
        
        state[neutral_word] ^= neutral_diff;
        alt_state[neutral_word] ^= neutral_diff;
        
        alg.decrypt_subrounds(final_state, state, intermediate_state, dec_subrounds, enc_subrounds);
        alg.decrypt_subrounds(alt_final_state, alt_state, alt_intermediate_state, dec_subrounds, enc_subrounds);
        xor_array(aux, intermediate_state, alt_intermediate_state, STATE_SIZE);
        new_diff = check_parity_of_equation(aux, od);
        if (diff == new_diff) 
            count++;
    }

    //printf("ntest %d count = %d\n",ntest, count);
    atomicAdd(&d_result[32 * word + neutral_bit], count);
}


/*This function computes \varepsilon_a from a PNB attack as presented in aumasson 2008*/
__global__ void compute_bias_of_g_for_random_key_kernel(
    unsigned long long seed, uint32_t enc_subrounds, uint32_t dec_subrounds,
    uint32_t *id_mask, uint32_t *od_mask,
    uint32_t *pnb, uint32_t number_of_pnb, int n_test_for_each_thread,
    unsigned long long int *d_result, int alg_type
)
{
    algorithm alg;
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    uint32_t k_with_zeros[8] = { 0 }, state[STATE_SIZE] = { 0 }, alt_state[STATE_SIZE] = { 0 };
    uint32_t final_state[STATE_SIZE] = { 0 }, alt_final_state[STATE_SIZE] = { 0 }, aux[STATE_SIZE];
    uint32_t intermediate_state[STATE_SIZE] = { 0 }, alt_intermediate_state[STATE_SIZE] = { 0 };
    uint32_t nonce[2] = { 0 }, ctr[2] = { 0 };
    curandState_t rng;
    uint32_t f_parity, g_parity;
    unsigned long long int sum_parity = 0;
    uint32_t mask;

    uint32_t k_rand[8];

    //printf("Ola sou a td %d\n", tid);
    define_alg(&alg, alg_type);
    curand_init(seed, tid, 0, &rng);

    for (int t = 0; t < n_test_for_each_thread; t++)
    {
        for (int i = 0; i < 8; i++)
        {
            k_rand[i] = curand(&rng);
            k_with_zeros[i] = k_rand[i];
        }

        for (uint32_t j = 0; j < number_of_pnb; j++)
        {
            mask = ~(1 << (pnb[j] % 32));
            k_with_zeros[pnb[j] / 32] = k_with_zeros[pnb[j] / 32] & mask;
        }

        nonce[0] = curand(&rng); nonce[1] = curand(&rng);
        ctr[0] = curand(&rng); ctr[1] = curand(&rng);

        //compute for f
        alg.init(state, k_rand, nonce, ctr);
        xor_array(alt_state, state, id_mask, STATE_SIZE);

        alg.encrypt_subrounds(final_state, state, enc_subrounds);
        alg.encrypt_subrounds(alt_final_state, alt_state, enc_subrounds);

        alg.decrypt_subrounds(final_state, state, intermediate_state, dec_subrounds, enc_subrounds);
        alg.decrypt_subrounds(alt_final_state, alt_state, alt_intermediate_state, dec_subrounds, enc_subrounds);

        xor_array(aux, intermediate_state, alt_intermediate_state, STATE_SIZE);
        f_parity = check_parity_of_equation(aux, od_mask);

        //compute for g
        alg.init(state, k_with_zeros, nonce, ctr);
        xor_array(alt_state, state, id_mask, STATE_SIZE);

        //use the same final state
        alg.decrypt_subrounds(final_state, state, intermediate_state, dec_subrounds, enc_subrounds);
        alg.decrypt_subrounds(alt_final_state, alt_state, alt_intermediate_state, dec_subrounds, enc_subrounds);

        xor_array(aux, intermediate_state, alt_intermediate_state, STATE_SIZE);
        g_parity = check_parity_of_equation(aux, od_mask);

        if (f_parity == g_parity)
            sum_parity++;
    }

    atomicAdd(d_result, sum_parity);
}



void compute_neutrality_vector(pnb_t *pnb, uint64_t number_of_trials)
{
    unsigned long long int *d_results;
    unsigned long long int h_results[KEY_SIZE_IN_BITS] = { 0 }, seed;
    uint64_t acc_results[KEY_SIZE_IN_BITS] = {0}, result[KEY_SIZE_IN_BITS] = {0};
    int lenght = KEY_SIZE_IN_BITS * sizeof(unsigned long long int);

    uint32_t *d_id, *d_od;
    srand_by_rank();
    uint64_t ntest = (1 << 10), nthreads = (1 << 8), numblocks = (1 << 1);
    uint64_t iterations = number_of_trials/ntest/nthreads/numblocks/(num_procs-1);

    if(my_rank!=0)
    {
        cudaSetDevice((my_rank-1)%NUMBER_OF_DEVICES_PER_MACHINE);
        cudaMalloc((void **)&d_results, lenght);
        cudaMalloc(&d_id, STATE_SIZE * sizeof(uint32_t));
        cudaMalloc(&d_od, STATE_SIZE * sizeof(uint32_t));
        cudaMemcpy(d_id, pnb->diff.input.mask, STATE_SIZE * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(d_od, pnb->la.output.mask, STATE_SIZE * sizeof(uint32_t), cudaMemcpyHostToDevice);

        for (int i = 0; i < iterations; i++)
        {
            memset(h_results, 0, lenght);
            cudaMemcpy(d_results, h_results, lenght, cudaMemcpyHostToDevice);
            seed = seed_by_rank();

            compute_neutrality_kernel <<< dim3(numblocks, 8, 32), nthreads >>> (seed, pnb->subrounds, 
                pnb->subrounds-pnb->la.output.subround, 
                d_id, d_od, ntest, d_results, pnb->alg_type);

            cudaMemcpy(h_results, d_results, lenght, cudaMemcpyDeviceToHost);

            for(int j=0;j<KEY_SIZE_IN_BITS; j++)
                acc_results[j]+= (uint64_t) h_results[j];
        }
        cudaFree(d_results);
        cudaFree(d_id);
        cudaFree(d_od);
    }

    MPI_Allreduce(&acc_results, &result, KEY_SIZE_IN_BITS, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);	
    
    for (int bit = 0; bit<KEY_SIZE_IN_BITS; bit++)
        pnb->neutrality_measure[bit] = 2 * ((double)result[bit]) / number_of_trials - 1;
}

void compute_correlation_of_g(pnb_t *pnb)
{
    int n_tests_for_each_thread = (1 << 10), n_threads = (1 << 8), n_blocks = (1 << 6);
    int executions_per_kernel = n_tests_for_each_thread * n_threads*n_blocks;
    uint64_t iterations;
    uint64_t result = 0, seed, local_sum=0;
    unsigned long long int *d_sum_parity;
    uint32_t *d_id, *d_od, *dPNB;
    unsigned long long int local_sum_parity = 0;

    srand_by_rank();
    iterations = pnb->correlation_of_g.number_of_trials / (executions_per_kernel) / (num_procs-1);

    if (my_rank != 0)
    {
        cudaSetDevice((my_rank-1)%NUMBER_OF_DEVICES_PER_MACHINE);

        cudaMalloc(&d_sum_parity, sizeof(unsigned long long int));
        cudaMalloc(&d_id, STATE_SIZE * sizeof(uint32_t));
        cudaMalloc(&d_od, STATE_SIZE * sizeof(uint32_t));
        cudaMalloc(&dPNB, pnb->number_of_pnb * sizeof(uint32_t));

        cudaMemcpy(d_id, pnb->diff.input.mask, STATE_SIZE * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(d_od, pnb->la.output.mask, STATE_SIZE * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(dPNB, pnb->pnb, pnb->number_of_pnb * sizeof(uint32_t), cudaMemcpyHostToDevice);

        for (int i = 0; i < iterations; i++)
        {
            seed = seed_by_rank();
            local_sum_parity = 0;
            cudaMemcpy(d_sum_parity, &local_sum_parity, sizeof(unsigned long long int), cudaMemcpyHostToDevice);

            compute_bias_of_g_for_random_key_kernel <<< n_blocks, n_threads >>> ((unsigned long long)seed,
                pnb->subrounds, pnb->subrounds - pnb->la.output.subround, d_id, d_od, dPNB, 
                pnb->number_of_pnb, n_tests_for_each_thread, d_sum_parity, pnb->alg_type);

            cudaMemcpy(&local_sum_parity, d_sum_parity, sizeof(unsigned long long int), cudaMemcpyDeviceToHost);
            local_sum += (uint64_t) local_sum_parity;
        }
        cudaFree(d_sum_parity);
        cudaFree(d_id);
        cudaFree(d_od); 
    }

    MPI_Allreduce(&local_sum, &result, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD); 

    pnb->correlation_of_g.correlation_count = result;	
    ct_compute_and_test_correlation(&(pnb->correlation_of_g));
}


void get_pnb_list(pnb_t *pnb)
{
    pnb->number_of_pnb = 0;
    for (int i = 0; i < KEY_SIZE_IN_BITS; i++)
    {
        if (pnb->neutrality_measure[i] >= pnb->threshold)
        {
            pnb->pnb[pnb->number_of_pnb] = i;
            pnb->number_of_pnb++;
        }
    }
}

double get_max(double x, double y)
{
    if (x>y)
        return x;
    else
        return y;
}


//compute complexity of pnb attack
void compute_complexity_of_the_attack(pnb_t *pnb)
{
    int alpha;
    int m = KEY_SIZE_IN_BITS - pnb->number_of_pnb;
    double N, tc, minN = KEY_SIZE_IN_BITS, minTC = KEY_SIZE_IN_BITS;
    double bias_of_g = pnb->diff.correlation.observed * pow(pnb->la.correlation.observed,2) * pnb->correlation_of_g.observed;

    for (alpha = 0; alpha < KEY_SIZE_IN_BITS; alpha++)
    {
        N = (sqrt(alpha*log(4)) + 3 * sqrt(1 - bias_of_g * bias_of_g)) / bias_of_g;
        N = N * N;
        tc = get_max(KEY_SIZE_IN_BITS - alpha, m + log2(N));

        if (tc < minTC)
        {
            minTC = tc;
            minN = N;
        }
    }

    pnb->data_complexity = log2(minN);
    pnb->time_complexity = minTC;
}



void pnb_attack_for_single_bit_differential(
    int idw, 
    int idb, 
    int odw, 
    int odb, 
    int subrounds,
    int differential_part_subrounds,
    int linear_part_subrounds, 
    double threshold,
    int alg_type,
    FILE *output_file
    )
{
    pnb_t pnb = {0};
    
    pnb.subrounds = subrounds;
    pnb.threshold = threshold;
    pnb.alg_type = alg_type;
    differential_compute_from_single_bit(&pnb.diff, idw, idb, odw, odb, differential_part_subrounds, alg_type);
    la_compute_from_differential(&pnb.la, pnb.diff, linear_part_subrounds);

    compute_neutrality_vector(&pnb, (1<<26));
    get_pnb_list(&pnb);

    pnb.correlation_of_g.number_of_trials = 1;
    pnb.correlation_of_g.number_of_trials <<= 34;
    compute_correlation_of_g(&pnb);
        
    compute_complexity_of_the_attack(&pnb);

    if(my_rank == 0)
        pnb_print(output_file, pnb);
}