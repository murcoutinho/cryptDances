#include "pnb_kernels.cuh"
#include <inttypes.h>
#include <stdint.h>
#include "arx_cryptanalysis.cuh"
#include <iostream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <cassert>
#include <iomanip>
#include <chrono>

using namespace std::chrono;
#include "median.cuh"

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

__global__ void compute_neutrality_kernel(
        unsigned long long seed, uint32_t enc_subrounds, uint32_t dec_subrounds,
        uint32_t *id, uint32_t *od, int ntest, unsigned long long *d_result, int alg_type, int my_rank, int iteration,
        int nthreads, int numblocks, unsigned long long int n_trials_per_proc
)
{
    int word = blockIdx.y, neutral_bit = blockIdx.z, neutral_word;

    uint32_t K[8] = { 0 }, state[MAXIMUM_STATE_SIZE] = { 0 }, alt_state[MAXIMUM_STATE_SIZE] = { 0 };
    uint32_t final_state[MAXIMUM_STATE_SIZE] = { 0 }, alt_final_state[MAXIMUM_STATE_SIZE] = { 0 };
    uint32_t intermediate_state[MAXIMUM_STATE_SIZE] = { 0 }, alt_intermediate_state[MAXIMUM_STATE_SIZE] = { 0 }, aux[MAXIMUM_STATE_SIZE];
    uint32_t nonce[2] = { 0 }, ctr[2] = { 0 };
    uint32_t diff, new_diff, neutral_diff;
    int count = 0;
    curandState_t rng;
    algorithm alg;

    const unsigned long long blockId = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;
    const unsigned long long tid = blockId * blockDim.x + threadIdx.x;

    define_alg(&alg, alg_type);
    neutral_word = alg.key_positions[word];

    curand_init(seed, n_trials_per_proc*my_rank + iteration*nthreads*numblocks*8*32+tid, 0, &rng);
    neutral_diff = 1;
    neutral_diff <<= neutral_bit;

    for (uint64_t i = 0; i < ntest; i++)
    {
        GENERATE_KTV();
        alg.init(state, K, nonce, ctr);
        xor_array(alt_state, state, id, MAXIMUM_STATE_SIZE);
        alg.encrypt_subrounds(final_state, state, enc_subrounds);
        alg.encrypt_subrounds(alt_final_state, alt_state, enc_subrounds);

        alg.decrypt_subrounds(final_state, state, intermediate_state, dec_subrounds, enc_subrounds);
        alg.decrypt_subrounds(alt_final_state, alt_state, alt_intermediate_state, dec_subrounds, enc_subrounds);
        xor_array(aux, intermediate_state, alt_intermediate_state, MAXIMUM_STATE_SIZE);
        diff = check_parity_of_equation(aux, od, alg.state_size);

        state[neutral_word] ^= neutral_diff;
        alt_state[neutral_word] ^= neutral_diff;

        alg.decrypt_subrounds(final_state, state, intermediate_state, dec_subrounds, enc_subrounds);
        alg.decrypt_subrounds(alt_final_state, alt_state, alt_intermediate_state, dec_subrounds, enc_subrounds);
        xor_array(aux, intermediate_state, alt_intermediate_state, MAXIMUM_STATE_SIZE);
        new_diff = check_parity_of_equation(aux, od, alg.state_size);
        if (diff == new_diff)
            count++;
    }
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
    uint32_t k_with_zeros[8] = { 0 }, state[MAXIMUM_STATE_SIZE] = { 0 }, alt_state[MAXIMUM_STATE_SIZE] = { 0 };
    uint32_t final_state[MAXIMUM_STATE_SIZE] = { 0 }, alt_final_state[MAXIMUM_STATE_SIZE] = { 0 }, aux[MAXIMUM_STATE_SIZE];
    uint32_t intermediate_state[MAXIMUM_STATE_SIZE] = { 0 }, alt_intermediate_state[MAXIMUM_STATE_SIZE] = { 0 };
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
        xor_array(alt_state, state, id_mask, MAXIMUM_STATE_SIZE);

        alg.encrypt_subrounds(final_state, state, enc_subrounds);
        alg.encrypt_subrounds(alt_final_state, alt_state, enc_subrounds);

        alg.decrypt_subrounds(final_state, state, intermediate_state, dec_subrounds, enc_subrounds);
        alg.decrypt_subrounds(alt_final_state, alt_state, alt_intermediate_state, dec_subrounds, enc_subrounds);

        xor_array(aux, intermediate_state, alt_intermediate_state, MAXIMUM_STATE_SIZE);
        f_parity = check_parity_of_equation(aux, od_mask, alg.state_size);

        //compute for g
        alg.init(state, k_with_zeros, nonce, ctr);
        xor_array(alt_state, state, id_mask, MAXIMUM_STATE_SIZE);

        //use the same final state
        alg.decrypt_subrounds(final_state, state, intermediate_state, dec_subrounds, enc_subrounds);
        alg.decrypt_subrounds(alt_final_state, alt_state, alt_intermediate_state, dec_subrounds, enc_subrounds);

        xor_array(aux, intermediate_state, alt_intermediate_state, MAXIMUM_STATE_SIZE);
        g_parity = check_parity_of_equation(aux, od_mask, alg.state_size);

        if (f_parity == g_parity)
            sum_parity++;
    }

    atomicAdd(d_result, sum_parity);
}

__global__ void compute_bias_of_g_for_random_key_kernel_using_median(
        unsigned long long seed, uint32_t enc_subrounds, uint32_t dec_subrounds,
        uint32_t *id_mask, uint32_t *od_mask,
        uint32_t *pnb, uint32_t number_of_pnb, int n_test_for_each_thread,
        int alg_type, unsigned long long int * medians, int iteration,
        unsigned long long int blocks, unsigned long long int threads
)
{
    algorithm alg;
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    uint32_t k_with_zeros[8] = { 0 }, state[MAXIMUM_STATE_SIZE] = { 0 }, alt_state[MAXIMUM_STATE_SIZE] = { 0 };
    uint32_t final_state[MAXIMUM_STATE_SIZE] = { 0 }, alt_final_state[MAXIMUM_STATE_SIZE] = { 0 }, aux[MAXIMUM_STATE_SIZE];
    uint32_t intermediate_state[MAXIMUM_STATE_SIZE] = { 0 }, alt_intermediate_state[MAXIMUM_STATE_SIZE] = { 0 };
    uint32_t nonce[2] = { 0 }, ctr[2] = { 0 };
    curandState_t rng;
    uint32_t f_parity, g_parity;
    unsigned long long int sum_parity = 0;
    uint32_t mask;

    uint32_t k_rand[8];

    //printf("Ola sou a td %d\n", tid);
    define_alg(&alg, alg_type);
    curand_init(seed, tid, 0, &rng);

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

    for (int t = 0; t < n_test_for_each_thread; t++)
    {
        nonce[0] = curand(&rng); nonce[1] = curand(&rng);
        ctr[0] = curand(&rng); ctr[1] = curand(&rng);

        //compute for f
        alg.init(state, k_rand, nonce, ctr);
        xor_array(alt_state, state, id_mask, MAXIMUM_STATE_SIZE);

        alg.encrypt_subrounds(final_state, state, enc_subrounds);
        alg.encrypt_subrounds(alt_final_state, alt_state, enc_subrounds);

        alg.decrypt_subrounds(final_state, state, intermediate_state, dec_subrounds, enc_subrounds);
        alg.decrypt_subrounds(alt_final_state, alt_state, alt_intermediate_state, dec_subrounds, enc_subrounds);

        xor_array(aux, intermediate_state, alt_intermediate_state, MAXIMUM_STATE_SIZE);
        f_parity = check_parity_of_equation(aux, od_mask, alg.state_size);

        //compute for g
        alg.init(state, k_with_zeros, nonce, ctr);
        xor_array(alt_state, state, id_mask, MAXIMUM_STATE_SIZE);

        //use the same final state
        alg.decrypt_subrounds(final_state, state, intermediate_state, dec_subrounds, enc_subrounds);
        alg.decrypt_subrounds(alt_final_state, alt_state, alt_intermediate_state, dec_subrounds, enc_subrounds);

        xor_array(aux, intermediate_state, alt_intermediate_state, MAXIMUM_STATE_SIZE);
        g_parity = check_parity_of_equation(aux, od_mask,alg.state_size);

        if (f_parity == g_parity)
            sum_parity++;
    }
    medians[iteration*(blocks*threads)+tid] = sum_parity;//fabs(2.0*sum_parity/n_test_for_each_thread-1);

}


void compute_neutrality_vector(pnb_t *pnb, uint64_t number_of_trials)
{

    unsigned long long int *d_results;
    unsigned long long int h_results[KEY_SIZE_IN_BITS] = { 0 };
    uint64_t acc_results[KEY_SIZE_IN_BITS] = {0}, result[KEY_SIZE_IN_BITS] = {0}, seed;
    int lenght = KEY_SIZE_IN_BITS * sizeof(unsigned long long int);

    uint32_t *d_id, *d_od;
    srand_by_rank();

    uint64_t ntest = (1 << 10), nthreads = (1 << 8), numblocks = (1 << 1);
    uint64_t iterations = number_of_trials/ntest/nthreads/numblocks/(num_procs-1);
    if (iterations % 2 != 0) {
        printf("For now (iterations % 2) should be 0\n");
        MPI_Finalize();
        exit(0);
    }
    if(my_rank!=0)
    {
        cudaSetDevice((my_rank-1)%NUMBER_OF_DEVICES_PER_MACHINE);
        cudaMalloc((void **)&d_results, lenght);
        cudaMalloc(&d_id, MAXIMUM_STATE_SIZE * sizeof(uint32_t));
        cudaMalloc(&d_od, MAXIMUM_STATE_SIZE * sizeof(uint32_t));
        cudaMemcpy(d_id, pnb->diff.input.mask, MAXIMUM_STATE_SIZE * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(d_od, pnb->la.output.mask, MAXIMUM_STATE_SIZE * sizeof(uint32_t), cudaMemcpyHostToDevice);

        for (int i = 0; i < iterations; i++)
        {
            memset(h_results, 0, lenght);
            cudaMemcpy(d_results, h_results, lenght, cudaMemcpyHostToDevice);
            seed = seed_by_rank(i);
            compute_neutrality_kernel <<< dim3(numblocks, 8, 32), nthreads >>> (seed, pnb->subrounds,
                                                                                pnb->subrounds-pnb->la.output.subround,
                                                                                d_id, d_od, ntest, d_results, pnb->alg_type, my_rank, i, nthreads, numblocks, number_of_trials/(num_procs-1));

            cudaMemcpy(h_results, d_results, lenght, cudaMemcpyDeviceToHost);
            for(int j=0; j<KEY_SIZE_IN_BITS; j++)
                acc_results[j]+= (uint64_t) h_results[j];
        }
        cudaFree(d_results);
        cudaFree(d_id);
        cudaFree(d_od);
    }
    MPI_Allreduce(&acc_results, &result, KEY_SIZE_IN_BITS, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

    for (int bit = 0; bit<KEY_SIZE_IN_BITS; bit++) {
        pnb->neutrality_measure[bit] = 2 * ((double) result[bit]) / number_of_trials - 1;
    }
}


void compute_correlation_of_g_using_mean(pnb_t *pnb)
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
        cudaMalloc(&d_id, MAXIMUM_STATE_SIZE * sizeof(uint32_t));
        cudaMalloc(&d_od, MAXIMUM_STATE_SIZE * sizeof(uint32_t));
        cudaMalloc(&dPNB, pnb->number_of_pnb * sizeof(uint32_t));

        cudaMemcpy(d_id, pnb->diff.input.mask, MAXIMUM_STATE_SIZE * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(d_od, pnb->la.output.mask, MAXIMUM_STATE_SIZE * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(dPNB, pnb->pnb, pnb->number_of_pnb * sizeof(uint32_t), cudaMemcpyHostToDevice);

        for (int i = 0; i < iterations; i++)
        {
            seed = seed_by_rank(i);
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

void compute_correlation_of_g_using_median(pnb_t *pnb)
{
    MPI_Group world_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    int ranks[num_procs];
    for (int i = 1; i < num_procs; i++)
        ranks[i-1] = i;

    // Construct a group containing all of the prime ranks in world_group
    MPI_Group new_group;
    MPI_Group_incl(world_group, num_procs - 1, ranks, &new_group);

    // Create a new communicator based on the group
    //MPI_Comm new_comm;
    MPI_Comm_create_group(MPI_COMM_WORLD, new_group, 0, &new_comm);

    int n_tests_for_each_thread = (1 << 15), n_threads = (1 << 7), n_blocks = (1 << 7);
    int executions_per_kernel = n_tests_for_each_thread * n_threads*n_blocks;
    uint64_t iterations, seed;
    uint32_t *d_id, *d_od, *dPNB;
    unsigned long long int * d_medians;

    srand_by_rank();
    // TODO:: Consider even num_procs
    iterations = pnb->correlation_of_g.number_of_trials / (executions_per_kernel) / (num_procs - 1);
    if (iterations % 2 != 0) {
        printf("For now (iterations % 2) should be 0.\n");
        exit(0);
    }
    unsigned long long parmed;
    if (my_rank != 0) {
        unsigned long long int *local_medians = (unsigned long long int *) malloc(
                iterations * n_threads * n_blocks * sizeof(unsigned long long int));
        memset(local_medians, 0, iterations * n_threads * n_blocks * sizeof(unsigned long long int));

        cudaSetDevice((my_rank-1) % NUMBER_OF_DEVICES_PER_MACHINE);
        cudaStream_t stream;
        cudaStreamCreate(&stream);

        cudaMalloc(&d_medians, sizeof(unsigned long long int) * iterations * n_threads * n_blocks);
        cudaMalloc(&d_id, MAXIMUM_STATE_SIZE * sizeof(uint32_t));
        cudaMalloc(&d_od, MAXIMUM_STATE_SIZE * sizeof(uint32_t));
        cudaMalloc(&dPNB, pnb->number_of_pnb * sizeof(uint32_t));

        cudaMemcpyAsync(d_id, pnb->diff.input.mask, MAXIMUM_STATE_SIZE * sizeof(uint32_t), cudaMemcpyHostToDevice, stream);
        cudaMemcpyAsync(d_od, pnb->la.output.mask, MAXIMUM_STATE_SIZE * sizeof(uint32_t), cudaMemcpyHostToDevice, stream);
        cudaMemcpyAsync(dPNB, pnb->pnb, pnb->number_of_pnb * sizeof(uint32_t), cudaMemcpyHostToDevice, stream);

        for (int i = 0; i < iterations; i++) {
            seed = seed_by_rank(i);

            compute_bias_of_g_for_random_key_kernel_using_median <<< n_blocks, n_threads, 0, stream >>>(
                    (unsigned long long) seed,
                    pnb->subrounds, pnb->subrounds - pnb->la.output.subround, d_id, d_od, dPNB,
                    pnb->number_of_pnb, n_tests_for_each_thread, pnb->alg_type,
                    d_medians, i, n_blocks, n_threads
            );
        }
        cudaFree(d_id);
        cudaFree(d_od);
        cudaMemcpyAsync(local_medians, d_medians, sizeof(unsigned long long int) * iterations * n_threads * n_blocks,
                        cudaMemcpyDeviceToHost, stream);
        cudaFree(d_medians);
        cudaStreamSynchronize(stream);
        cudaStreamDestroy(stream);
        vector<unsigned long long int> v(local_medians, local_medians + iterations * n_threads * n_blocks);
        v.pop_back();
        parmed = par::median(v);
        MPI_Send(&parmed, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        free(local_medians);
    } else {
        MPI_Recv(&parmed, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    pnb->correlation_of_g.correlation_count = parmed;
    ct_compute_and_test_correlation_using_median(&(pnb->correlation_of_g), n_tests_for_each_thread);
    MPI_Group_free(&world_group);
    MPI_Group_free(&new_group);

    if (MPI_COMM_NULL != new_comm) {
        MPI_Comm_free(&new_comm);
    }
}

void compute_correlation_of_g(pnb_t *pnb)
{
    if(pnb->statistic_type == STAT_MEAN)
        compute_correlation_of_g_using_mean(pnb);
    else
        compute_correlation_of_g_using_median(pnb);
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

        //Note: uncomment to get new formula from euro2020
        //tc = get_max(tc, KEY_SIZE_IN_BITS - m);

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
    uint64_t number_of_trials_for_neutrality,
    uint64_t number_of_trials_for_bias_of_g,
    int alg_type,
    int statistic_type,
    FILE *output_file
    )
{
    pnb_t pnb = {0};
    pnb.subrounds = subrounds;
    pnb.threshold = threshold;
    pnb.alg_type = alg_type;
    pnb.statistic_type = statistic_type;
    differential_compute_from_single_bit(&pnb.diff, idw, idb, odw, odb, differential_part_subrounds, alg_type);
    la_compute_from_differential(&pnb.la, pnb.diff, linear_part_subrounds);
    compute_neutrality_vector(&pnb, number_of_trials_for_neutrality);
    get_pnb_list(&pnb);
    pnb.correlation_of_g.number_of_trials = number_of_trials_for_bias_of_g;
    compute_correlation_of_g(&pnb);
    compute_complexity_of_the_attack(&pnb);

    if(my_rank == 0)
        pnb_print(output_file, pnb);
}


void pnb_remove(pnb_t *pnb, int position)
{
    int count = 0, current;
    for(int i=0;i<pnb->number_of_pnb;i++)
    {
        current = pnb->pnb[i];
        if(current != position)
            pnb->pnb[count++] = current;
    }
    pnb->number_of_pnb = count;
}

void pnb_add(pnb_t *pnb, int position)
{
    pnb->pnb[pnb->number_of_pnb++] = position;
}

void pnb_divide_groups(pnb_t *current, pnb_t *reserve, double threshold)
{
    int aux[KEY_SIZE_IN_BITS] = {0}, n = current->number_of_pnb;
    memcpy(aux, current->pnb, current->number_of_pnb*sizeof(int));

    current->number_of_pnb = 0;
    reserve->number_of_pnb = 0;

    for(int i=0; i< n; i++)
    {
        if(current->neutrality_measure[aux[i]] > threshold)
            current->pnb[current->number_of_pnb++] = aux[i];
        else
            reserve->pnb[reserve->number_of_pnb++] = aux[i];
    }
}

//This method is the step 2 of the technique proposed in Eurocrypt 2022. 
//In cryptdances we assume good computational resources. 
//Therefore, step 3 is ignored, i.e., all bits are selected with step 2.
void pnb_iteractive_selection(
    pnb_t *pnb, 
    double threshold_direct, 
    uint64_t number_of_trials_for_bias_of_g, 
    int number_of_bits_selected
    )
{
    pnb_t pnb_reserve = {0};
    int maxPos = 0;
    double maxCorrelation = 0;

    pnb->correlation_of_g.number_of_trials = number_of_trials_for_bias_of_g;
    memcpy(&pnb_reserve, pnb, sizeof(pnb_t));
    pnb_divide_groups(pnb, &pnb_reserve, threshold_direct);

    for(int n2 = 0; n2<number_of_bits_selected; n2++)
    {
        maxPos = 0; maxCorrelation = 0;
        for(int i=0; i<pnb_reserve.number_of_pnb; i++)
        {
            pnb_add(pnb, pnb_reserve.pnb[i]);
            compute_correlation_of_g(pnb);
            pnb_remove(pnb, pnb_reserve.pnb[i]);
            if(fabs(pnb->correlation_of_g.observed) > maxCorrelation)
            {
                maxPos = i;
                maxCorrelation = fabs(pnb->correlation_of_g.observed);
            }
        }
        pnb_add(pnb, pnb_reserve.pnb[maxPos]);
        pnb_remove(&pnb_reserve, pnb_reserve.pnb[maxPos]);
        //if(my_rank == 0) {printf("Max is %d\n", maxPos); pnb_print(NULL, *pnb); pnb_print(NULL, pnb_reserve);}
    }
}

