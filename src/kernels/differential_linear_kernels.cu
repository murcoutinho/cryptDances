#include "differential_linear_kernels.cuh"
#include <inttypes.h>
#include "types.cuh"
#include "arx_cryptanalysis.cuh"

void get_differential_from_position(int position, differential_t *diff)
{
    uint32_t iv_positions[4];
    uint32_t state_size_in_bits = sizeof(uint32_t)*get_state_size(diff->alg_type);

    diff->input.words[0] = position/(32*state_size_in_bits);
    position -= diff->input.words[0] * 32 * state_size_in_bits;

    diff->input.bits[0] = position/state_size_in_bits;
    position -= diff->input.bits[0] * state_size_in_bits;

    diff->output.words[0] = position/32;
    position -= diff->output.words[0] * 32;

    diff->output.bits[0] = position;

    get_alg_iv_positions(iv_positions, diff->alg_type);
    diff->input.words[0] = iv_positions[diff->input.words[0]];

    diff->input.number_of_bits = 1;
    diff->output.number_of_bits = 1;

    lob_compute_mask_from_list_of_bits(&(diff->input));
    lob_compute_mask_from_list_of_bits(&(diff->output));
}

/*
    Computes the differential correlation for every output bit for every 128 possible single input differentials.
    The test executes I iterations for n keys.

    The organization is the following:
        - blockIdx.y - indicate which word for the input differential, from 0 to 3.
        - blockIdx.z - indicate which bit for the input differential, from 0 to 31.
        - blockIdx.x * thread.x * ntest_for_each_key - number of tests

 * . For increased performance we tried to minimize all unnecessary operations, 
 * including some that were not very straightforward. 
 * . For instance, the initialization of the state is bypassed and a completely 
 * random state is considered. This is useful to avoid copying memory from one place to another.
 * . Additionally, the CUDA PRNG is called only once to initialize the first state. 
 * After that, the previous resulting state is considered as a starting point since it is itself random.
 * . These choices do not seem to affect the results, as the replication of previous papers shows.
 * . Therefore, this kernel is not useful for researchers that are trying to analyze some particular 
 * behaviours of the constants of the algorithms.
*/
__global__ void differential_correlation_exhaustion_kernel(unsigned long long seed, 
int subrounds, int last_subround, int n_test_for_each_thread, unsigned long long *d_result, int alg_type)
{
    uint32_t id[MAXIMUM_STATE_SIZE] = { 0 };
    algorithm alg;
    uint32_t observed_od[MAXIMUM_STATE_SIZE];
    uint8_t observed_od_bits[MAXIMUM_STATE_SIZE_IN_BITS];
    curandState_t rng;
    int sum[MAXIMUM_STATE_SIZE_IN_BITS] = { 0 };
    uint32_t state[MAXIMUM_STATE_SIZE] = { 0 }, alt_state[MAXIMUM_STATE_SIZE] = { 0 };

    int word = blockIdx.y, bit = blockIdx.z;
    const unsigned long long blockId = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;
    const unsigned long long tid = blockId * blockDim.x + threadIdx.x;

    uint32_t state_size_in_bits = sizeof(uint32_t)*get_state_size(alg_type);
    define_alg(&alg, alg_type);
    curand_init(seed, tid, 0, &rng);

    //comput id - each block may test a different id
    id[alg.iv_positions[word]] = 1 << bit;

    for (int t = 0; t < n_test_for_each_thread; t++)
    {
        GENERATE_RANDOM_STATE(state,alg.state_size);
        xor_array(alt_state, state, id, MAXIMUM_STATE_SIZE);
        alg.subrounds(state, subrounds,last_subround);
        alg.subrounds(alt_state, subrounds,last_subround);
        xor_array(observed_od, state, alt_state, MAXIMUM_STATE_SIZE);
        transform_state_to_bits(observed_od, observed_od_bits);
        update_result(sum, observed_od_bits);
    }

    for (int i = 0; i < state_size_in_bits; i++)
        atomicAdd(&d_result[32 * state_size_in_bits * word + state_size_in_bits * bit + i], sum[i]);
}


void compute_all_single_bit_differential_correlation(int alg_type, int subrounds, int last_subround, 
uint64_t number_of_trials, const char *out_file_name)
{
    unsigned long long int *d_results;
    uint64_t numblocks = NUMBER_OF_CUDA_BLOCKS/(4*32);
    uint64_t iterations = number_of_trials / NUMBER_OF_CUDA_THREADS / numblocks / NUMBER_OF_TESTS_PER_THREAD / (num_procs-1);
    unsigned long long int h_results[MAXIMUM_NUMBER_OF_POSSIBLE_SINGLE_BIT_DIFFERENTIALS] = { 0 }, seed;
    uint64_t acc_results[MAXIMUM_NUMBER_OF_POSSIBLE_SINGLE_BIT_DIFFERENTIALS] = {0}, result[MAXIMUM_NUMBER_OF_POSSIBLE_SINGLE_BIT_DIFFERENTIALS] = {0};
    uint32_t state_size_in_bits = sizeof(uint32_t)*get_state_size(alg_type);
    uint32_t iv_size_in_bits = sizeof(uint32_t)*get_state_size(alg_type);
    uint32_t number_of_possible_single_bit_differentials = state_size_in_bits * iv_size_in_bits;
    int L = number_of_possible_single_bit_differentials * sizeof(unsigned long long int);

    srand_by_rank(); //initialize prng with different internal state for each MPI process

    memset(result,0x00,MAXIMUM_NUMBER_OF_POSSIBLE_SINGLE_BIT_DIFFERENTIALS*sizeof(uint64_t));
    memset(acc_results,0x00,MAXIMUM_NUMBER_OF_POSSIBLE_SINGLE_BIT_DIFFERENTIALS*sizeof(uint64_t));
    
    if (my_rank != 0)
    {
        cudaSetDevice((my_rank-1)%NUMBER_OF_DEVICES_PER_MACHINE);

        cudaMalloc((void **)&d_results, L);
        cudaMemcpy(d_results, h_results, L, cudaMemcpyHostToDevice);

        for (int i = 0; i < iterations; i++)
        {
            memset(h_results, 0, L);
            cudaMemcpy(d_results, h_results, L, cudaMemcpyHostToDevice);
            seed = seed_by_rank();

            differential_correlation_exhaustion_kernel <<<dim3(numblocks, 4, 32), NUMBER_OF_CUDA_THREADS>>> (seed, 
                subrounds, last_subround, NUMBER_OF_TESTS_PER_THREAD, d_results, alg_type);

            cudaMemcpy(h_results, d_results, L, cudaMemcpyDeviceToHost);

            for(int j=0;j<number_of_possible_single_bit_differentials; j++)
                acc_results[j]+= (uint64_t) h_results[j];
        }

        cudaFree(d_results);
    }

    MPI_Allreduce(&acc_results, &result, number_of_possible_single_bit_differentials, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);	

    if(my_rank == 0)
    {
        differential_t *diff = NULL;
        diff = (differential_t * ) malloc(sizeof(differential_t) * number_of_possible_single_bit_differentials);
        if(diff == NULL)
        {
            printf("Not enough memory\n");
            return;
        }
        for(int position = 0; position < number_of_possible_single_bit_differentials; position++)
        {
            memset(&diff[position],0x00, sizeof(differential_t));
            diff[position].alg_type = alg_type;
            get_differential_from_position(position, &diff[position]);
            diff[position].correlation.number_of_trials = number_of_trials;
            diff[position].correlation.correlation_count = number_of_trials - result[position];
            ct_compute_and_test_correlation(&(diff[position].correlation));
        }

        update_single_bit_differentials_from_file(out_file_name, diff, alg_type);

        free(diff);
    }
}


/** 
 * Computes the differential bias given the number of subrounds, ID and OD mask. Some important observations:
 * . For increased performance we tried to minimize all unnecessary operations, 
 * including some that were not very straightforward. 
 * . For instance, the initialization of the state is bypassed and a completely 
 * random state is considered. This is useful to avoid copying memory from one place to another.
 * . Additionally, the CUDA PRNG is called only once to initialize the first state. 
 * After that, the previous resulting state is considered as a starting point since it is itself random.
 * . These choices do not seem to affect the results, as the replication of previous papers shows.
 * . Therefore, this kernel is not useful for researchers that are trying to analyze some particular 
 * behaviours of the constants of the algorithms.
 * */
__global__ void differential_correlation_kernel(unsigned long long seed, int subrounds, int last_subround, uint32_t *id,
    uint32_t *od, int n_test_for_each_thread, unsigned long long int *d_result, int alg_type)
{
    algorithm alg;
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    uint32_t state[MAXIMUM_STATE_SIZE] = { 0 }, alt_state[MAXIMUM_STATE_SIZE] = { 0 }, observed_od[MAXIMUM_STATE_SIZE] = {0};
    curandState_t rng;
    unsigned long long int sum_parity = 0;

    define_alg(&alg, alg_type);
    curand_init(seed, tid, 0, &rng);

    for (int t = 0; t < n_test_for_each_thread; t++)
    {	
        GENERATE_RANDOM_STATE(state,alg.state_size);
        xor_array(alt_state, state, id, alg.state_size);
        alg.subrounds(state, subrounds,last_subround);
        alg.subrounds(alt_state, subrounds,last_subround);
        xor_array(observed_od, state, alt_state, alg.state_size);
        sum_parity += check_parity_of_equation(observed_od, od, alg.state_size);
    }

    // if(tid == 0)
    // {
    //     for(int i=0;i<alg.state_size;i++){
    //         printf("state[i]=%08X\n", state[i]);
    //         printf("seed = %016lx\n", seed);
    //         printf("alt_state[i]=%08X\n", alt_state[i]);
    //         printf("observed_od[i]=%08X\n", observed_od[i]);
    //         printf("id[i]=%08X\n", id[i]);
    //         printf("od[i]=%08X\n", od[i]);
    //     }

    // }

    atomicAdd(d_result, sum_parity);
}


/*
Computes the linear bias given the number of rounds, ID and OD mask.
*/
__global__ void linear_correlation_kernel(unsigned long long seed, int subrounds, 
    int last_subround, uint32_t *id_mask, uint32_t *od_mask, int n_test_for_each_thread, 
    unsigned long long int *d_result, int alg_type)
{
    algorithm alg;
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    uint32_t input_state[MAXIMUM_STATE_SIZE] = { 0 }, output_state[MAXIMUM_STATE_SIZE] = { 0 };
    curandState_t rng;
    unsigned long long int sum_parity = 0;

    define_alg(&alg, alg_type);
    curand_init(seed, tid, 0, &rng);

    for (int t = 0; t < n_test_for_each_thread; t++)
    {
        GENERATE_RANDOM_STATE(output_state,alg.state_size);
        for(int i=0;i<alg.state_size;i++)
            input_state[i] = output_state[i];
        alg.subrounds(output_state, subrounds, last_subround);
        sum_parity += check_parity_of_linear_relation(id_mask, input_state, od_mask, output_state, alg.state_size);
    }

    atomicAdd(d_result, sum_parity);
}



void compute_differential_or_linear_correlation(diff_lin_t *diff_lin, int type)
{
    uint64_t result = 0, seed, local_sum=0;
    uint64_t iterations;
    unsigned long long int *d_sum_parity;
    uint32_t *d_id, *d_od;
    unsigned long long int local_sum_parity = 0;
    uint32_t state_size = get_state_size(diff_lin->alg_type);

    int subrounds = diff_lin->output.subround - diff_lin->input.subround;
    int last_subround = diff_lin->input.subround;
    
    srand_by_rank(); //initialize prng with different internal state for each MPI process
    iterations = diff_lin->correlation.number_of_trials / TOTAL_EXECUTIONS_PER_KERNEL / (num_procs-1);
    
    if (my_rank != 0)
    {
        cudaSetDevice((my_rank-1)%NUMBER_OF_DEVICES_PER_MACHINE);
 
        cudaMalloc(&d_sum_parity, sizeof(unsigned long long int));
        cudaMalloc(&d_id, state_size * sizeof(uint32_t));
        cudaMalloc(&d_od, state_size * sizeof(uint32_t));

        cudaMemcpy(d_id, diff_lin->input.mask, state_size * sizeof(uint32_t), cudaMemcpyHostToDevice);
        cudaMemcpy(d_od, diff_lin->output.mask, state_size * sizeof(uint32_t), cudaMemcpyHostToDevice); 

        for (int i = 0; i < iterations; i++)
        {
            seed = seed_by_rank();
            local_sum_parity = 0;
            cudaMemcpy(d_sum_parity, &local_sum_parity, 
                sizeof(unsigned long long int), cudaMemcpyHostToDevice);

            if(type == TYPE_DIFFERENTIAL)
                differential_correlation_kernel <<< NUMBER_OF_CUDA_BLOCKS, NUMBER_OF_CUDA_THREADS >>> ((unsigned long long)seed, 
                    subrounds, last_subround, d_id, d_od, NUMBER_OF_TESTS_PER_THREAD, d_sum_parity, diff_lin->alg_type);
            else
                linear_correlation_kernel <<< NUMBER_OF_CUDA_BLOCKS, NUMBER_OF_CUDA_THREADS >>> ((unsigned long long)seed, 
                    subrounds, last_subround, d_id, d_od, NUMBER_OF_TESTS_PER_THREAD, d_sum_parity, diff_lin->alg_type);

            cudaMemcpy(&local_sum_parity, d_sum_parity, 
                sizeof(unsigned long long int), cudaMemcpyDeviceToHost);
            local_sum += (uint64_t) local_sum_parity;
        }
    
        cudaFree(d_sum_parity);
        cudaFree(d_id);
        cudaFree(d_od);
    }

    printf("local sum = %ld\n", local_sum);
    MPI_Allreduce(&local_sum,&result,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);   

    diff_lin->correlation.correlation_count = diff_lin->correlation.number_of_trials-result;	
    ct_compute_and_test_correlation(&(diff_lin->correlation));
}


void search_until_find_correlation(diff_lin_t *diff_lin, int type)
{
    uint64_t count = 0, total = 0;
    int level = 34;
    double correlation = 0;

    diff_lin->correlation.number_of_trials = 1;
    diff_lin->correlation.number_of_trials <<= level;

    while (1)
    {
        compute_differential_or_linear_correlation(diff_lin, type);
        count += diff_lin->correlation.correlation_count;
        total += diff_lin->correlation.number_of_trials;

        correlation = ((double)count) / total;  
        correlation = 2 * correlation - 1;

        if(test_significance_of_correlation(correlation, total))
            break;
        if(level>MAX_LEVEL)
        {
            correlation = 0;
            break;
        }
        
        diff_lin->correlation.number_of_trials = 1;
        diff_lin->correlation.number_of_trials <<= level; //first repeat previous level since 2^level + 2^level = 2^{level+1}
        level++;
    }

    diff_lin->correlation.observed = correlation;
    diff_lin->correlation.number_of_trials = total;
    diff_lin->correlation.correlation_count = count;
    diff_lin->correlation.is_significant = TRUE;
}