#include "tests.cuh"

int my_rank, num_procs;

__global__ void setup_kernel(curandState *state, int seed)
{
    int id = threadIdx.x + blockIdx.x * blockDim.x;
    /* Each thread gets same seed, a different sequence
       number, no offset */
    curand_init(seed, id, 0, &state[id]);
}

int main()
{
    int seed;
    curandState *dev_states;
    srand((unsigned int)time(NULL));

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    //Host do not have GPUs...
    if(my_rank == 0)
    {
        MPI_Finalize();
        return RV_SUCESS;
    }

    //Each node has 8 GPUs...
    cudaSetDevice((my_rank-1)%8);
    cudaMalloc((void **)&dev_states, NUMBER_OF_THREADS_TESTS * NUMBER_OF_BLOCKS_TESTS *
                sizeof(curandState));

    //Initializes the PRNG for each cuda thread.
    seed = my_rank * random(); //makes sure that the seed is different in each MPI process.
    setup_kernel<<<NUMBER_OF_THREADS_TESTS, NUMBER_OF_BLOCKS_TESTS>>>(dev_states, seed);

    //Tests the implementation of each algorithm in cuda.
    test_algorithms_on_gpu(dev_states, ALG_TYPE_FORRO);
    test_algorithms_on_gpu(dev_states, ALG_TYPE_CHACHA);
    test_algorithms_on_gpu(dev_states, ALG_TYPE_SALSA);

    cudaFree(dev_states);

    MPI_Finalize();
    return RV_SUCESS;
}