#ifndef TESTS_CUH
#define TESTS_CUH

#include "arx_cryptanalysis.cuh"

#define NUMBER_OF_THREADS_TESTS 512
#define NUMBER_OF_BLOCKS_TESTS 32
#define TOTAL_THREADS_TESTS (NUMBER_OF_THREADS_TESTS * NUMBER_OF_BLOCKS_TESTS)

__global__ void ker_test_algorithms(int *rv);
int test_algorithms_on_gpu(curandState *dev_states, int alg_type);

#endif

