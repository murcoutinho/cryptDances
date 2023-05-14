#ifndef CONSTANTS_CUH
#define CONSTANTS_CUH

#define RV_SUCESS 0
#define RV_ERROR 1

#define TRUE 1
#define FALSE 0

//CUDA
#define NUMBER_OF_DEVICES_PER_MACHINE 1
#define NUMBER_OF_TESTS_PER_THREAD (1<<12)//(1<<15) 
#define NUMBER_OF_CUDA_THREADS (1<<6) //(1<<8)
#define NUMBER_OF_CUDA_BLOCKS (1<<5) //Must be >= (1<<7)!!!
#define TOTAL_EXECUTIONS_PER_KERNEL (NUMBER_OF_TESTS_PER_THREAD * NUMBER_OF_CUDA_THREADS * NUMBER_OF_CUDA_BLOCKS)

#define MAX_LEVEL 38 //We call LEVEL the log2(N), where N is the number of trials when computing a correlation.
#define MIN_LEVEL 34 //We call LEVEL the log2(N), where N is the number of trials when computing a correlation.
#define MAX_BITS_IN_LIST_OF_BITS 128

#define MAXIMUM_STATE_SIZE 16 //Maximum state size
#define MAXIMUM_STATE_SIZE_IN_BYTES (MAXIMUM_STATE_SIZE * sizeof(uint32_t))
#define MAXIMUM_STATE_SIZE_IN_BITS (MAXIMUM_STATE_SIZE_IN_BYTES * 8 )
#define KEY_SIZE 8 //Maximum key size
#define KEY_SIZE_IN_BITS (256)
#define NONCE_SIZE 2
#define CTR_SIZE 2
#define MAXIMUM_NUMBER_OF_POSSIBLE_SINGLE_BIT_DIFFERENTIALS 65536

#endif