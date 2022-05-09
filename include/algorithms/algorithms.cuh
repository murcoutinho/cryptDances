#ifndef ALGO_H
#define ALGO_H

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <string.h>
#include <curand.h>
#include <curand_kernel.h>
#include <mpi.h>
#include "util.cuh"
#include "gputimer.cuh"
#include "cryptanalysis_types.cuh"

#define U32C(v) (v##U)
#define U32V(v) ((uint32_t)(v) &U32C(0xFFFFFFFF))

#define ROTATE(v, c) ( (c && (32-c)) ? ((v<<c)^(v>>(32-c))) : v )
#define XOR(v, w) ((v) ^ (w))
#define PLUS(v, w) (((v) + (w)))
#define MINUS(v, w) (((v) - (w)))
#define PLUSONE(v) (PLUS((v), 1))
#define U32TO8_LITTLE(p,v) (((uint32_t *) (p))[0]=v)
#define U8TO32_LITTLE(p) (((uint32_t *) (p))[0])

typedef struct {
    int alg_type;
    uint32_t key_positions[8];
    uint32_t iv_positions[4];
    int number_of_rounds;
    int number_of_subrounds_in_one_round;
    void(*init)(uint32_t *, uint32_t *, uint32_t *, uint32_t *);
    void(*rounds)(uint32_t *, uint32_t, uint32_t);
    void(*subrounds)(uint32_t *, uint32_t, uint32_t);
    void(*invert_rounds)(uint32_t *, uint32_t, uint32_t);
    void(*invert_subrounds)(uint32_t *, uint32_t, uint32_t);
    void(*encrypt_rounds)(uint32_t *, uint32_t *, uint32_t);
    void(*encrypt_subrounds)(uint32_t *, uint32_t *, uint32_t);
    void(*decrypt_rounds)(uint32_t *, uint32_t *, uint32_t *, uint32_t, uint32_t);
    void(*decrypt_subrounds)(uint32_t *, uint32_t *, uint32_t *, uint32_t, uint32_t);
    void(*expand_bit)(linear_approximation_t *, int, int, int);
    void(*special_expansion_cases)(linear_approximation_t *, int, int);
    int(*get_letter)(int, int);
    void(*differential_update)(uint32_t *, int, int *);
    char name[10];
} algorithm;

#define ALG_TYPE_CHACHA 0
#define ALG_TYPE_SALSA 1
#define ALG_TYPE_FORRO 2

__host__ __device__ void define_alg(algorithm *alg, uint32_t type);
void get_alg_name(char name[10], int alg_type);
void get_alg_iv_positions(uint32_t iv_positions[4], int alg_type);
#endif