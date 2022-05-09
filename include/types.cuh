#ifndef TYPESCUH
#define TYPESCUH

#include <stdint.h>
#include "constants.cuh"

#define R1 10
#define R2 27
#define R3 8

#define NUMBER_OF_SPECIAL_CASES 1

#define PRINT_DO_NOT 0
#define PRINT_SIGNIFICANT 1
#define PRINT_ALWAYS 2

#define LetterA 0
#define LetterB 1
#define LetterC 2
#define LetterD 3
#define LetterE 4

#define EXPANSION_LINEAR 0
#define EXPANSION_NON_LINEAR 1

#define SPECIAL_EXPANSION_EDD 1

typedef struct expansion
{
    uint32_t list_of_words[10];
    uint32_t list_of_bits[10];
    uint32_t number_of_terms;
    uint32_t number_of_linear_terms;
    uint32_t special_cases_bits[NUMBER_OF_SPECIAL_CASES][2];
} expansion;

#endif