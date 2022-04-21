#ifndef DLAH
#define DLAH

#include "types.cuh"
#include "cryptanalysis_types.cuh"

void expand_linear_equation(linear_approximation_t *la, int subrounds);

void expand_linear_equation_for_single_subround(linear_approximation_t *la);
#endif