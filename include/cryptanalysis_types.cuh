#ifndef DIFFERENTIAL_CUH
#define DIFFERENTIAL_CUH

#include "arx_cryptanalysis.cuh"

typedef struct{
    uint32_t mask[STATE_SIZE];
    int words[MAX_BITS_IN_LIST_OF_BITS];
    int bits[MAX_BITS_IN_LIST_OF_BITS];
    int subround;
    int number_of_bits;
} list_of_bits_t;

typedef list_of_bits_t difference_t;
typedef list_of_bits_t linear_mask_t;

typedef struct{
    char paper[40];
    double expected;
    double observed;
    uint64_t correlation_count;
    uint64_t number_of_trials;
    uint8_t is_significant;
} correlation_t;

typedef struct{
    int alg_type;
    list_of_bits_t input;
    list_of_bits_t output;
    correlation_t correlation;
} diff_lin_t;

typedef diff_lin_t differential_t;
typedef diff_lin_t linear_approximation_t;

typedef struct{
    int alg_type;
    int statistic_type;
    int subrounds;
    differential_t diff;
    linear_approximation_t la;
    correlation_t correlation_of_g;
    double threshold;
    double neutrality_measure[KEY_SIZE_IN_BITS];
    int pnb[KEY_SIZE_IN_BITS];
    int number_of_pnb;
    double time_complexity;
    double data_complexity;
} pnb_t;

void lob_compute_mask_from_list_of_bits(list_of_bits_t *lob);
int lob_compute_list_of_bits_from_mask(list_of_bits_t *lob);
void lob_print(FILE *p, list_of_bits_t lob);
void lob_set_bit(list_of_bits_t *lob, int word, int bit);
void lob_define_single_bit(list_of_bits_t *lob, int word, int bit, int subround);
void differential_print(FILE *p, differential_t diff);
void la_print(FILE *p, linear_approximation_t lin_approx);
void ct_compute_and_test_correlation(correlation_t *corr);
void ct_compute_and_test_correlation_using_median(correlation_t *corr, unsigned long long int number_of_threads);
void print_latex_linear_relation(linear_approximation_t *lin_approx);
void differential_compute_from_single_bit(differential_t *diff, int idw, int idb, int odw, int odb, int output_subround, int alg_type);
void la_compute_from_differential(linear_approximation_t *la, differential_t diff, int subrounds);
void pnb_print(FILE *p, pnb_t pnb);
void pnb_define_alg(pnb_t *pnb, int alg_type);

#endif