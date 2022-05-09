/*
* This is an example of functionalities that might be useful
* to the researcher trying to find differentials or 
* linear approximations to ChaCha, Salsa or Forro.
*/
#include "arx_cryptanalysis.cuh"

int my_rank, num_procs;

void coutinho_2022_chacha_linear_approximations()
{
    uint64_t N = 1;
    uint32_t mask[STATE_SIZE] = {0};
    int words[4] = {0,4,8,12};
    linear_approximation_t la = {0};

    N<<=36;
    la.input.subround = 6;
    la.correlation.expected = 0.006942;
    strcpy(la.correlation.paper,"[Coutinho 2021]");
    la.alg_type = ALG_TYPE_CHACHA;
    la.correlation.number_of_trials = N;
    lob_set_bit(&la.input, 3, 0);
    lob_set_bit(&la.input, 4, 0);
    expand_linear_equation(&la, 6);
    compute_differential_or_linear_correlation(&la, TYPE_LINEAR);
    if(la.correlation.is_significant & my_rank==0)
      la_print(NULL, la);

    //For each group of proof of lemma 11
    linear_approximation_t group_la;
    strcpy(group_la.correlation.paper,"[Coutinho 2022]");
    group_la.alg_type = ALG_TYPE_CHACHA;
    group_la.correlation.number_of_trials = N;
    double expected[4] = {0.000307, 0.000383, 0.000435, 0.0625};
    for(int g=0; g<4; g++)
    {
        memset(mask,0x00,STATE_SIZE*sizeof(uint32_t));
        for(int i=0; i<4; i++)
            mask[words[i]+g] = 0xFFFFFFFF;

        memset(&group_la.output,0x00,sizeof(list_of_bits_t));
        group_la.correlation.expected = expected[g];
        group_la.input.subround = la.output.subround;
        and_array(group_la.input.mask, mask, la.output.mask, STATE_SIZE);
        lob_compute_list_of_bits_from_mask(&(group_la.input));
        expand_linear_equation(&group_la, 1);

        memcpy(&group_la.input,&group_la.output,sizeof(list_of_bits_t));
        memset(&group_la.output,0x00,sizeof(list_of_bits_t));

        if(g == 1) //in this case we use a extra rule to improve correlation
        {
            lob_set_bit(&group_la.input, 1,11);
            lob_set_bit(&group_la.input, 5,23);
            lob_set_bit(&group_la.input, 9,11);
        }
        if(g == 2) //in this case we use a extra rule to improve correlation
        {
            lob_set_bit(&group_la.input, 2,6);
            lob_set_bit(&group_la.input, 2,18);
            lob_set_bit(&group_la.input, 6,18);
            lob_set_bit(&group_la.input, 6,30);
            lob_set_bit(&group_la.input, 10,6);
            lob_set_bit(&group_la.input, 10,18);
        }
        expand_linear_equation(&group_la, 1);

        group_la.input.subround = la.output.subround;
        and_array(group_la.input.mask, mask, la.output.mask, STATE_SIZE);
        lob_compute_list_of_bits_from_mask(&(group_la.input));

        compute_differential_or_linear_correlation(&group_la, TYPE_LINEAR);

        if(la.correlation.is_significant & my_rank==0)
            la_print(NULL, group_la);
    }
}

void example_differential_correlation()
{
    differential_t diff = {
            ALG_TYPE_CHACHA,
            {{0}, {14},{6},0, 1}, //id
            {{0}, {3,4},{0,0},6,2}, //od
            {"Sec 3.1 of [Coutinho 20]",0.00048, 0, 0, 0}
        };

    lob_compute_mask_from_list_of_bits(&(diff.input));
    lob_compute_mask_from_list_of_bits(&(diff.output));
    
    diff.correlation.number_of_trials = 1;
    diff.correlation.number_of_trials <<= 35;
    compute_differential_or_linear_correlation(&diff, TYPE_DIFFERENTIAL);
    if(my_rank == 0)
        differential_print(NULL, diff);
}

int main()
{
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    example_differential_correlation();
    coutinho_2022_chacha_linear_approximations();

    //The following code was used to find all single bit differentials for 4 rounds of Salsa in order to apply BLE
    uint64_t number_of_trials = 1;
    number_of_trials<<=28;
    compute_all_single_bit_differential_correlation(ALG_TYPE_SALSA, 4, 0, 
        number_of_trials, "all_single_bit_differentials_for_4_rounds_of_salsa.dat");

    MPI_Finalize();
    return 0;
}