#include "arx_cryptanalysis.cuh"
#include "papers.cuh"
#include "util.cuh"

int my_rank, num_procs;

void forro_4round_linear_approximation(FILE *output_file)
{
    uint64_t N = 1;
    linear_approximation_t la = {0};

    N<<=38;
    la.correlation.expected = 0.0476;
    strcpy(la.correlation.paper,"[Coutinho 2022]");
    la.alg_type = ALG_TYPE_FORRO;
    la.correlation.number_of_trials = N;
    lob_define_single_bit(&la.input, 10,0,8);
    
    expand_linear_equation(&la, 8);

    compute_differential_or_linear_correlation(&la, TYPE_LINEAR);
    if(la.correlation.is_significant & my_rank==0)
      la_print(output_file, la);
}

void forro_5round_linear_approximation(FILE *output_file)
{
    uint64_t N = 1;
    linear_approximation_t L[15] = {0};
    double expected[5] = {0.0278, 0.1667, 0.0046, 0.0046, 0.000284};

    N<<=38;
    L[0].input.subround = 8;
    set_bit(L[0].input.mask, 10, 0);
    
    for(int i=0;i<13;i++)
    {
        L[i].alg_type = ALG_TYPE_FORRO;
        strcpy(L[i].correlation.paper,"[Coutinho 2022]");
        expand_linear_equation(&L[i],1);
        if(i>=8)
        {
            L[i].correlation.number_of_trials = N;
            L[i].correlation.expected = expected[i-8];
            compute_differential_or_linear_correlation(&L[i], TYPE_LINEAR);
            if(L[i].correlation.is_significant & my_rank == 0)
                la_print(output_file, L[i]);
        }
        memcpy(L[i+1].input.mask, L[i].output.mask, sizeof(uint32_t)*STATE_SIZE);
        L[i+1].input.subround = L[i].output.subround;
    }
}


void salsa_pnb_attack_using_ble(FILE *output_file)
{
    int differential_E1_part_subrounds = 1, differential_E2_part_subrounds = 3;
    int k = 0, alg_type = ALG_TYPE_SALSA;
    int idw = 7;
    int idb = 31;
    int odw = 4;
    int odb = 7; 
    pnb_t pnb = {0};
    difference_t od = {{0}, {0,4,12},{0,7,0},4, 3};
    lob_compute_mask_from_list_of_bits(&od);
    algorithm alg;
    define_alg(&alg, alg_type);

    pnb.subrounds = 7;
    pnb.threshold = 0.3;
    pnb.statistic_type = STAT_MEDIAN;
    pnb.correlation_of_g.number_of_trials = 1;
    pnb.correlation_of_g.number_of_trials <<= 34;
    pnb_define_alg(&pnb, alg_type);
    
    //Use Lipmaa and Moriai to compute first round differential
    lob_define_single_bit(&pnb.diff.input, idw, idb, 0);
    alg.differential_update(pnb.diff.input.mask, differential_E1_part_subrounds, &k);
    lob_compute_list_of_bits_from_mask(&pnb.diff.input);
    pnb.diff.input.subround = differential_E1_part_subrounds;

    //Compute differential correlation for each bit from the backward expansion
    double differential_correlation = 1;
    for(int i=0;i<od.number_of_bits;i++)
    {
        memset(&pnb.diff.output, 0, sizeof(difference_t));
        lob_define_single_bit(&pnb.diff.output, od.words[i], od.bits[i], 
            differential_E2_part_subrounds + differential_E1_part_subrounds);
        search_until_find_correlation(&pnb.diff, TYPE_DIFFERENTIAL);
        differential_correlation *= pnb.diff.correlation.observed;
    }
    pnb.diff.correlation.observed = differential_correlation;
    memcpy(&pnb.diff.output, &od, sizeof(difference_t));

    //Set linear expansion with correlation 1
    memcpy(&pnb.la.input, &pnb.diff.output, sizeof(diff_lin_t));
    pnb.la.correlation.observed = 1;
    lob_define_single_bit(&pnb.la.output, odw, odb,pnb.la.input.subround+1);

    //Attack 7 and 8 rounds
    for(int sr=7; sr<=8; sr++)
    {
        pnb.subrounds = sr;
        compute_neutrality_vector(&pnb, (1<<26));
        get_pnb_list(&pnb);
        compute_correlation_of_g(&pnb);
        compute_complexity_of_the_attack(&pnb);
        pnb.data_complexity += k;
        pnb.time_complexity += k;
        if(my_rank == 0)
            pnb_print(output_file, pnb);
    }
}

void euro2022(FILE *output_file)
{
    pnb_t pnb = {0};
    algorithm alg;
    int diff_E1_subrounds = 2, diff_E2_subrounds = 5, linear_part_subrounds = 1;
    int k=0, alg_type = ALG_TYPE_CHACHA;

    int idw = 13, idb = 6, odw = 2, odb = 0;
    pnb.subrounds = 14;
    pnb.threshold = 0.2;
    pnb.alg_type = alg_type;
    pnb.statistic_type = STAT_MEDIAN;
    uint64_t number_of_trials_for_neutrality = (1<<28), number_of_trials_for_bias_of_g = 1;
    number_of_trials_for_bias_of_g<<=31;

    //Use Lipmaa and Moriai to compute first round differential
    define_alg(&alg, alg_type);
    lob_define_single_bit(&pnb.diff.input, idw, idb, 0);
    alg.differential_update(pnb.diff.input.mask, diff_E1_subrounds, &k);
    lob_compute_list_of_bits_from_mask(&pnb.diff.input);
    pnb.diff.input.subround = diff_E1_subrounds;

    lob_define_single_bit(&pnb.diff.output, odw, odb, diff_E2_subrounds + diff_E1_subrounds);
    search_until_find_correlation(&pnb.diff, TYPE_DIFFERENTIAL);
    la_compute_from_differential(&pnb.la, pnb.diff, linear_part_subrounds);

    compute_neutrality_vector(&pnb, number_of_trials_for_neutrality);
    get_pnb_list(&pnb);

    pnb_remove(&pnb, 49); pnb_remove(&pnb, 51); pnb_remove(&pnb, 52);
    pnb_remove(&pnb, 59); pnb_remove(&pnb, 169);

    pnb_iteractive_selection(&pnb, 0.45, number_of_trials_for_bias_of_g, 12);
    pnb.correlation_of_g.number_of_trials = number_of_trials_for_bias_of_g;

    compute_correlation_of_g(&pnb);
    compute_complexity_of_the_attack(&pnb);

    if(my_rank == 0)
        pnb_print(output_file, pnb);
}


void differential_results(FILE *output_file)
{
    int count = 0;
    while(1)
    {
        differential_t diff = paperdiff[count];

        if(strcmp(diff.correlation.paper, "Stop")==0)
            break;

        lob_compute_mask_from_list_of_bits(&(diff.input));
        lob_compute_mask_from_list_of_bits(&(diff.output));
        search_until_find_correlation(&diff, TYPE_DIFFERENTIAL);
        if(my_rank == 0)
            differential_print(output_file, diff);

        count++;
    }
}

void linear_results(FILE *output_file)
{
    int count = 0;
    while(1)
    {
        linear_approximation_t lin_approx = paperlin[count];

        if(strcmp(lin_approx.correlation.paper, "Stop")==0)
            break;

        lob_compute_mask_from_list_of_bits(&(lin_approx.input));
        lob_compute_mask_from_list_of_bits(&(lin_approx.output));
        search_until_find_correlation(&lin_approx, TYPE_LINEAR);
        if(my_rank == 0)
            la_print(output_file, lin_approx);

        count++;
    }

#ifdef COUTINHO_2022_FORRO_LINEAR_APPROXIMATIONS
    forro_4round_linear_approximation(output_file);
    forro_5round_linear_approximation(output_file);
#endif
}

void pnb_results(FILE *output_file)
{
    uint64_t number_of_trials_for_neutrality = (1<<30), number_of_trials_for_neutrality_for_bias_of_g = 1;
    number_of_trials_for_neutrality_for_bias_of_g <<= 34;

#ifdef AUMASSON_2008_PNB
    //TABLE 2
    pnb_attack_for_single_bit_differential(7,31,1,14,7,4,0,0.9, number_of_trials_for_neutrality,
        number_of_trials_for_neutrality_for_bias_of_g, ALG_TYPE_SALSA, STAT_MEDIAN,output_file);
    pnb_attack_for_single_bit_differential(7,31,1,14,7,4,0,0.8, number_of_trials_for_neutrality,
        number_of_trials_for_neutrality_for_bias_of_g, ALG_TYPE_SALSA, STAT_MEDIAN, output_file);
    pnb_attack_for_single_bit_differential(7,31,1,14,7,4,0,0.7, number_of_trials_for_neutrality,
        number_of_trials_for_neutrality_for_bias_of_g, ALG_TYPE_SALSA, STAT_MEDIAN, output_file);
    pnb_attack_for_single_bit_differential(7,31,1,14,7,4,0,0.6, number_of_trials_for_neutrality,
        number_of_trials_for_neutrality_for_bias_of_g, ALG_TYPE_SALSA, STAT_MEDIAN, output_file);
    pnb_attack_for_single_bit_differential(7,31,1,14,7,4,0,0.5, number_of_trials_for_neutrality,
        number_of_trials_for_neutrality_for_bias_of_g, ALG_TYPE_SALSA, STAT_MEDIAN, output_file);

    //Attacks sec 3.5
    pnb_attack_for_single_bit_differential(7,31,1,14,8,4,0,0.12, number_of_trials_for_neutrality,
        number_of_trials_for_neutrality_for_bias_of_g, ALG_TYPE_SALSA, STAT_MEDIAN, output_file);
    pnb_attack_for_single_bit_differential(13,13,11,0,12,6,0,0.6, number_of_trials_for_neutrality,
        number_of_trials_for_neutrality_for_bias_of_g, ALG_TYPE_CHACHA, STAT_MEDIAN, output_file);
    pnb_attack_for_single_bit_differential(13,13,11,0,14,6,0,0.5, number_of_trials_for_neutrality,
        number_of_trials_for_neutrality_for_bias_of_g, ALG_TYPE_CHACHA, STAT_MEDIAN, output_file);
#endif

#ifdef CHOUDHURI_2016_PNB
    pnb_attack_for_single_bit_differential(12,21,2,0,12,6,2,0.4, number_of_trials_for_neutrality,
        number_of_trials_for_neutrality_for_bias_of_g, ALG_TYPE_CHACHA, STAT_MEDIAN, output_file);
    pnb_attack_for_single_bit_differential(12,21,2,0,12,6,0,0.4, number_of_trials_for_neutrality,
        number_of_trials_for_neutrality_for_bias_of_g, ALG_TYPE_CHACHA, STAT_MEDIAN, output_file);
#endif

#ifdef COUTINHO_2020_PNB
    pnb_attack_for_single_bit_differential(14,6,3,0,12,7,1,0.4, number_of_trials_for_neutrality,
        number_of_trials_for_neutrality_for_bias_of_g, ALG_TYPE_CHACHA, STAT_MEDIAN, output_file);
    pnb_attack_for_single_bit_differential(14,6,3,0,14,7,1,0.35, number_of_trials_for_neutrality,
        number_of_trials_for_neutrality_for_bias_of_g, ALG_TYPE_CHACHA, STAT_MEDIAN, output_file);
#endif
#ifdef COUTINHO_2022_FORRO_PNB
    pnb_attack_for_single_bit_differential(5,11,10,0,20,8,0,0.25, number_of_trials_for_neutrality,
        number_of_trials_for_neutrality_for_bias_of_g, ALG_TYPE_FORRO, STAT_MEDIAN, output_file);
#endif
#ifdef COUTINHO_2022_SALSA_PNB
    salsa_pnb_attack_using_ble(output_file);
#endif
#ifdef DEY_2022_CHACHA_PNB
    euro2022(output_file);
#endif
}

int main()
{
    FILE *p = NULL;
    
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    if(my_rank == 0) {
        create_folder_if_doesnt_exist("results/");
        p = fopen("results/paper_results.dat", "w");
    }
    differential_results(p);
    linear_results(p);
    pnb_results(p);

    if(my_rank == 0)
        fclose(p);

    MPI_Finalize();
    return 0;
}
