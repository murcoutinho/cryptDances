#include "arx_cryptanalysis.cuh"
#include "papers.cuh"

int my_rank, num_procs;

void forro_4round_linear_approximation(FILE *output_file)
{
	uint64_t N = 1;
	linear_approximation_t la = {0};

	N<<=38;
	la.input.subround = 8;
	la.correlation.expected = 0.0476;
	strcpy(la.correlation.paper,"[Coutinho 2022]");
	la.alg_type = ALG_TYPE_FORRO;
	la.correlation.number_of_trials = N;
	set_bit(la.input.mask, 10, 0);

	lob_compute_list_of_bits_from_mask(&(la.input));
	
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
#ifdef AUMASSON_2008_PNB
	//TABLE 2
	pnb_attack_for_single_bit_differential(7,31,1,14,7,4,0,0.9, ALG_TYPE_SALSA, output_file);
	pnb_attack_for_single_bit_differential(7,31,1,14,7,4,0,0.8, ALG_TYPE_SALSA, output_file);
	pnb_attack_for_single_bit_differential(7,31,1,14,7,4,0,0.7, ALG_TYPE_SALSA, output_file);
	pnb_attack_for_single_bit_differential(7,31,1,14,7,4,0,0.6, ALG_TYPE_SALSA, output_file);
	pnb_attack_for_single_bit_differential(7,31,1,14,7,4,0,0.5, ALG_TYPE_SALSA, output_file);

	//Attacks sec 3.5
	pnb_attack_for_single_bit_differential(7,31,1,14,8,4,0,0.12, ALG_TYPE_SALSA, output_file);
	pnb_attack_for_single_bit_differential(13,13,11,0,12,6,0,0.6, ALG_TYPE_CHACHA, output_file);
	pnb_attack_for_single_bit_differential(13,13,11,0,14,6,0,0.5, ALG_TYPE_CHACHA, output_file);
#endif

#ifdef CHOUDHURI_2016_PNB
	pnb_attack_for_single_bit_differential(12,21,2,0,12,6,2,0.4, ALG_TYPE_CHACHA, output_file);
	pnb_attack_for_single_bit_differential(12,21,2,0,12,6,0,0.4, ALG_TYPE_CHACHA, output_file);
#endif

#ifdef COUTINHO_2020_PNB
	pnb_attack_for_single_bit_differential(14,6,3,0,12,7,1,0.4, ALG_TYPE_CHACHA, output_file);
	pnb_attack_for_single_bit_differential(14,6,3,0,14,7,1,0.35, ALG_TYPE_CHACHA, output_file);
#endif
#ifdef COUTINHO_2022_PNB
	pnb_attack_for_single_bit_differential(5,11,10,0,20,8,0,0.25, ALG_TYPE_FORRO, output_file);
#endif
}

int main()
{
	FILE *p = NULL;
	
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	if(my_rank == 0)
		p = fopen("results/paper_results.dat", "w");

	differential_results(p);
	linear_results(p);
	pnb_results(p);

	if(my_rank == 0)
		fclose(p);

	MPI_Finalize();
	return 0;
}