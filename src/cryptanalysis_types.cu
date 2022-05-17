#include "arx_cryptanalysis.cuh"

//Prefix dt stands for "differential type".
void lob_compute_mask_from_list_of_bits(list_of_bits_t *lob)
{
    memset(lob->mask, 0x00, sizeof(uint32_t) * STATE_SIZE);
    for(int i=0; i<lob->number_of_bits; i++)
        lob->mask[lob->words[i]] |= 1<<(lob->bits[i]);
}

int lob_compute_list_of_bits_from_mask(list_of_bits_t *lob)
{
    lob->number_of_bits = 0;
    
    for (int w = 0; w < STATE_SIZE; w++)
    {
        for (int b = 0; b < 32; b++)
        {
            if (get_bit_from_word_and_bit(lob->mask, w, b))
            {
                if(lob->number_of_bits == MAX_BITS_IN_LIST_OF_BITS)
                {
                    memset(lob->mask, 0x00, sizeof(uint32_t) * STATE_SIZE);
                    lob->number_of_bits = 0;
                    return RV_ERROR;
                }
                lob->words[lob->number_of_bits] = w;
                lob->bits[lob->number_of_bits] = b;
                lob->number_of_bits++;
            }
        }
    }
    return RV_SUCESS;
}

void lob_print(FILE *p, list_of_bits_t lob)
{
    int current_word = 0;

    dprintf(p, "x^{[%d]}_{%d}[%d", lob.subround, lob.words[0], lob.bits[0]);
    current_word = lob.words[0];

    for(int i=1; i<lob.number_of_bits; i++)
    {
        if(lob.words[i] == current_word)
        {
            dprintf(p, ",%d", lob.bits[i]);
        }
        else
        {
            dprintf(p, "] + x^{[%d]}_{%d}[%d", lob.subround, lob.words[i], lob.bits[i]);
            current_word = lob.words[i];
        }
    }
    dprintf(p, "]\n");
}

void lob_set_bit(list_of_bits_t *lob, int word, int bit)
{
    set_bit(lob->mask, word, bit);
    lob_compute_list_of_bits_from_mask(lob);
}

void lob_define_single_bit(list_of_bits_t *lob, int word, int bit, int subround)
{
    memset(lob, 0x00, sizeof(list_of_bits_t));
    lob->number_of_bits = 1;
    lob->subround = subround;
    lob_set_bit(lob, word, bit);
}

void differential_print(FILE *p, differential_t diff)
{
    char name[10];

    get_alg_name(name ,diff.alg_type);
    if(diff.correlation.paper[0] == 0)
    {
        dprintf(p, "Differential correlation.\n");
    }
    else
    {
        dprintf(p, "Differential correlation from %s.\n",diff.correlation.paper);
    }

    dprintf(p, "Algorithm = %s.\n",name);
    dprintf(p, "ID = ");
    lob_print(p, diff.input);
    dprintf(p, "OD = ");
    lob_print(p, diff.output);
    if(diff.correlation.expected)
        dprintf(p,"Expected corr = %.15f.\n",diff.correlation.expected);
    dprintf(p,"Correlation   = %.15f.\n",diff.correlation.observed);
    dprintf(p,"Number of trials = 2^{%d}.\n\n",(int)log2((double)diff.correlation.number_of_trials));
}

void la_print(FILE *p, linear_approximation_t lin_approx)
{
    char name[10];

    get_alg_name(name ,lin_approx.alg_type);
    if(lin_approx.correlation.paper[0] == 0)
    {
        dprintf(p, "Linear correlation.\n");
    }
    else
    {
        dprintf(p, "Linear correlation from %s.\n",lin_approx.correlation.paper);
    }
    dprintf(p, "Algorithm = %s.\n",name);
    dprintf(p, "Input  Mask = ");
    lob_print(p, lin_approx.input);
    dprintf(p, "Output Mask = ");
    lob_print(p, lin_approx.output);
    if(lin_approx.correlation.expected)
        dprintf(p,"Expected corr = %.15f.\n",lin_approx.correlation.expected);
    dprintf(p,"Correlation   = %.15f.\n",lin_approx.correlation.observed);
    if(lin_approx.correlation.number_of_trials == 0)
    {
        dprintf(p,"Number of trials = 0.\n\n");
    }
    else
    {
        dprintf(p,"Number of trials = 2^{%d}.\n\n",(int)log2((double)lin_approx.correlation.number_of_trials));
    }
}

void pnb_print(FILE *p, pnb_t pnb)
{
    char name[10];
    get_alg_name(name ,pnb.alg_type);
    dprintf(p, "--------------------------------------------\n");
    if(pnb.correlation_of_g.paper[0] == 0)
    {
        dprintf(p, "PNB attack.\n");
    }
    else
    {
        dprintf(p, "PNB attack from %s.\n",pnb.correlation_of_g.paper);
    }
    dprintf(p, "Algorithm = %s.\n",name);

    differential_print(p, pnb.diff);
    la_print(p,pnb.la);

    dprintf(p, "Neutrality measure:\n");
    for(int i=0;i<KEY_SIZE_IN_BITS;i++)
    {
        dprintf(p, "%f ", pnb.neutrality_measure[i]);
    }
    dprintf(p, "\n\nNumber of PNBs = %d.\n", pnb.number_of_pnb);
    dprintf(p, "List of PNBs:\n");
    for(int i=0;i<pnb.number_of_pnb;i++)
    {
        dprintf(p, "%d ", pnb.pnb[i]);
    }

    if(pnb.correlation_of_g.is_significant)
    {
        dprintf(p, "\n\nCorrelation e_a = %f.\n", pnb.correlation_of_g.observed);
        dprintf(p, "Data complexity = %f.\n", pnb.data_complexity);
        dprintf(p, "Time complexity = %f.\n", pnb.time_complexity);
    }
    else
    {
        dprintf(p, "\n\nCorrelation is not significant...\n")
    }
    
    dprintf(p, "--------------------------------------------\n");
}

void ct_compute_and_test_correlation(correlation_t *corr)
{
    corr->observed = ((double)corr->correlation_count) / corr->number_of_trials;  
    corr->observed = 2 * corr->observed - 1;
    corr->is_significant = test_significance_of_correlation(corr->observed, corr->number_of_trials);
}

void ct_compute_and_test_correlation_using_median(correlation_t *corr, unsigned long long int number_of_threads)
{
    corr->observed = ((double)corr->correlation_count) / number_of_threads;
    corr->observed = 2 * corr->observed - 1;
    corr->is_significant = test_significance_of_correlation(corr->observed, corr->number_of_trials);
}


void print_latex_linear_relation(linear_approximation_t *lin_approx)
{
    double prob;
    printf("$$ \\begin{array}{cl}\n");
    if (fabs(lin_approx->correlation.observed) > 0)
    {
        printf("%f = (1+%f)/2 = &\\\\ \\Pr(", (1 + lin_approx->correlation.observed) / 2, lin_approx->correlation.observed);
    }
    else
    {
        prob = ((double)lin_approx->correlation.number_of_trials - lin_approx->correlation.correlation_count) / lin_approx->correlation.number_of_trials;
        printf("%f = (1+%f)/2 = \\Pr(", prob, 2 * prob - 1);
    }
    for (int w = 0; w < STATE_SIZE; w++)
    {
        for (int b = 0; b < 32; b++)
        {
            if (get_bit_from_word_and_bit(lin_approx->input.mask, w, b))
                printf("x^{[%d]}_{%d,%d} \\oplus ", lin_approx->input.subround, w, b);
        }
    }
    printf(" = & ");

    int count = 0;
    for (int w = 0; w < STATE_SIZE; w++)
    {
        for (int b = 0; b < 32; b++)
        {
            if (get_bit_from_word_and_bit(lin_approx->output.mask, w, b))
            {
                printf("x^{[%d]}_{%d,%d} \\oplus ", lin_approx->output.subround, w, b);
                count++;
                if (count == 8)
                {
                    count = 0;
                    printf("\\\\ &");
                }
            }
        }
    }
    printf(") \\\\ & ");

    count = 0;
    int flag = 0;
    for (int w = 0; w < STATE_SIZE; w++)
    {
        count++;
        if (count == 4)
        {
            count = 0;
            printf("\\\\ &");
        }
        flag = 0;
        for (int b = 0; b < 32; b++)
        {
            if (get_bit_from_word_and_bit(lin_approx->output.mask, w, b))
            {
                if (flag)
                    printf(",");
                if (flag == 0)
                {
                    printf("x^{[%d]}_{%d}[", lin_approx->output.subround, w);
                    flag = 1;
                }
                printf("%d",  b);
            }
        }
        if(flag)
            printf("] \\oplus ");
    }

    printf("\\end{array} $$ \n\n");
}


void differential_compute_from_single_bit(differential_t *diff, int idw, int idb, int odw, int odb, int output_subround, int alg_type)
{
    memset(diff, 0x00, sizeof(differential_t));
    diff->alg_type = alg_type;
    lob_define_single_bit(&(diff->input), idw, idb, 0);
    lob_define_single_bit(&(diff->output), odw, odb, output_subround);
    search_until_find_correlation(diff, TYPE_DIFFERENTIAL);
}

void la_compute_from_differential(linear_approximation_t *la, differential_t diff, int subrounds)
{
    la->alg_type = diff.alg_type;
    memcpy(&(la->input), &diff.output, sizeof(list_of_bits_t));
    if(subrounds == 0)
    {
        memcpy(&(la->output), &diff.output, sizeof(list_of_bits_t));
        la->correlation.observed = 1;
        return;
    }
    expand_linear_equation(la, subrounds);
    search_until_find_correlation(la, TYPE_LINEAR);
}

void pnb_define_alg(pnb_t *pnb, int alg_type)
{
    pnb->alg_type = alg_type;
    pnb->diff.alg_type = alg_type;
    pnb->la.alg_type = alg_type;
}