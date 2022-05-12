#include "arx_cryptanalysis.cuh"

void rule_word_not_from_subround(algorithm alg, uint32_t bitsRemaining[16], linear_approximation_t *la)
{
    int subround = la->input.subround + 1;
    for(int w=0; w<16; w++)
    {
        for(int bit=0; bit<32; bit++)
        {
            if(get_bit_from_word_and_bit(bitsRemaining, w, bit))
            {
                if(alg.get_letter(w, subround) == -1)
                {
                    set_bit(la->output.mask, w, bit);
                    set_bit(bitsRemaining, w, bit);
                }
            }
        }
    }
}

void rule_expand_linear_terms(algorithm alg, uint32_t bitsRemaining[16], linear_approximation_t *la)
{
    int bit = 0;
    for(int w=0; w<16; w++)
    {
        if(get_bit_from_word_and_bit(bitsRemaining, w, bit))
        {
            alg.expand_bit(la, w, bit, EXPANSION_LINEAR);
            set_bit(bitsRemaining, w, bit);
        }
    }
}


void rule_expand_adjacent_terms(algorithm alg, uint32_t bitsRemaining[16], linear_approximation_t *la)
{
    for(int w=0; w<16; w++)
    {
        for(int bit=0; bit<32; bit++)
        {
            if(get_bit_from_word_and_bit(bitsRemaining, w, bit))
            {
                if(get_bit_from_word_and_bit(bitsRemaining, w, bit+1))
                {
                    int letter = alg.get_letter(w, la->input.subround);
                    alg.expand_bit(la, w, bit, EXPANSION_LINEAR);
                    set_bit(bitsRemaining, w, bit);
                    alg.expand_bit(la, w, bit+1, EXPANSION_LINEAR);
                    set_bit(bitsRemaining, w, bit+1);
                }
            }
        }
    }
}


void rule_expand_nonlinear_terms(algorithm alg, uint32_t bitsRemaining[16], linear_approximation_t *la)
{   
    for(int w=0; w<16; w++)
    {
        for(int bit=0; bit<32; bit++)
        {
            if(get_bit_from_word_and_bit(bitsRemaining, w, bit))
            {
                alg.expand_bit(la, w, bit, EXPANSION_NON_LINEAR);
                set_bit(bitsRemaining, w, bit);
            }
        }
    }
}


int get_position_of_letter(algorithm alg, int letter, int subround)
{
    for(int w=0; w<16; w++)
    {
        if(letter == alg.get_letter(w, subround))
            return w;
    }

    return -1;
}

void special_rule_forro_expand_EDD(algorithm alg, uint32_t bitsRemaining[16], linear_approximation_t *la)
{
    int subround = la->input.subround + 1;
    int wD = get_position_of_letter(alg, LetterD, subround);
    int wE = get_position_of_letter(alg, LetterE, subround);

    for(int bit=0; bit<32; bit++)
    {
        if(get_bit_from_word_and_bit(bitsRemaining, wE, bit) && get_bit_from_word_and_bit(bitsRemaining, wD, bit))
        {
            if(get_bit_from_word_and_bit(bitsRemaining, wD, bit-1))
            {
                alg.expand_bit(la, wE, bit, EXPANSION_LINEAR);
                set_bit(bitsRemaining, wE, bit);
                alg.expand_bit(la, wD, bit, EXPANSION_LINEAR);
                set_bit(bitsRemaining, wD, bit);
                alg.expand_bit(la, wD, bit-1, EXPANSION_LINEAR);
                set_bit(bitsRemaining, wD, bit-1);
                alg.special_expansion_cases(la, bit-1, SPECIAL_EXPANSION_EDD);
            }
            else if(get_bit_from_word_and_bit(bitsRemaining, wD, bit+1))
            {
                alg.expand_bit(la, wE, bit, EXPANSION_LINEAR);
                set_bit(bitsRemaining, wE, bit);
                alg.expand_bit(la, wD, bit, EXPANSION_LINEAR);
                set_bit(bitsRemaining, wD, bit);
                alg.expand_bit(la, wD, bit+1, EXPANSION_LINEAR);
                set_bit(bitsRemaining, wD, bit+1);
                alg.special_expansion_cases(la, bit+1, SPECIAL_EXPANSION_EDD);
            }
        }
    }
}


void special_rule_chacha_expand_AAB(algorithm alg, uint32_t bitsRemaining[16], linear_approximation_t *la)
{
    int pattern[4][4][4] = { //[subrounds][number of words used] 
        {
            {0,4,8,12},
            {1,5,9,13},
            {2,6,10,14},
            {3,7,11,15},
        },
        {
            {0,4,8,12},
            {1,5,9,13},
            {2,6,10,14},
            {3,7,11,15},
        },
        {
            {0,5,10,15},
            {1,6,11,12},
            {2,7,8,13},
            {3,4,9,14}
        },
        {
            {0,5,10,15},
            {1,6,11,12},
            {2,7,8,13},
            {3,4,9,14}
        }
    };
    int subround = la->input.subround + 1;

    for(int wA=0; wA<4; wA++)
    {
        for(int bit=1; bit<32; bit++)
        {
            int wB = pattern[subround%4][wA][1];

            if(get_bit_from_word_and_bit(bitsRemaining, wA, bit) && 
                get_bit_from_word_and_bit(bitsRemaining, wA, bit-1) &&
                get_bit_from_word_and_bit(bitsRemaining, wB, bit))
                {
                    set_bit(bitsRemaining, wA, bit);
                    set_bit(bitsRemaining, wA, bit-1);
                    set_bit(bitsRemaining, wB, bit);
                    set_bit(la->output.mask, wA, bit);
                }
        }
    }
}

void expand_linear_equation_for_single_subround(linear_approximation_t *la)
{
    algorithm alg;
    uint32_t bitsRemaining[16];

    define_alg(&alg,la->alg_type);
    
    la->output.subround = la->input.subround + 1;
    memcpy(bitsRemaining, la->input.mask, sizeof(uint32_t)* 16);
    
    //NOTICE: rules must be in order, otherwise it will fail!!!!!!
    rule_word_not_from_subround(alg, bitsRemaining, la);
    //if(la->alg_type == ALG_TYPE_CHACHA)
    //    special_rule_chacha_expand_AAB(alg, bitsRemaining, la);
    rule_expand_linear_terms(alg, bitsRemaining, la);
    if(la->alg_type == ALG_TYPE_FORRO)
        special_rule_forro_expand_EDD(alg, bitsRemaining, la);
    rule_expand_adjacent_terms(alg, bitsRemaining, la);
    rule_expand_nonlinear_terms(alg, bitsRemaining, la);
}


void expand_linear_equation(linear_approximation_t *la, int subrounds)
{
    linear_approximation_t aux_la = {0};

    lob_compute_list_of_bits_from_mask(&(la->input));
    memcpy(&aux_la, la, sizeof(linear_approximation_t));
    for(int i=0; i<subrounds; i++)
    {
        expand_linear_equation_for_single_subround(&aux_la);
        lob_compute_list_of_bits_from_mask(&(aux_la.output));
        memcpy(&aux_la.input, &aux_la.output, sizeof(list_of_bits_t));
        memset(&aux_la.output, 0x00, sizeof(list_of_bits_t));
    }

    memcpy(&(la->output), &aux_la.input, sizeof(list_of_bits_t));
}