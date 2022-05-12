#include "arx_cryptanalysis.cuh"
#include "forro.cuh"
#include "chacha.cuh"
#include "salsa.cuh"

__host__ __device__ void define_alg(algorithm *alg, uint32_t type)
{
    uint32_t forro_iv_positions[4] = { 4,5,12,13 };
    uint32_t forro_key_positions[8] = { 0,1,2,3,8,9,10,11 };

    uint32_t chacha_iv_positions[4] = { 12,13,14,15 };
    uint32_t chacha_key_positions[8] = { 4,5,6,7,8,9,10,11 };

    uint32_t salsa_iv_positions[4] = { 6,7,8,9 };
    uint32_t salsa_key_positions[8] = { 1,2,3,4,11,12,13,14 };

    switch (type)
    {
    case ALG_TYPE_FORRO:
        for (int i = 0; i < 4; i++)
            alg->iv_positions[i] = forro_iv_positions[i];
        for (int i = 0; i < 8; i++)
            alg->key_positions[i] = forro_key_positions[i];
        alg->alg_type = ALG_TYPE_FORRO;
        alg->number_of_rounds = FORRO_NUMBER_OF_ROUNDS;
        alg->number_of_subrounds_in_one_round = FORRO_NUMBER_OF_SUBROUNDS_IN_EACH_ROUND;
        alg->init = &forro_init;
        alg->rounds = &forro_rounds;
        alg->subrounds = &forro_subrounds;
        alg->invert_rounds = &forro_invert_rounds;
        alg->invert_subrounds = &forro_invert_subrounds;
        alg->encrypt_rounds = &forro_encrypt_rounds;
        alg->encrypt_subrounds = &forro_encrypt_subrounds;
        alg->decrypt_rounds = &forro_decrypt_rounds;
        alg->decrypt_subrounds = &forro_decrypt_subrounds;
        alg->expand_bit = &forro_expand_bit;
        alg->special_expansion_cases = &forro_special_cases;
        alg->get_letter = &forro_get_letter;
        alg->differential_update = NULL;
        alg->name[0] = 'F'; alg->name[1] = 'o'; alg->name[2] = 'r';
        alg->name[3] = 'r'; alg->name[4] = 'o'; alg->name[5] = 0;
        break;

    case ALG_TYPE_CHACHA:
        for (int i = 0; i < 4; i++)
            alg->iv_positions[i] = chacha_iv_positions[i];
        for (int i = 0; i < 8; i++)
            alg->key_positions[i] = chacha_key_positions[i];

        alg->alg_type = ALG_TYPE_CHACHA;
        alg->number_of_rounds = CHACHA_NUMBER_OF_ROUNDS;
        alg->number_of_subrounds_in_one_round = CHACHA_NUMBER_OF_SUBROUNDS_IN_EACH_ROUND;
        alg->init = &chacha_init;
        alg->rounds = &chacha_rounds;
        alg->subrounds = &chacha_subrounds;
        alg->invert_rounds = &chacha_invert_rounds;
        alg->invert_subrounds = &chacha_invert_subrounds;
        alg->encrypt_rounds = &chacha_encrypt_rounds;
        alg->encrypt_subrounds = &chacha_encrypt_subrounds;
        alg->decrypt_rounds = &chacha_decrypt_rounds;
        alg->decrypt_subrounds = &chacha_decrypt_subrounds;
        alg->expand_bit = &chacha_expand_bit;
        alg->special_expansion_cases = NULL;
        alg->get_letter = &chacha_get_letter;
        alg->differential_update = &chacha_differential_update;

        alg->name[0] = 'C'; alg->name[1] = 'h'; alg->name[2] = 'a';
        alg->name[3] = 'C'; alg->name[4] = 'h'; alg->name[5] = 'a';
        alg->name[6] = 0;
        break;

    case ALG_TYPE_SALSA:
        for (int i = 0; i < 4; i++)
            alg->iv_positions[i] = salsa_iv_positions[i];
        for (int i = 0; i < 8; i++)
            alg->key_positions[i] = salsa_key_positions[i];

        alg->alg_type = ALG_TYPE_SALSA;
        alg->number_of_rounds = SALSA_NUMBER_OF_ROUNDS;
        alg->number_of_subrounds_in_one_round = SALSA_NUMBER_OF_SUBROUNDS_IN_EACH_ROUND;
        alg->init = &salsa_init;
        alg->rounds = &salsa_rounds;
        alg->subrounds = &salsa_subrounds;
        alg->invert_rounds = &salsa_invert_rounds;
        alg->invert_subrounds = &salsa_invert_subrounds;
        alg->encrypt_rounds = &salsa_encrypt_rounds;
        alg->encrypt_subrounds = &salsa_encrypt_subrounds;
        alg->decrypt_rounds = &salsa_decrypt_rounds;
        alg->decrypt_subrounds = &salsa_decrypt_subrounds;
        alg->expand_bit = NULL;
        alg->special_expansion_cases = NULL;
        alg->get_letter = NULL;
        alg->differential_update = &salsa_differential_update;
        alg->name[0] = 'S'; alg->name[1] = 'a'; alg->name[2] = 'l';
        alg->name[3] = 's'; alg->name[4] = 'a'; alg->name[5] = 0;
        break;
    default:
        break;
    }
}

void get_alg_name(char name[10], int alg_type)
{
    algorithm alg;
    define_alg(&alg, alg_type);
    strcpy(name, alg.name);
}

void get_alg_iv_positions(uint32_t iv_positions[4], int alg_type)
{
    algorithm alg;
    define_alg(&alg, alg_type);
    memcpy(iv_positions, alg.iv_positions, 4*sizeof(uint32_t));
}


