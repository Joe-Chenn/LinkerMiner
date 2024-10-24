// pch.cpp: 与预编译标头对应的源文件

#include "pch.h"
#include "malloc.h"
#include "stdio.h"
#include "math.h"

// 当使用预编译的头时，需要使用此源文件，编译才能成功。
double calculate_peptide_mass(char* p_sq, int sq_len, double input_mass, double* AA_mass)
{
    int i, j, k;
    double out_mass = input_mass;
    for (i = 0; i < sq_len; i++)
    {
        out_mass = out_mass + AA_mass[p_sq[i] - 'A'];
    }
    return out_mass;
}

double calculate_gdm(char* p_sq, int sq_len, int* p_mod_site_list, int fix_mod_num, double* p_gdm_matrix, int AA_num, int gdm_max_sq_len)
{
    int i, j, k;
    double gdm_value;

    // 一串是一个氨基酸的哥德尔编码值，位置索引是氨基酸在肽段的位置

    gdm_value = 0.0;

    if (p_mod_site_list[0] >= fix_mod_num)
    {
        gdm_value = gdm_value + p_gdm_matrix[(AA_num + p_mod_site_list[0] - fix_mod_num) * gdm_max_sq_len + 0];
    }
    else
    {
        gdm_value = gdm_value + 0.0;
    }

    for (i = 1; i < sq_len + 1; i++)
    {
        if (p_mod_site_list[i] >= fix_mod_num)
        {
            //printf("i = %d, j = %d\n", AA_num + p_mod_site_list[i] - fix_mod_num, i);
            gdm_value = gdm_value + p_gdm_matrix[(AA_num + p_mod_site_list[i] - fix_mod_num) * gdm_max_sq_len + i];
            //printf("%lf\n", p_gdm_matrix[(AA_num + p_mod_site_list[i] - fix_mod_num) * gdm_max_sq_len + i]);
        }
        else
        {
            //printf("i = %d, j = %d\n", p_sq[i - 1] - 'A', i);
            gdm_value = gdm_value + p_gdm_matrix[(p_sq[i - 1] - 'A') * gdm_max_sq_len + i];
            //printf("%lf\n", p_gdm_matrix[(p_sq[i - 1] - 'A') * gdm_max_sq_len + i]);
        }
    }

    if (p_mod_site_list[sq_len + 1] >= fix_mod_num)
    {
        gdm_value = gdm_value + p_gdm_matrix[(AA_num + p_mod_site_list[sq_len] - fix_mod_num) * gdm_max_sq_len + sq_len + 1];
    }
    else
    {
        gdm_value = gdm_value + 0.0;
    }
    return gdm_value;
}

long long tool_C_get_ion_code(long long pkl_file_index, long long ind_file_index, long long one_ind_num, long long link_site_combination_num)
{
    long long tmp_link_site_combination_num, tmp_one_ind_num, tmp_ind_file_index, tmp_pkl_file_index, out;
    tmp_link_site_combination_num = link_site_combination_num;
    tmp_one_ind_num = one_ind_num << 6;
    tmp_ind_file_index = ind_file_index << 50;
    tmp_pkl_file_index = pkl_file_index << 57;
    out = tmp_link_site_combination_num | tmp_one_ind_num | tmp_ind_file_index | tmp_pkl_file_index;
    return out;
}

long* tool_C_new_get_ion_data(long long input_data)
{
    long out[4];
    out[0] = input_data & 0x3f;
    out[1] = input_data >> 6 & 0x7ffffffffff;
    out[2] = input_data >> 50 & 0x7f;
    out[3] = input_data >> 57 & 0x7f;
    return out;
}

long long tool_C_peptide_info(int protein_index, int peptide_start, int peptide_end)
{
    return 0;
}

int tool_C_put_ion_to_ion_index(ION_INDEX_MZ* ion_index_start, char* pep_sq, int pep_len, int* mod_site_list, double b_start_mass, double y_start_mass, double* mod_mass_list, double* AA_mass_list, int link_site, int multi, int max_mass, long pep_index_num)
{
    int i, j, k;
    double b_mass = b_start_mass;
    double y_mass = y_start_mass;

    double b_mz;
    double y_mz;

    for (i = 0; i < link_site; i++)
    {
        if (mod_site_list[i] != -1)
        {
            b_mass += mod_mass_list[mod_site_list[i]];
        }
        if (i > 0 and i < pep_len + 1)
        {
            b_mass += AA_mass_list[pep_sq[i - 1] - 65];
            put_peptide_to_ion_index(ion_index_start, b_mass, multi, max_mass, pep_index_num);
        }
    }

    for (i = pep_len + 1; i > link_site; i--)
    {
        if (mod_site_list[i] != -1)
        {
            y_mass += mod_mass_list[mod_site_list[i]];
        }
        if (i > 0 and i < pep_len + 1)
        {
            y_mass += AA_mass_list[pep_sq[i - 1] - 65];
            put_peptide_to_ion_index(ion_index_start, y_mass, multi, max_mass, pep_index_num);
        }
    }
    return 0;
}

int put_peptide_to_ion_index(ION_INDEX_MZ* ion_index_start, double input_mz, int multi, int max_mz, long pep_index_num)
{
    int i, j, k;
    int mz;
    long* one_new_pep_start, * one_old_pep_start, * one_pep_start;
    mz = (int)(input_mz * multi + 0.5);
    if (mz > max_mz)
    {
    }
    else
    {
        if (ion_index_start[mz].pep_num == 0)
        {
            one_new_pep_start = (long*)malloc(sizeof(long) * 100);
            one_new_pep_start[ion_index_start[mz].pep_num] = pep_index_num;
            ion_index_start[mz].pep_num++;
            ion_index_start[mz].pep_size = 100;
            ion_index_start[mz].pep = one_new_pep_start;
        }
        else
        {
            if (ion_index_start[mz].pep_num >= ion_index_start[mz].pep_size)
            {
                // 需要重新建一个，搬移数据
                one_new_pep_start = (long*)malloc(sizeof(long) * ion_index_start[mz].pep_size * 10);
                one_old_pep_start = ion_index_start[mz].pep;
                for (i = 0; i < ion_index_start[mz].pep_num; i++)
                {
                    one_new_pep_start[i] = one_old_pep_start[i];
                }
                // 释放旧的内存
                free(one_old_pep_start);
                one_new_pep_start[ion_index_start[mz].pep_num] = pep_index_num;
                ion_index_start[mz].pep_num++;
                ion_index_start[mz].pep_size = ion_index_start[mz].pep_size * 10;
                ion_index_start[mz].pep = one_new_pep_start;
            }
            else
            {
                ion_index_start[mz].pep[ion_index_start[mz].pep_num] = pep_index_num;
                ion_index_start[mz].pep_num++;
            }
        }
    }
    return 0;
}

int tool_C_put_ion_to_ion_index_singlepeptide(ION_INDEX_MZ* ion_index_start, char* pep_sq, int pep_len, int* mod_site_list, double b_start_mass, double y_start_mass, double* mod_mass_list, double* AA_mass_list, int multi, int max_mass, long pep_index_num)
{
    int i, j, k;
    double b_mass = b_start_mass;
    double y_mass = y_start_mass;

    double b_mz;
    double y_mz;

    double add_aa_mass = 0;

    for (i = 0; i < pep_len; i++)
    {
        add_aa_mass = 0;
        if (mod_site_list[i] != -1)
        {
            add_aa_mass += mod_mass_list[mod_site_list[i]];
        }
        if (i > 0 and i < pep_len + 1)
        {
            add_aa_mass += AA_mass_list[pep_sq[i - 1] - 65];
        }
        if (fabs(add_aa_mass) > 1e-4)
        {
            b_mass += add_aa_mass;
            put_peptide_to_ion_index(ion_index_start, b_mass, multi, max_mass, pep_index_num);
        }
    }
    for (i = pep_len + 1; i > 0; i--)
    {
        add_aa_mass = 0;
        if (mod_site_list[i] != -1)
        {
            add_aa_mass += mod_mass_list[mod_site_list[i]];
        }
        if (i > 0 and i < pep_len + 1)
        {
            add_aa_mass += AA_mass_list[pep_sq[i - 1] - 65];
        }
        if (fabs(add_aa_mass) > 1e-4)
        {
            y_mass += add_aa_mass;
            put_peptide_to_ion_index(ion_index_start, y_mass, multi, max_mass, pep_index_num);
        }
    }
    return 0;
}

int tool_C_put_ion_to_ion_index_crosslink(ION_INDEX_MZ* ion_index_start, char* pep_sq, int pep_len, int* mod_site_list, double b_start_mass, double y_start_mass, double* mod_mass_list, double* AA_mass_list, int link_site, int multi, int max_mass, long pep_index_num)
{
    int i, j, k;
    double b_mass = b_start_mass;
    double y_mass = y_start_mass;

    double b_mz;
    double y_mz;

    double add_aa_mass = 0;

    for (i = 0; i < link_site; i++)
    {
        add_aa_mass = 0;
        if (mod_site_list[i] != -1)
        {
            add_aa_mass += mod_mass_list[mod_site_list[i]];
            if (i == 0)
            {
                b_mass += add_aa_mass;
            }
        }
        if (i > 0 and i < pep_len + 1)
        {
            add_aa_mass += AA_mass_list[pep_sq[i - 1] - 65];
        }
        if (fabs(add_aa_mass) > 1e-4)
        {
            b_mass += add_aa_mass;
            
            put_peptide_to_ion_index(ion_index_start, b_mass, multi, max_mass, pep_index_num);
        }
    }

    for (i = pep_len + 1; i > link_site; i--)
    {
        add_aa_mass = 0;
        if (mod_site_list[i] != -1)
        {
            add_aa_mass += mod_mass_list[mod_site_list[i]];
            if (i == pep_len + 1)
            {
                y_mass += add_aa_mass;
            }
        }
        if (i > 0 and i < pep_len + 1)
        {
            add_aa_mass += AA_mass_list[pep_sq[i - 1] - 65];
        }
        if (fabs(add_aa_mass) > 1e-4)
        {
            y_mass += add_aa_mass;
            put_peptide_to_ion_index(ion_index_start, y_mass, multi, max_mass, pep_index_num);
        }
    }
    return 0;
}

int tool_C_put_ion_to_ion_index_looplink(ION_INDEX_MZ* ion_index_start, char* pep_sq, int pep_len, int* mod_site_list, double b_start_mass, double y_start_mass, double* mod_mass_list, double* AA_mass_list, int link_site_1, int link_site_2, double loop_mass, int multi, int max_mass, long pep_index_num)
{
    // 要求link_site_1 小于link_site_1
    int i, j, k;
    double b_mass = b_start_mass;
    double y_mass = y_start_mass;

    double b_mz;
    double y_mz;

    double add_aa_mass = 0;

    for (i = 0; i < pep_len; i++)
    {
        add_aa_mass = 0;
        if (mod_site_list[i] != -1)
        {
            add_aa_mass += mod_mass_list[mod_site_list[i]];
            if (i == 0)
            {
                b_mass += add_aa_mass;
            }
        }
        if ((i >= link_site_1) && (i < link_site_2))
        {
            add_aa_mass += AA_mass_list[pep_sq[i - 1] - 65];
            if (i == link_site_1)
            {
                add_aa_mass += loop_mass;
            }
            b_mass += add_aa_mass;
        }
        else if (i > 0)
        {
            add_aa_mass += AA_mass_list[pep_sq[i - 1] - 65];
            if (fabs(add_aa_mass) > 1e-4)
            {
                b_mass += add_aa_mass;
                put_peptide_to_ion_index(ion_index_start, b_mass, multi, max_mass, pep_index_num);
            }
        }

    }

    for (i = pep_len + 1; i > 1; i--)
    {
        add_aa_mass = 0;
        if (mod_site_list[i] != -1)
        {
            add_aa_mass += mod_mass_list[mod_site_list[i]];
            if (i == pep_len + 1)
            {
                y_mass += add_aa_mass;
            }
        }
        if ((i > link_site_1) && (i <= link_site_2))
        {
            add_aa_mass += AA_mass_list[pep_sq[i - 1] - 65];
            if (i == link_site_2)
            {
                add_aa_mass += loop_mass;
            }
            y_mass += add_aa_mass;
        }
        else if (i < pep_len + 1)
        {
            add_aa_mass += AA_mass_list[pep_sq[i - 1] - 65];
            if (fabs(add_aa_mass) > 1e-4) 
            {
                y_mass += add_aa_mass;
                put_peptide_to_ion_index(ion_index_start, y_mass, multi, max_mass, pep_index_num);
            }
        }

    }
    return 0;
}

int tool_C_put_ion_to_ion_index_crosscross(ION_INDEX_MZ* ion_index_start, char* pep_sq, int pep_len, int* mod_site_list, double b_start_mass, double y_start_mass, double* mod_mass_list, double* AA_mass_list, int link_site_1, int link_site_2, int multi, int max_mass, long pep_index_num)
{
    // 要求link_site_1 小于link_site_1
    int i, j, k;
    double b_mass = b_start_mass;
    double y_mass = y_start_mass;

    double b_mz;
    double y_mz;

    double add_aa_mass = 0;

    for (i = 0; i < link_site_1; i++)
    {
        add_aa_mass = 0;
        if (mod_site_list[i] != -1)
        {
            add_aa_mass += mod_mass_list[mod_site_list[i]];
            if (i == 0)
            {
                b_mass += add_aa_mass;
            }
        }
        if (i > 0)
        {
            add_aa_mass += AA_mass_list[pep_sq[i - 1] - 65];
            if (fabs(add_aa_mass) > 1e-4) 
            {
                b_mass += add_aa_mass;
                put_peptide_to_ion_index(ion_index_start, b_mass, multi, max_mass, pep_index_num);
            }

        }
    }

    for (i = pep_len + 1; i >= link_site_2; i--)
    {
        add_aa_mass = 0;
        if (mod_site_list[i] != -1)
        {
            add_aa_mass += mod_mass_list[mod_site_list[i]];
            if (i == pep_len + 1)
            {
                y_mass += add_aa_mass;
            }
        }
        if (i < pep_len + 1) 
        {
            add_aa_mass += AA_mass_list[pep_sq[i - 1] - 65];
            if (fabs(add_aa_mass) > 1e-4)
            {
                y_mass += add_aa_mass;
                put_peptide_to_ion_index(ion_index_start, y_mass, multi, max_mass, pep_index_num);
            }
        }
    }
    return 0;
}

int tool_C_put_ion_to_ion_index_crossloop(ION_INDEX_MZ* ion_index_start, char* pep_sq, int pep_len, int* mod_site_list, double b_start_mass, double y_start_mass, double* mod_mass_list, double* AA_mass_list, int link_site_1, int link_site_2, double loop_mass, int cross_site, int multi, int max_mass, long pep_index_num)
{
    // 要求link_site_1 小于link_site_1
    int i, j, k;
    double b_mass = b_start_mass;
    double y_mass = y_start_mass;

    double b_mz;
    double y_mz;

    int left = 0, right = 0;
    double add_aa_mass = 0;
    for (i = 0; i <= pep_len + 1; i++)
    {
        add_aa_mass = 0;
        if (mod_site_list[i] != -1)
        {
            add_aa_mass += mod_mass_list[mod_site_list[i]];
            if (i == 0)
            {
                b_mass += add_aa_mass;
            }
        }
        if (i > 0)
        {
            add_aa_mass += AA_mass_list[pep_sq[i - 1] - 65];
            if (i == cross_site) 
            {
                break;
            }
            if ((i >= link_site_1) && (i < link_site_2))
            {
                if (i == link_site_1)
                {
                    add_aa_mass += loop_mass;
                }
                b_mass += add_aa_mass;
            }
            else 
            {
                if (fabs(add_aa_mass) > 1e-4) 
                {
                    b_mass += add_aa_mass;
                    put_peptide_to_ion_index(ion_index_start, b_mass, multi, max_mass, pep_index_num);
                }
            }
        }
    }

    for (i = pep_len + 1; i > 0; i--)
    {
        add_aa_mass = 0;
        if (mod_site_list[i] != -1)
        {
            add_aa_mass += mod_mass_list[mod_site_list[i]];
            if (i == pep_len + 1) 
            {
                y_mass += add_aa_mass;
            }
        }
        if (i < pep_len + 1)
        {
            if (i == cross_site)
            {
                break;
            }
            add_aa_mass += AA_mass_list[pep_sq[i - 1] - 65];
            if ((i >= link_site_1) && (i < link_site_2))
            {
                if (i == link_site_2)
                {
                    add_aa_mass += loop_mass;
                }
                y_mass += add_aa_mass;
            }
            else
            {
                if (fabs(add_aa_mass) > 1e-4)
                {
                    y_mass += add_aa_mass;
                    put_peptide_to_ion_index(ion_index_start, y_mass, multi, max_mass, pep_index_num);
                }
            }
        }
    }
    return 0;
}