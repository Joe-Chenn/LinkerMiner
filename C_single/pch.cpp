// pch.cpp: 与预编译标头对应的源文件

#include "pch.h"
# include<time.h>

// 当使用预编译的头时，需要使用此源文件，编译才能成功。

// 第一步申请离子索引的内存，反正就把肽段的号填进去，用记录离子索引分布的索引来划分属于那些质荷比上的东西

ION_INDEX_MZ* malloc_ion_index(int memory_size)
{
    // 反正就申请所有肽段可产生的所有碎片离子个数，超过就会被舍弃
    long i;
    ION_INDEX_MZ* ion_index_start;
    ion_index_start = (ION_INDEX_MZ*)malloc(sizeof(ION_INDEX_MZ) * memory_size);
    printf("[Info] Finish Malloc ion_dic_mz, size = %ld\n", memory_size);
    return ion_index_start;
};

int free_ion_index(int memory_size, ION_INDEX_MZ* ion_index_start)
{
    int i, j, k;
    long* pep_node;
    for (i = 0; i < memory_size; i++)
    {
        pep_node = ion_index_start[i].pep;
        free(pep_node);
    }
    free(ion_index_start);
    printf("[Info] Finish free ion_dic_mz !\n");
    return 0;
}

double* malloc_record_score_list(long max_pep_score_len, int all_num)
{
    int i, j, k;
    double* p_record;
    if (max_pep_score_len == 0)
    {
        max_pep_score_len = 1;
    }
    p_record = (double*)malloc(sizeof(double) * all_num * max_pep_score_len);
    // 看着申请用，咋填都行，是候选肽段数目的整数倍内存空间
    for (i = 0; i < max_pep_score_len; i++)
    {
        for (j = 0; j < all_num; j++)
        {
            p_record[i + max_pep_score_len * j] = 0.0;
        }
    }
    return p_record;
}

int free_record_score_list(double* record_score_list)
{
    free(record_score_list);
    return 0;
}

int get_ion_index_match_num(int* input_mz_list, int input_mz_list_len, ION_INDEX_MZ* mz_index_start, int max_mz, long min_pep_index_num, long max_pep_index_num, int is_ppm, double tol_fragement, double* record_match_score_list, int record_match_ion_num_index)
{
    int i, j, k;
    int one_mz, one_mz_start, one_mz_end;
    int one_mz_index_num;
    int match_num_add;

    match_num_add = max_pep_index_num * record_match_ion_num_index;
    for (i = 0; i < input_mz_list_len; i++)
    {
        one_mz = input_mz_list[i];
        one_mz_start = input_mz_list[i];
        one_mz_end = input_mz_list[i];
        // 确定偏差限制下可查找的范围
        if (is_ppm == 1)
        {
            one_mz_start = floor(one_mz_start - (one_mz_start * tol_fragement));
            one_mz_end = ceil(one_mz_end + (one_mz_end * tol_fragement));
        }
        else
        {
            one_mz_start = floor(one_mz_start - tol_fragement);
            one_mz_end = ceil(one_mz_end + tol_fragement);
        }
        if (one_mz >= max_mz)
        {
            break;
        }
        for (j = one_mz_start; j < one_mz_end + 1; j++)
        {
            for (k = 0; k < mz_index_start[j].pep_num; k++)
            {
                one_mz_index_num = mz_index_start[j].pep[k];
                if (one_mz_index_num >= max_pep_index_num)
                {
                    break;
                }
                if (one_mz_index_num < min_pep_index_num)
                {
                    continue;
                }
                record_match_score_list[one_mz_index_num + match_num_add - min_pep_index_num] += 1;
            }
        }
    }
    return 0;
}

int get_ion_index_match_num_percent(long* input_pep_ion_num_list, double* record_match_score_list, int record_match_score_list_num_index, int record_match_score_list_percent_index, long min_pep_index_num, long max_pep_index_num)
{
    int i, j, k;
    int match_num_add, percent_add;
    match_num_add = (max_pep_index_num - min_pep_index_num) * record_match_score_list_num_index;
    percent_add = (max_pep_index_num - min_pep_index_num) * record_match_score_list_percent_index;
    for (i = 0; i < max_pep_index_num - min_pep_index_num; i++)
    {
        if (record_match_score_list[i + match_num_add] < 0.0001)
        {
            record_match_score_list[i + percent_add] = 10000;
        }
        else
        {
            record_match_score_list[i + percent_add] = record_match_score_list[i + match_num_add] / input_pep_ion_num_list[i + min_pep_index_num];
        }
    }
    return 0;
}

int get_ion_index_match_score(int* input_mz_list, int input_mz_list_len, double* input_mz_intensity_list, ION_INDEX_MZ* mz_index_start, int max_mz, long min_pep_index_num, long max_pep_index_num, int is_ppm, double tol_fragement, double* record_match_score_list, int record_match_score_list_percent_index, int record_match_score_list_index)
{
    int i, j, k;
    double k1 = 0.001, b = -25;
    int one_mz, one_mz_start, one_mz_end;
    int one_mz_index_num;
    double one_mz_bias_score, one_mz_inten_score, one_mz_score, one_mz_bias, one_mz_start_bias, one_mz_end_bias, max_mz_bias;
    int percent_add, score_add;
    percent_add = (max_pep_index_num - min_pep_index_num) * record_match_score_list_percent_index;
    score_add = (max_pep_index_num - min_pep_index_num) * record_match_score_list_index;
    for (i = 0; i < input_mz_list_len; i++)
    {
        //强度分数
        one_mz_inten_score = sin(input_mz_intensity_list[i] * 1.57075);
        //计算偏差
        one_mz = input_mz_list[i];
        one_mz_start = input_mz_list[i];
        one_mz_end = input_mz_list[i];
        // 确定偏差限制下可查找的范围
        if (is_ppm == 1)
        {
            one_mz_start = floor(one_mz_start - (one_mz_start * tol_fragement));
            one_mz_end = ceil(one_mz_end + (one_mz_end * tol_fragement));
        }
        else
        {
            one_mz_start = floor(one_mz_start - tol_fragement);
            one_mz_end = ceil(one_mz_end + tol_fragement);
        }
        one_mz_start_bias = abs(one_mz_start - one_mz);
        one_mz_end_bias = abs(one_mz_end - one_mz);

        if (one_mz_start_bias > one_mz_end_bias)
        {
            max_mz_bias = one_mz_start_bias;
        }
        else
        {
            max_mz_bias = one_mz_end_bias;
        }
        for (j = one_mz_start; j < one_mz_end + 1; j++)
        {
            one_mz_bias = double(abs(j - one_mz));
            one_mz_bias_score = log(2.718281828459 - one_mz_bias / max_mz_bias);
            for (k = 0; k < mz_index_start[j].pep_num; k++)
            {
                //one_mz_index_num = mz_index_start[j].pep[k].pep_index_num;
                one_mz_index_num = mz_index_start[j].pep[k];
                if (one_mz_index_num >= max_pep_index_num)
                {
                    break;
                }
                if (one_mz_index_num < min_pep_index_num)
                {
                    continue;
                }
                one_mz_score = one_mz_bias_score * one_mz_inten_score * (1 + k1) / (one_mz_inten_score + k1 * (record_match_score_list[one_mz_index_num + percent_add - min_pep_index_num] * b + 1 - b));
                record_match_score_list[one_mz_index_num + score_add - min_pep_index_num] = record_match_score_list[one_mz_index_num + score_add - min_pep_index_num] + one_mz_score;
            }
        }
    }
    return 0;
}

int get_sum_score(double* record_match_score_list, long max_pep_index_num, int record_match_no_linker_index, int record_match_linker_index, int record_match_sum_index)
{
    int i;
    int sum_add, no_linker_add, linker_add;
    sum_add = max_pep_index_num * record_match_sum_index;
    no_linker_add = max_pep_index_num * record_match_no_linker_index;
    linker_add = max_pep_index_num * record_match_linker_index;

    for (i = 0; i < max_pep_index_num; i++)
    {
        record_match_score_list[i + sum_add] = record_match_score_list[i + no_linker_add] + record_match_score_list[i + linker_add];
    }
    return 0;
}

int read_ion_index_bin(char* path_write, ION_INDEX_MZ* data, int input_mz_list_len)
{
    clock_t start, end;
    int i, j, k;
    FILE* fp;
    long* one_mz_pep_index_num;
    long* mz_matrix;
    mz_matrix = (long*)malloc(sizeof(long) * input_mz_list_len);
    fopen_s(&fp, path_write, "rb");
    fread(mz_matrix, sizeof(long) * input_mz_list_len, 1, fp);
    start = clock();
    for (i = 0; i < input_mz_list_len; i++)
    {
        data[i].pep_num = mz_matrix[i];
        if (data[i].pep_num > 0)
        {
            one_mz_pep_index_num = (long*)malloc(sizeof(long) * data[i].pep_num);
            fread(one_mz_pep_index_num, sizeof(long), data[i].pep_num, fp);
            data[i].pep = one_mz_pep_index_num;
        }
    }
    end = clock();
    printf("[Info] Reading fragment index time : %.2lf s\n", (double)(end - start) / CLOCKS_PER_SEC);
    fclose(fp);
    return 0;
}
