// pch.cpp: 与预编译标头对应的源文件

#include "pch.h"
#include<time.h>
# include <stdlib.h>
// 当使用预编译的头时，需要使用此源文件，编译才能成功。

ION_INDEX_MZ* malloc_ion_index(int memory_size)
{
    // 反正就申请所有肽段可产生的所有碎片离子个数，超过就会被舍弃
    int i;
    ION_INDEX_MZ* ion_index_start;
    ion_index_start = (ION_INDEX_MZ*)malloc(sizeof(ION_INDEX_MZ) * memory_size);
    for (i = 0; i < memory_size; i++)
    {
        ion_index_start[i].pep_num = 0;
        ion_index_start[i].pep = NULL;
        ion_index_start[i].pep_size = 0;
    }
    printf("[Info] Finish Malloc ion_dic_mz, size = %ld\n", memory_size);
    return ion_index_start;
};

int free_ion_index(int memory_size, ION_INDEX_MZ* ion_index_start)
{
    clock_t start, end;
    int i, j, k;
    long* pep_start;
    start = clock();
    for (i = 0; i < memory_size; i++)
    {
        pep_start = ion_index_start[i].pep;
        free(pep_start);
    }
    free(ion_index_start);
    end = clock();
    printf("[Info] Finish free ion_dic_mz !\n");
    printf("[Info] Free fragment index time : %.2lf s\n", (double)(end - start) / CLOCKS_PER_SEC);
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

int write_ion_index_bin(const char* path_write, ION_INDEX_MZ* data, int input_mz_list_len)
{
    clock_t start, end;
    clock_t tmp_1, tmp_2;

    FILE* fp;
    int i, j, k;
    long one_mz_pep_num; // 一个质荷比候选的肽段数目
    long one_mz_pep_index_num; // 一个肽段编号
    long* p_write_long;

    fopen_s(&fp, path_write, "wb");
    //先写入每个质荷比数目的这个结构，方便读取的时候申请空间存
    start = clock();
    for (i = 0; i < input_mz_list_len; i++)
    {
        fwrite(&(data[i].pep_num), sizeof(long), 1, fp);
    }
    end = clock();
    start = clock();
    for (i = 0; i < input_mz_list_len; i++)
    {
        tmp_1 = clock();
        if (data[i].pep_num > 0)
        {
            fwrite(data[i].pep, sizeof(long) * data[i].pep_num, 1, fp);
        }
    }
    end = clock();
    printf("[Info] Writting fragment index time : %.2lf s\n", (double)(end - start) / CLOCKS_PER_SEC);
    fclose(fp);
    return 0;
}

