// pch.h: 这是预编译标头文件。
// 下方列出的文件仅编译一次，提高了将来生成的生成性能。
// 这还将影响 IntelliSense 性能，包括代码完成和许多代码浏览功能。
// 但是，如果此处列出的文件中的任何一个在生成之间有更新，它们全部都将被重新编译。
// 请勿在此处添加要频繁更新的文件，这将使得性能优势无效。

#ifndef PCH_H
#define PCH_H

// 添加要在此处预编译的标头
#include "framework.h"
#include<stdio.h>
#include<malloc.h>
#include<math.h>
#include"data.h"
#include<stdlib.h>

extern "C" _declspec(dllexport) ION_INDEX_MZ * malloc_ion_index(int memory_size);
// 申请离子索引内存空间
extern "C" _declspec(dllexport) int free_ion_index(int memory_size, ION_INDEX_MZ * ion_index_start);
// 释放离子索引空间
extern "C" _declspec(dllexport) double* malloc_record_score_list(long max_pep_score_len, int all_num);
// 申请计分列表空间
extern "C" _declspec(dllexport) int free_record_score_list(double* record_score_list);
// 释放计分列表空间
extern "C" _declspec(dllexport) int get_ion_index_match_num(int* input_mz_list, int input_mz_list_len, ION_INDEX_MZ * mz_index_start, int max_mz, long max_pep_index_num, int is_ppm, double tol_fragement, double* record_match_score_list, int record_match_score_list_index);
// 计算匹配离子数目
extern "C" _declspec(dllexport) int get_ion_index_match_num_percent(long* input_pep_ion_num_list, int charge_num, double* record_match_score_list, int record_match_score_list_num_index, int record_match_score_list_percent_index, long max_pep_index_num);
// 计算匹配离子比列
extern "C" _declspec(dllexport) int get_ion_index_match_score(int* input_mz_list, int input_mz_list_len, double* input_mz_intensity_list, ION_INDEX_MZ * mz_index_start, int max_mz, long max_pep_index_num, int is_ppm, double tol_fragement, double* record_match_score_list, int record_match_score_list_percent_index, int record_match_score_list_index);
// 计算匹配离子分数
extern "C" _declspec(dllexport) int get_sum_score(double* record_match_score_list, long max_pep_index_num, int record_match_no_linker_index, int record_match_linker_index, int record_match_sum_index, double max_pep_score, long max_pep_score_index);
// 计算计分列表的分数和
extern "C" _declspec(dllexport) int get_peptide_mass(double* p_pep_mass, double* p_search_left_pep_mass, double* p_search_rigth_pep_mass, long* p_candidate_pep_index, int candidate_peppep_num, int is_ppm, double tol_precursor, double precursor_mass, double link_mass);
// 计算每条肽段的互补质量
extern "C" _declspec(dllexport) int find_max_value_index(double* record_score_list, long max_candidate_pep_num, int search_index);
// 查找最大值的索引号
extern "C" _declspec(dllexport) int read_ion_index_bin(char* path_write, ION_INDEX_MZ * data, int input_mz_list_len);
// 读入离子索引
extern "C" _declspec(dllexport) int filter_record_score_list(double* record_match_score_list, long max_pep_index_num, int filter_num_index, int filter_score_index, double filter_num_value, double filter_score_value, double* return_score, long* return_index, long* return_num);

extern "C" _declspec(dllexport) long* malloc_long_p(long max_candidate_pep_num, int all_num);

extern "C" _declspec(dllexport) int free_long_p(long* p_long);

extern "C" _declspec(dllexport) int quicksort(double* p_list, int start, int end);

extern "C" _declspec(dllexport) int bubblesort(double* p_list, int start, int end);

#endif //PCH_H
