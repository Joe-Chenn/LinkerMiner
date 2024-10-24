// pch.h: 这是预编译标头文件。
// 下方列出的文件仅编译一次，提高了将来生成的生成性能。
// 这还将影响 IntelliSense 性能，包括代码完成和许多代码浏览功能。
// 但是，如果此处列出的文件中的任何一个在生成之间有更新，它们全部都将被重新编译。
// 请勿在此处添加要频繁更新的文件，这将使得性能优势无效。

#ifndef PCH_H
#define PCH_H

// 添加要在此处预编译的标头
#include "framework.h"
#include "data.h"
#endif //PCH_H

int put_peptide_to_ion_index(ION_INDEX_MZ* ion_index_start, double input_mz, int multi, int max_mz, long pep_index_num);

extern "C" _declspec(dllexport) double calculate_peptide_mass(char* p_sq, int sq_len, double input_mass, double* AA_mass);

extern "C" _declspec(dllexport) double calculate_gdm(char* p_sq, int sq_len, int* p_mod_site_list, int fix_mod_num, double* p_gdm_matrix, int AA_num, int gdm_max_sq_len);

extern "C" _declspec(dllexport) long long tool_C_get_ion_code(long long pkl_file_index, long long ind_file_index, long long one_ind_num, long long link_site_combination_num);

extern "C" _declspec(dllexport) long* tool_C_new_get_ion_data(long long input_data);

extern "C" _declspec(dllexport) int tool_C_put_ion_to_ion_index(ION_INDEX_MZ * ion_index_start, char* pep_sq, int pep_len, int* mod_site_list, double b_start_mass, double y_start_mass, double* mod_mass_list, double* AA_mass_list, int link_site, int multi, int max_mass, long pep_index_num);

extern "C" _declspec(dllexport) int tool_C_put_ion_to_ion_index_crosslink(ION_INDEX_MZ * ion_index_start, char* pep_sq, int pep_len, int* mod_site_list, double b_start_mass, double y_start_mass, double* mod_mass_list, double* AA_mass_list, int link_site, int multi, int max_mass, long pep_index_num);

extern "C" _declspec(dllexport) int tool_C_put_ion_to_ion_index_singlepeptide(ION_INDEX_MZ * ion_index_start, char* pep_sq, int pep_len, int* mod_site_list, double b_start_mass, double y_start_mass, double* mod_mass_list, double* AA_mass_list, int multi, int max_mass, long pep_index_num);

extern "C" _declspec(dllexport) int tool_C_put_ion_to_ion_index_looplink(ION_INDEX_MZ * ion_index_start, char* pep_sq, int pep_len, int* mod_site_list, double b_start_mass, double y_start_mass, double* mod_mass_list, double* AA_mass_list, int link_site_1, int link_site_2, double loop_mass, int multi, int max_mass, long pep_index_num);

extern "C" _declspec(dllexport) int tool_C_put_ion_to_ion_index_crosscross(ION_INDEX_MZ * ion_index_start, char* pep_sq, int pep_len, int* mod_site_list, double b_start_mass, double y_start_mass, double* mod_mass_list, double* AA_mass_list, int link_site_1, int link_site_2, int multi, int max_mass, long pep_index_num);

extern "C" _declspec(dllexport) int tool_C_put_ion_to_ion_index_crossloop(ION_INDEX_MZ * ion_index_start, char* pep_sq, int pep_len, int* mod_site_list, double b_start_mass, double y_start_mass, double* mod_mass_list, double* AA_mass_list, int link_site_1, int link_site_2, double loop_mass, int cross_site, int multi, int max_mass, long pep_index_num);