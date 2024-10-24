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
#include<data.h>
#include<malloc.h>

extern "C" _declspec(dllexport) ION_INDEX_MZ * malloc_ion_index(int memory_size);

extern "C" _declspec(dllexport) int free_ion_index(int memory_size, ION_INDEX_MZ * ion_index_start);

extern "C" _declspec(dllexport) int put_peptide_to_ion_index(ION_INDEX_MZ * ion_index_start, double input_mz, int multi, int max_mz, long pep_index_num);

extern "C" _declspec(dllexport) int write_ion_index_bin(const char* path_write, ION_INDEX_MZ * data, int input_mz_list_len);

#endif //PCH_H
