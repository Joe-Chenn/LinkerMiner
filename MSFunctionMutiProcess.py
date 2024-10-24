# -*- mode: python ; coding: utf-8 -*-
import ctypes
import math
from multiprocessing.pool import Pool
from multiprocessing.pool import ThreadPool

from CTypesStructure import ION_INDEX_MZ
from MSTool import tool_binary_search, tool_binary_search_index, flatten_and_dump
from MSLogging import logGetError
from MSSysterm import DLL_CROSS, FOLDER_INDEX, FILE_NAME_ION_INDEX
from MSDataC import ION_INDEX_MZ_LOCK, ION_INDEX_PEP_LOCK
from ctypes import *

import os
import numpy as np
import pickle
import operator
import time


class CFunctionMultiProcessSortPKL:

    def __init__(self):
        pass

    def MutiProcessSortPKL(self, pkl_file_list, number_thread, type_thread, output_ind_file=False):

        if type_thread == 0:
            pool = Pool(processes=number_thread)  # multiply processes
        else:
            pool = ThreadPool(number_thread)  # multiply threads
        # data = pool.map(self.process_sort_one_pkl, pkl_file_list)

        '''20220826新版本测试'''
        if output_ind_file:
            data = pool.map(self.process_sort_var_pep_index, pkl_file_list)
        else:
            data = pool.map(self.process_sort_fix_pep_index, pkl_file_list)
        '''20220826新版本测试'''

        pool.close()
        pool.join()

        pep_num = 0
        for one_data in data:
            pep_num += one_data
        return pep_num

    def process_sort_fix_pep_index(self, one_pkl_file_data):
        # (pkl路径，config）
        pkl_file, ind_file = one_pkl_file_data
        # 这里是提取pkl文件名而已
        # 不要后缀
        all_peptides = []
        # 一行位置存一个类
        # 所属蛋白索引，肽段在蛋白中的起始位置，肽段在蛋白中的结束位置，肽段质量，修饰list修饰(一行是一个位点pkl类：修饰信息位置索引，位点# mod_index, site），哥德尔编码值
        f = open(pkl_file, 'rb')
        while True:
            try:
                one_peptide = pickle.load(f)
                all_peptides.append(one_peptide)
            except Exception:  # the end of the file
                break
        f.close()
        print(pkl_file + str(len(all_peptides)))
        # pkl中的文件是很多个类
        # 所属蛋白索引，肽段在蛋白中的起始位置，肽段在蛋白中的结束位置，肽段质量，修饰list修饰(一行是一个位点pkl类：修饰信息位置索引，位点# mod_index, site），哥德尔编码值
        cmpfun = operator.attrgetter('mass', 'gdm')
        all_peptides.sort(key=cmpfun)  # 按mass和gdm排序这个pkl文件内的类内容
        old_num = len(all_peptides)  # 酶切后的数量
        new_all_peptides = []
        self.__captain_update_peptide(all_peptides, new_all_peptides)  # 整理去重的肽段
        del all_peptides
        print("[Info]#Peptides in {0} file: {1} to {2}".format(pkl_file, old_num, len(new_all_peptides)))
        # 所属蛋白索引，肽段在蛋白中的起始位置，肽段在蛋白中的结束位置，肽段质量，修饰list修饰(一行是一个位点pkl类：修饰信息位置索引，位点# mod_index, site），哥德尔编码值
        # 相同肽段的归类，在所属蛋白索引里会有多个值，代表不同的蛋白
        fw = open(pkl_file, 'wb')
        fw.close()
        fw = open(pkl_file, 'ab')
        for write_all_peptides in new_all_peptides:
            pickle.dump(write_all_peptides, fw, protocol=pickle.HIGHEST_PROTOCOL)
        fw.close()
        return len(new_all_peptides)

    def process_sort_var_pep_index(self, one_pkl_file_data):
        # (pkl路径，config）
        pkl_file, ind_file = one_pkl_file_data
        start_mass, end_mass = 0.0, 0.0
        if pkl_file.find('/') >= 0:
            start_end_str = pkl_file[pkl_file.rfind('/') + 1:]
        else:
            start_end_str = pkl_file[pkl_file.rfind('\\') + 1:]
        # 返回字符串最后一次出现的位置(从右向左查询)，如果没有匹配项则返回-1
        # 这里是提取pkl文件名而已
        start_end_str = start_end_str[:start_end_str.find('.')]
        # 不要后缀
        index_str = start_end_str.find('_')
        # _这个是质量分割符
        # 表示这个pkl文件对应的质量区间
        start_mass = int(start_end_str[:index_str]) * 1.0
        end_mass = int(start_end_str[index_str + 1:]) * 1.0

        all_peptides = []
        # 一行位置存一个类
        # 所属蛋白索引，肽段在蛋白中的起始位置，肽段在蛋白中的结束位置，肽段质量，修饰list修饰(一行是一个位点pkl类：修饰信息位置索引，位点# mod_index, site），哥德尔编码值
        f = open(pkl_file, 'rb')
        while True:
            try:
                one_peptide = pickle.load(f)
                if type(one_peptide) is list:
                    all_peptides += one_peptide
                else:
                    all_peptides.append(one_peptide)
            except Exception:  # the end of the file
                break
        f.close()

        # pkl中的文件是很多个类
        # 所属蛋白索引，肽段在蛋白中的起始位置，肽段在蛋白中的结束位置，肽段质量，修饰list修饰(一行是一个位点pkl类：修饰信息位置索引，位点# mod_index, site），哥德尔编码值
        cmpfun = operator.attrgetter('mass', 'gdm')
        all_peptides.sort(key=cmpfun)  # 按mass和gdm排序这个pkl文件内的类内容
        old_num = len(all_peptides)  # 酶切后的数量
        new_all_peptides = []
        self.__captain_update_peptide(all_peptides, new_all_peptides)  # 整理去重的肽段
        del all_peptides
        print("[Info]#Peptides in {0} file: {1} to {2}".format(pkl_file, old_num, len(new_all_peptides)))
        # 所属蛋白索引，肽段在蛋白中的起始位置，肽段在蛋白中的结束位置，肽段质量，修饰list修饰(一行是一个位点pkl类：修饰信息位置索引，位点# mod_index, site），哥德尔编码值
        # 相同肽段的归类，在所属蛋白索引里会有多个值，代表不同的蛋白
        # print("[Info]#Peptides in {0} file: {1} to {2}".format(pkl_file, old_num, len(all_peptides)))
        index_list = self.__captain_create_peptide_or_spectrum_index(new_all_peptides, start_mass, end_mass,
                                                                     1)  # config.D12_MULTI_MASS)
        # 去重后的的list，一行是一个肽段的类：所属蛋白索引，肽段在蛋白中的起始位置，肽段在蛋白中的结束位置，肽段质量，修饰list修饰(一行是一个位点pkl类：修饰信息位置索引，位点# mod_index, site），哥德尔编码值

        # 起始质量，结束质量
        # index_list每一行表示到这个小区间已经有多少个肽段，整个大质量区间又被分成很多个小区间
        fw = open(pkl_file, 'wb')
        pickle.dump(new_all_peptides, fw)  # 把去重后的按质量排好序的肽段写入原始pkl文件，这是质量对肽段
        '''20220825更新'''
        # for write_all_peptides in new_all_peptides:
        #     pickle.dump(write_all_peptides, fw, protocol=pickle.HIGHEST_PROTOCOL)
        '''20220825更新'''
        fw.close()

        fw = open(ind_file, 'wb')
        pickle.dump(index_list, fw)  # 新pkl文件加了_ind.pkl这个存的是对应的大区间中的小区间肽段分分布索引，哪个质量能返回哪些肽段，质量到肽段的索引
        fw.close()

        return len(new_all_peptides)

    def __captain_update_peptide(self, all_peptides, new_all_peptides):
        # list 一行位置存一个类，表示一种修饰的肽段，这个已经按质量和哥德尔编码排好序了
        # 所属蛋白索引，肽段在蛋白中的起始位置，肽段在蛋白中的结束位置，肽段质量，修饰list修饰(一行是一个位点pkl类：修饰信息位置索引，位点# mod_index, site），哥德尔编码值
        # pro_index, start_pos, end_pos, mass, mods, gdm
        # mods : mod_index, site
        # all_peptides is the list of peptides sorted by masses and gdm values
        # return the list of peptides after removing the peptides with the same gdm values and update the protein indexes
        start_i = 0  # 索引all_peptides这个list，相当于索引肽段
        while start_i < len(all_peptides):

            one_pep2pro_is_C_term = 0
            one_pep2pro_is_N_term = 0
            is_target = 0

            new_all_peptides.append(all_peptides[start_i])
            cur_p = all_peptides[start_i]  # 一个肽段的信息
            cur_mass = cur_p.mass  # 一个肽段质量
            cur_gdm = cur_p.gdm  # 一个肽段哥德尔编码值
            start_j = start_i + 1  # 下一个肽段索引
            cur_pro_index_list = [cur_p.pro_index]  # 上一个肽段对应的蛋白索引
            cur_pro_index_start_list = [cur_p.start_pos]
            # 因为排好序，下一个和这个不同那么后面都不会有与这个相同的了

            # 只要候选结果里存在就是可以认为是C端或者N端能交联上
            if all_peptides[start_i].start_pos == 0:
                one_pep2pro_is_N_term = 1
            if all_peptides[start_i].end_pos == all_peptides[start_i].len_pro_sq:
                one_pep2pro_is_C_term = 1
            if all_peptides[start_i].is_target == 1:
                is_target = 1
            while start_j < len(all_peptides) and all_peptides[start_j].mass == cur_mass and all_peptides[
                start_j].gdm == cur_gdm:
                # 下一个肽段的质量以及哥德尔编码值与上一个肽段值一样
                if all_peptides[start_j].start_pos == 0:
                    one_pep2pro_is_N_term = 1
                if all_peptides[start_j].end_pos == all_peptides[start_j].len_pro_sq:
                    one_pep2pro_is_C_term = 1
                if all_peptides[start_j].is_target == 1:
                    is_target = 1
                cur_pro_index_list.append(all_peptides[start_j].pro_index)  # 记录一样的肽段以及修饰一样的肽段来组不同的蛋白的索引号
                cur_pro_index_start_list.append(all_peptides[start_j].start_pos)
                start_j += 1  # 再下一个肽段
            # 只要候选结果里存在就是可以认为是C端或者N端能交联上
            new_all_peptides[-1].pro_index_list = cur_pro_index_list  # 记录酶解修饰相同的肽段不同的蛋白
            new_all_peptides[-1].pro_index_start_pos = cur_pro_index_start_list  # 记录酶解修饰相同的肽段不同的蛋白起始位点
            new_all_peptides[-1].pep2pro_is_C_term = one_pep2pro_is_C_term
            new_all_peptides[-1].pep2pro_is_N_term = one_pep2pro_is_N_term
            new_all_peptides[-1].is_target = is_target
            start_i = start_j
        # 整理好去重后的list

    def __captain_create_peptide_or_spectrum_index(self, all_peptides, start_mass, end_mass, mul=1):
        # 去重后的的list，一行是一个肽段的类：所属蛋白索引，肽段在蛋白中的起始位置，肽段在蛋白中的结束位置，肽段质量，修饰list修饰(一行是一个位点pkl类：修饰信息位置索引，位点# mod_index, site），哥德尔编码值
        # 起始质量，结束质量，
        # 也按质量排好序了
        # all_peptides is the list of peptides sorted by masses and gdm values
        # return the list of mass index to the all_peptides
        num = int((end_mass - start_mass) * mul) + 1  # the number of list 区间按成0.01Da质量差进行分的份数
        index_list = [-1 for i in range(num)]  # index_list[X] is the min index in all_peptides whose mass is >= X
        for i in range(len(all_peptides)):  # 遍历一个质量区间的所有肽段，i表示第几个肽段
            p = all_peptides[i]  # 这个肽段的质量相对于区间起始质量的差*mul
            index = int((p.mass - start_mass) * mul)  # 计算属于哪个区间的质量
            if index >= num: index = num - 1
            if index_list[index] == -1:
                index_list[index] = i
        # 这样得到的index_list表示这质量区间内，再分的小区间，其中对应的值是这个区间所在的质量在all_peptides中的右区间边界
        end_val = len(all_peptides)
        for i in range(num)[::-1]:
            if index_list[i] == -1:
                index_list[i] = end_val  # index_list每一行表示到这个小区间已经有多少个肽段
            else:
                end_val = index_list[i]
        return index_list


class CFunctionThreadCoarseSearch:

    def __init__(self, inputDP,
                 list_ion_num, list_pep_mass, list_pep_data_pkl, list_link_site,
                 linker_mass,
                 list_pepcode2resultindex):

        self.dp = inputDP

        # 这里是做函数的输入输出类型的定义
        self.list_ion_num = list_ion_num
        self.list_pep_mass = list_pep_mass
        self.list_pep_data_pkl = list_pep_data_pkl
        self.list_link_site = list_link_site

        self.linker_mass = linker_mass

        self.list_pepcode2resultindex = list_pepcode2resultindex  # 记录记结果所在位置

    def MutiThreadOnlyCrossCoarseSearch(self, list_coarse_data, coarse_res, path_ion_index, multithread=True):

        multithread = self.dp.myCFG.D2_TYPE_THREAD

        time_1 = time.time()
        # 读离子索引
        bird = self.import_DLL()
        self.__captain_malloc_C_data(bird, path_ion_index)
        print("[Info] Start Coarse searching ... ")
        if multithread:
            # 分发谱图
            split_coarse_MS2_data = []
            self.__captian_split_coarse_MS2_data(list_coarse_data, split_coarse_MS2_data,
                                                 self.dp.myCFG.D1_NUMBER_THREAD)

            pool = ThreadPool(self.dp.myCFG.D1_NUMBER_THREAD)  # multiply threads
            pool.map(self.coarse_only_cross_search, split_coarse_MS2_data)
            pool.close()
            pool.join()
            for one_data in split_coarse_MS2_data:
                coarse_res += one_data[2]
        else:

            split_coarse_MS2_data = []
            self.__captian_split_coarse_MS2_data(list_coarse_data, split_coarse_MS2_data, 1)
            for part_ms2 in split_coarse_MS2_data:
                self.coarse_only_cross_search(part_ms2)
                coarse_res += part_ms2[2]
        # 释放离子索引
        self.__captain_free_C_data(bird)
        time_2 = time.time()
        print("[Info] Finish crosslink search ! using time %.4f s" % (time_2 - time_1))
    def dll_test(self):
        dll = CDLL(os.getcwd() + '\\' + "CrossDLL.dll")
        dll.print.argtypes = (c_int, c_int)
        dll.print.restype = (c_int)
        test = dll.print(1, 4)
        print(test)
    def import_DLL(self):

        # 导入动态链接库
        bird = CDLL(os.getcwd() + '\\' + DLL_CROSS)

        # 这里是做函数的输入输出类型的定义
        bird.malloc_ion_index.argtype = (c_int)
        bird.malloc_ion_index.restype = (POINTER(ION_INDEX_MZ_LOCK))

        bird.free_ion_index.argtypes = (c_int, POINTER(ION_INDEX_MZ_LOCK))
        bird.free_ion_index.restype = (c_int)

        bird.malloc_record_score_list.argtype = (c_long, c_int)
        bird.malloc_record_score_list.restype = (POINTER(c_double))

        bird.free_record_score_list.argtype = (POINTER(c_double))
        bird.free_record_score_list.restype = (c_int)

        bird.get_ion_index_match_num.argtype = (
            POINTER(c_int), c_int, POINTER(ION_INDEX_MZ_LOCK), c_int, c_long, c_int, c_double, POINTER(c_double), c_int)
        bird.get_ion_index_match_num.restype = (c_int)

        bird.get_ion_index_match_num_percent.argtype = (POINTER(c_long), POINTER(c_double), c_int, c_int, c_int)
        bird.get_ion_index_match_num_percent.restype = (c_int)

        bird.get_ion_index_match_score.argtype = (
            POINTER(c_int), c_int, POINTER(c_double), POINTER(ION_INDEX_MZ_LOCK), c_int, c_long, c_int, c_double,
            POINTER(c_double), c_int, c_int)
        bird.get_ion_index_match_score.restype = (c_int)

        bird.get_sum_score.argtype = (POINTER(c_double), c_long, c_int, c_int, c_int, c_double, c_long)
        bird.get_sum_score.restype = (c_int)

        bird.get_peptide_mass.argtype = (
            POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_long), c_int, c_int, c_double, c_double,
            c_double)
        bird.get_peptide_mass.restype = (c_int)

        bird.find_max_value_index.argtype = (POINTER(c_double), c_int, c_int)
        bird.find_max_value_index.restype = (c_int)

        bird.read_ion_index_bin.argtype = (POINTER(c_char_p), POINTER(ION_INDEX_MZ_LOCK), c_int)
        bird.read_ion_index_bin.restype = (c_int)

        bird.filter_record_score_list.argtype = (
            POINTER(c_double), c_long, c_int, c_int, c_double, c_double, POINTER(c_double), POINTER(c_long),
            POINTER(c_int))
        bird.filter_record_score_list.restype = (c_int)

        bird.malloc_long_p.argtype = (c_long, c_int)
        bird.malloc_long_p.restype = (POINTER(c_long))

        bird.free_long_p.argtype = (POINTER(c_long))
        bird.free_long_p.restype = (c_int)

        bird.quicksort.argtype = (POINTER(c_double), c_int, c_int)
        bird.quicksort.restype = (c_int)

        bird.bubblesort.argtype = (POINTER(c_double), c_int, c_int)
        bird.bubblesort.restype = (c_int)

        return bird

    def __captain_malloc_C_data(self, bird, path_ion_index):

        # 提前处理的数据
        if self.dp.myCFG.is_ppm_pre:
            self.is_ppm_pre = 1
            self.fragment_ppm = c_double(self.dp.myCFG.C2_PPM_TOL_FRAGMENT)
            self.precursor_ppm = c_double(self.dp.myCFG.C1_PPM_TOL_PRECURSOR)
        else:
            self.is_ppm_pre = 0
            self.fragment_ppm = c_double(self.dp.myCFG.C2_PPM_TOL_FRAGMENT)
            self.precursor_ppm = c_double(self.dp.myCFG.C1_PPM_TOL_PRECURSOR)
        # 提前处理的数据

        # 需要用到的指针
        c_double_list_ion_num = (c_long * len(self.list_ion_num))()
        for n in range(len(self.list_ion_num)):
            c_double_list_ion_num[n] = self.list_ion_num[n]
        p_C_list_ion_num = byref(c_double_list_ion_num)
        self.p_C_list_ion_num = p_C_list_ion_num

        p_C_ion_index_start = bird.malloc_ion_index(self.dp.myCFG.max_mz_index_size)
        bird.read_ion_index_bin(path_ion_index, p_C_ion_index_start, self.dp.myCFG.max_mz_index_size)

        self.p_C_ion_index_start = p_C_ion_index_start
        # self.__captain_get_numpy_data()

    # def __captain_dump_ion_num(self):

    def __captain_get_numpy_data(self):

        ivf_array = ctypes.cast(self.p_C_ion_index_start,
                                ctypes.POINTER(ION_INDEX_MZ * self.dp.myCFG.max_mz_index_size)).contents
        ivf_dic = []
        for item in ivf_array:
            print(item.pep_num)
            if item.pep_num != 0:
                nums_array = ctypes.cast(item.pep, ctypes.POINTER(ctypes.c_long * item.pep_num)).contents
                nums_list = [nums_array[i] for i in range(item.pep_num)]
            else:
                nums_list = []
            ivf_dic.append(nums_list)
        dump_path = self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[2][1] + '\\' + \
                    FILE_NAME_ION_INDEX[2][1] + '.np.ivf'
        flatten_and_dump(ivf_dic, dump_path)

    def __captain_free_C_data(self, bird):

        bird.free_ion_index(self.dp.myCFG.max_mz_index_size, self.p_C_ion_index_start)

    def __captian_split_coarse_MS2_data(self, list_coarse_data, split_coarse_MS2_data, part_num):

        # 记粗打分结果，coarse_alpha_res和coarse_beta_res是记录结果来自于所有结果里的哪个位置
        # 建这两个list只是为了记这个肽段有没有被用过
        split_num = int(len(list_coarse_data) // part_num)
        list_start_num = [n * split_num for n in range(1, part_num)]
        list_start_num.append(len(list_coarse_data))

        spectrum_num = 0
        part_index = 0
        part_spectrum_data = []
        cur_part_spectrum = 0

        while spectrum_num < len(list_coarse_data):

            if list_start_num[part_index] == spectrum_num:
                tmp_res = [[] for i in range(cur_part_spectrum)]
                split_coarse_MS2_data.append([part_spectrum_data, part_index, tmp_res])
                cur_part_spectrum = 0

            part_spectrum_data.append(list_coarse_data[spectrum_num])
            cur_part_spectrum += 1
            spectrum_num += 1

        tmp_res = [[] for i in range(cur_part_spectrum)]
        split_coarse_MS2_data.append([part_spectrum_data, part_index, tmp_res])

    def coarse_only_cross_search(self, input_data):

        spectrum_data = input_data[0]
        part_index = input_data[1]
        coarse_res = input_data[2]

        finish_spectrum_search_num = 0
        # self.dll_test()
        bird = self.import_DLL()
        time_1 = time.time()

        # 记录每个区间的肽段数目
        mass_index = [0 for i in range(int(self.dp.myCFG.D7_MASS_PEP_UP) + 1)]
        cur_sum_pep_num = 0
        one_pep_index = 0
        i = self.dp.myCFG.D6_MASS_PEP_LOW
        for i in range(int(self.dp.myCFG.D6_MASS_PEP_LOW), int(self.dp.myCFG.D7_MASS_PEP_UP), 1):
            while True:
                if one_pep_index == len(self.list_pep_mass):
                    break
                if self.list_pep_mass[one_pep_index] > i + 1:
                    break
                one_pep_index += 1
                cur_sum_pep_num += 1
            mass_index[i] = cur_sum_pep_num

        for one_pep_index, one_pep_mass in enumerate(self.list_pep_mass):
            if one_pep_mass > i + 1:
                mass_index[i] = cur_sum_pep_num
                i = i + 1
            cur_sum_pep_num += 1
        mass_index[i] = cur_sum_pep_num
        filepath = self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[2][1] + '\\' + "mass_index.pk"
        with open(filepath, 'wb') as f:
            pickle.dump(mass_index, f)
        # 记录每个区间的肽段数目

        # 把所有肽段的质量列表转成Ctypes形式
        data_pep_mass = (c_double * len(self.list_pep_mass))()
        for peptide_index, mass in enumerate(self.list_pep_mass):
            data_pep_mass[peptide_index] = mass
        p_C_pep_mass = byref(data_pep_mass)
        # 把所有肽段的质量列表转成Ctypes形式

        # self.__captain_read_ion_index_data("F:\wjq0619\pkl\only_cross\ion_index_only_cross.txt")
        print("process spectrum")
        self.__captain_flatten_spectrum_data(spectrum_data)
        # 搜索谱图
        get_index_time = 0
        search_time = 0
        filter_time = 0
        get_mass_time = 0
        get_candidate_time = 0
        write_file_time = 0
        for res_index, one_spe in enumerate(spectrum_data):
            if (res_index + 1) % 2000 == 0:
                print("[Info] Part : %d" % (part_index + 1) + " Finished crosslink coarse search : %.2f%%" % (
                        (res_index + 1) / len(spectrum_data) * 100) + " %d" % (res_index + 1) + " / %d" % (
                          len(spectrum_data)))

                if (res_index + 1) % 20000 == 0:
                    with open(fr"F:\wjq0619\result\coarse_res_{res_index}", 'wb') as f:
                        pickle.dump(coarse_res, f)
                # 打印耗时
                print("Time: {}\t{}\t{}\t{}\t{}\t{}"
                      .format(get_index_time, search_time, filter_time, get_mass_time, get_candidate_time,
                              write_file_time))
                get_index_time = 0
                search_time = 0
                filter_time = 0
                get_mass_time = 0
                get_candidate_time = 0
                write_file_time = 0

            if one_spe is None:
                continue
            if res_index > self.dp.myCFG.D4_NUMBER_SPECTRUM:
                break

            one_precursor_moz = one_spe[0]
            one_precursor_charge = one_spe[1]
            one_precursor_mass = one_precursor_moz - self.dp.myINI.MASS_PROTON_MONO
            # 谱图的谱峰信息
            no_linker_mz = one_spe[2]
            no_linker_intensity = one_spe[3]
            linker_mz = one_spe[4]
            linker_intensity = one_spe[5]

            # 不对强度加和所以原始的就是最高的
            raw_max_inten = 0.0
            for one_inten in no_linker_intensity:
                if one_inten > raw_max_inten:
                    raw_max_inten = one_inten
            if raw_max_inten <= 0:
                continue
            # 匹配的肽段不会超过这个肽段质量列表里的这个位置
            # 因为这个moz多了1.007，所以比允许的误差大得多了就不用算那么仔细了
            start_time = time.time()
            index_end_candidate_peptide = tool_binary_search(self.list_pep_mass,
                                                             one_precursor_moz - self.dp.myCFG.D6_MASS_PEP_LOW - self.linker_mass)
            max_candidate_pep_num = index_end_candidate_peptide + 1
            end_time = time.time()
            get_index_time += end_time - start_time

            # 索引6是最终分数
            # print("idx:{}: moz: {}\t max_pep: {}".format(res_index, one_precursor_moz, index_end_candidate_peptide))
            start_time = time.time()
            p_C_record_list = bird.malloc_record_score_list(max_candidate_pep_num, 7)
            # print("spectrum idx: {}, candidate: {}".format(res_index, max_candidate_pep_num))

            max_score_index = self.__captain_get_peptide_match_score(no_linker_mz, no_linker_intensity, linker_mz,
                                                                     linker_intensity, raw_max_inten,
                                                                     max_candidate_pep_num, one_precursor_charge, bird,
                                                                     p_C_record_list)
            end_time = time.time()
            search_time += end_time - start_time
            limit_score = p_C_record_list[max_candidate_pep_num * 6 + max_score_index] / 2
            # print("max_score_index: {}\tlimit_score: {}".format(max_score_index, limit_score))
            # 保留有匹配的结果

            start_time = time.time()
            valid_num = c_long(0)
            p_valid_index = bird.malloc_long_p(max_candidate_pep_num, 1)
            p_valid_score = bird.malloc_record_score_list(max_candidate_pep_num, 1)
            bird.filter_record_score_list(p_C_record_list, max_candidate_pep_num, 2, 6, c_double(1),
                                          c_double(limit_score), p_valid_score, p_valid_index, byref(valid_num))
            # if res_index == 130:
            #     print("matched: {}".format(p_C_record_list[max_candidate_pep_num * 2 + 0]))

            end_time = time.time()
            filter_time += end_time - start_time
            # 获取互补的肽段质量
            start_time = time.time()
            p_valid_search_left_mass = bird.malloc_record_score_list(max_candidate_pep_num, 1)
            p_valid_search_right_mass = bird.malloc_record_score_list(max_candidate_pep_num, 1)
            # list_pep_mass done
            #
            bird.get_peptide_mass(p_C_pep_mass, p_valid_search_left_mass, p_valid_search_right_mass, p_valid_index,
                                  valid_num, self.is_ppm_pre, self.precursor_ppm, c_double(one_precursor_mass),
                                  c_double(self.linker_mass))
            end_time = time.time()
            get_mass_time += end_time - start_time

            # 记录候结果信息
            candidate_crosslink_index_list = []
            candidate_crosslink_score_list = []
            candidate_crosslink_mass_list = []
            candidate_crosslink_separate_score_list = []
            # 找候选交联肽段

            start_time = time.time()
            if valid_num.value < 1:
                bird.free_long_p(p_valid_index)
                bird.free_record_score_list(p_valid_score)
                bird.free_record_score_list(p_valid_search_left_mass)
                bird.free_record_score_list(p_valid_search_right_mass)

                bird.free_record_score_list(p_C_record_list)

            else:
                candidate_pep_index = p_valid_index[0: valid_num.value]
                candidate_pep_BM25_score = p_valid_score[0: valid_num.value]
                candidate_pep_search_left_mass = p_valid_search_left_mass[0: valid_num.value]
                candidate_pep_search_right_mass = p_valid_search_right_mass[0: valid_num.value]

                # 释放内存
                bird.free_long_p(p_valid_index)
                bird.free_record_score_list(p_valid_score)
                bird.free_record_score_list(p_valid_search_left_mass)
                bird.free_record_score_list(p_valid_search_right_mass)

                list_BM25_score = p_C_record_list[max_candidate_pep_num * 6: max_candidate_pep_num * 7]
                bird.free_record_score_list(p_C_record_list)

                self.__captain_get_candidate_peptide(
                    candidate_pep_index,
                    candidate_pep_BM25_score,
                    candidate_pep_search_left_mass,
                    candidate_pep_search_right_mass,
                    mass_index,
                    index_end_candidate_peptide,
                    list_BM25_score,
                    max_score_index,
                    candidate_crosslink_index_list,
                    candidate_crosslink_score_list,
                    candidate_crosslink_mass_list,
                    candidate_crosslink_separate_score_list,
                    res_index
                )
            if res_index == 62462:
                print(1)

            end_time = time.time()
            get_candidate_time += end_time - start_time

            start_time = time.time()
            if len(candidate_crosslink_score_list) == 0:
                pass
            else:
                #  第一条肽是按分数排序，把它们换成按质量排序，从小到大顺序，索引号是质量排第几，值是对应pep_1_list第几个，
                #  pep_1_list[pep_1_sort_index_list[0]]  质量最小的那个肽段它在整个肽段列表的位置
                # 这里是把alpha肽段对应序列相同位点不同的去掉
                candidate_crosslink_score_sort_index_list = np.argsort(candidate_crosslink_score_list)
                get_index = len(candidate_crosslink_score_sort_index_list) - 1
                get_num = 0
                while get_index >= 0:

                   coarse_res[res_index].append(
                        (candidate_crosslink_index_list[candidate_crosslink_score_sort_index_list[get_index]][0],
                         candidate_crosslink_index_list[candidate_crosslink_score_sort_index_list[get_index]][1],
                         candidate_crosslink_separate_score_list[candidate_crosslink_score_sort_index_list[get_index]]))
                   self.list_pepcode2resultindex[
                        candidate_crosslink_index_list[candidate_crosslink_score_sort_index_list[get_index]][0]] = 1
                   self.list_pepcode2resultindex[
                        candidate_crosslink_index_list[candidate_crosslink_score_sort_index_list[get_index]][1]] = 1
                   get_num += 1
                   get_index -= 1
                   if get_num > 100:
                        break
            end_time = time.time()
            write_file_time += end_time - start_time
            finish_spectrum_search_num += 1
        time_2 = time.time()
        print("[Info] Part : %d" % (part_index + 1) + " Finished crosslink coarse search : 100%")
        print("[Info] Part : %d" % (part_index + 1) + " Using time : %f s" % (time_2 - time_1))

    def __captain_search_by_cuda(self):
        pass

    def __captain_flatten_spectrum_data(self, spectrum_data):
        list_precursor = []
        charge = []
        no_linker_mz = []
        no_linker_mz_prefix = [0]
        no_linker_intensity = []
        no_linker_intensity_prefix = [0]
        linker_mz = []
        linker_mz_prefix = [0]
        linker_intensity = []
        linker_intensity_prefix = [0]
        file_dir = "F:\wjq0619\pkl\spectrum2\\"
        file_path = file_dir + "spectrum.npz"
        for res_index, one_spe in enumerate(spectrum_data):
            list_precursor.append(one_spe[0])
            charge.append(one_spe[1])
            no_linker_mz += one_spe[2]
            no_linker_mz_prefix.append(no_linker_mz_prefix[-1] + len(one_spe[2]))
            no_linker_intensity += one_spe[3]
            no_linker_intensity_prefix.append(no_linker_intensity_prefix[-1] + len(one_spe[3]))
            linker_mz += one_spe[4]
            linker_mz_prefix.append(linker_mz_prefix[-1] + len(one_spe[4]))
            linker_intensity += one_spe[5]
            linker_intensity_prefix.append(linker_intensity_prefix[-1] + len(one_spe[5]))

        print("save file : {}".format(file_path))
        np.savez(file_path, no_linker_mz=np.array(no_linker_mz, dtype=np.int32),
                 no_linker_mz_prefix=np.array(no_linker_mz_prefix, dtype=np.int32),
                 no_linker_intensity=np.array(no_linker_intensity, dtype=np.float32),
                 no_linker_intensity_prefix=np.array(no_linker_intensity_prefix, dtype=np.int32),
                 linker_mz=np.array(linker_mz, dtype=np.int32),
                 linker_mz_prefix=np.array(linker_mz_prefix, dtype=np.int32),
                 linker_intensity=np.array(linker_intensity, dtype=np.float32),
                 linker_intensity_prefix=np.array(linker_intensity_prefix, dtype=np.int32),
                 list_precursor=np.array(list_precursor, dtype=np.float64),
                 charge=np.array(charge, dtype=np.int32))
    def __captain_read_ion_index_data(self, path):
        ivf_dic = []
        with open(path, 'r') as f:
            for line in f:
                values = line.strip().split('\t')
                pep_num = values[0]
                ivf_dic.append(values[1:])
        flatten_and_dump(ivf_dic, path.split('.')[0] + '_bin.npz')

    def __captain_get_peptide_match_score(self,
                                          no_linker_mz,
                                          no_linker_intensity,
                                          linker_mz,
                                          linker_intensity,
                                          raw_max_inten,
                                          max_candidate_pep_num,
                                          one_precursor_charge,
                                          bird,
                                          p_C_record_list):

        # 每一张谱图都要是C能读的，要转换一下
        spe_no_linker_mz = (c_long * len(no_linker_mz))()
        spe_no_linker_inten = (c_double * len(no_linker_intensity))()
        for mz_index, mz in enumerate(no_linker_mz):
            spe_no_linker_mz[mz_index] = mz
            spe_no_linker_inten[mz_index] = no_linker_intensity[mz_index] / raw_max_inten
        p_C_spe_no_linker_mz = byref(spe_no_linker_mz)
        p_C_spe_no_linker_inten = byref(spe_no_linker_inten)

        # 每一张谱图都要是C能读的，要转换一下
        spe_linker_mz = (c_long * len(linker_mz))()
        spe_linker_inten = (c_double * len(linker_intensity))()
        for mz_index, mz in enumerate(linker_mz):
            spe_linker_mz[mz_index] = mz
            spe_linker_inten[mz_index] = linker_intensity[mz_index] / raw_max_inten
        p_C_spe_linker_mz = byref(spe_linker_mz)
        p_C_spe_linker_inten = byref(spe_linker_inten)

        max_candidate_score = c_double(0.0)
        max_candidate_score_index = -1

        # step_result = []
        # 不跨交联匹配数
        # print("[MAIN] bird.get_ion_index_match_num1")
        bird.get_ion_index_match_num(p_C_spe_no_linker_mz,
                                     len(no_linker_intensity),
                                     self.p_C_ion_index_start,
                                     self.dp.myCFG.max_mz_index_size,
                                     max_candidate_pep_num,
                                     self.is_ppm_pre,
                                     self.fragment_ppm,
                                     p_C_record_list,
                                     0)
        # max_score_index = bird.find_max_value_index(p_C_record_list, max_candidate_pep_num, 0)
        # # print(p_C_record_list[max_score_index])
        # # # step_result.append([max_score_index, p_C_record_list[max_candidate_pep_num * 0 + max_score_index]])
        # if 11080 == max_candidate_pep_num:
        #     print(max_score_index, p_C_record_list[max_score_index])
        #     # step_result.append([9436, p_C_record_list[max_candidate_pep_num * 0 + 9436]])
        # # # print(p_C_record_list[max_score_index])
        # # 跨交联匹配数
        # print("[MAIN] bird.get_ion_index_match_num2")
        bird.get_ion_index_match_num(p_C_spe_linker_mz,
                                     len(linker_mz),
                                     self.p_C_ion_index_start,
                                     self.dp.myCFG.max_mz_index_size,
                                     max_candidate_pep_num,
                                     self.is_ppm_pre,
                                     self.fragment_ppm,
                                     p_C_record_list,
                                     1)
        # max_score_index = bird.find_max_value_index(p_C_record_list, max_candidate_pep_num, 1)
        # if 11080 == max_candidate_pep_num:
        #     print(max_score_index, p_C_record_list[max_candidate_pep_num * 1 + max_score_index])
        # # # step_result.append([max_score_index, p_C_record_list[max_candidate_pep_num * 1 + max_score_index]])
        # # if 9436 < max_candidate_pep_num:
        # #     step_result.append([9436, p_C_record_list[max_candidate_pep_num * 1 + 9436]])
        # # 两者得分求和
        # print("[MAIN] bird.get_sum_score")

        bird.get_sum_score(p_C_record_list, max_candidate_pep_num, 0, 1, 2,
                           max_candidate_score, max_candidate_score_index)

        # print("[MAIN] bird.find_max_value_index")
        # max_score_index = bird.find_max_value_index(p_C_record_list, max_candidate_pep_num, 2)
        # #
        # if 11080 == max_candidate_pep_num:
        #     print(max_score_index, p_C_record_list[max_candidate_pep_num * 2 + max_score_index])
        # # # step_result.append([max_score_index, p_C_record_list[max_candidate_pep_num * 2 + max_score_index]])
        # if 9436 < max_candidate_pep_num:
        #     step_result.append([9436, p_C_record_list[max_candidate_pep_num * 2 + 9436]])
        # 匹配比例
        # print(self.p_C_list_ion_num[])
        # print("[MAIN] bird.get_ion_index_match_num_percent")
        bird.get_ion_index_match_num_percent(self.p_C_list_ion_num, one_precursor_charge, p_C_record_list, 2, 3,
                                             max_candidate_pep_num)
        max_score_index = bird.find_max_value_index(p_C_record_list, max_candidate_pep_num, 3)
        # if 11080 == max_candidate_pep_num:
        #     print(2513, p_C_record_list[max_candidate_pep_num * 3 + 2513])
        # # # step_result.append([max_score_index, p_C_record_list[max_candidate_pep_num * 3 + max_score_index]])
        # # if 9436 < max_candidate_pep_num:
        # #     step_result.append([9436, p_C_record_list[max_candidate_pep_num * 3 + 9436]])
        # # 不跨交联得分
        # print("[MAIN] bird.get_ion_index_match_score")
        bird.get_ion_index_match_score(p_C_spe_no_linker_mz,
                                       len(no_linker_mz),
                                       p_C_spe_no_linker_inten,
                                       self.p_C_ion_index_start,
                                       self.dp.myCFG.max_mz_index_size,
                                       max_candidate_pep_num,
                                       self.is_ppm_pre,
                                       self.fragment_ppm,
                                       p_C_record_list,
                                       3,
                                       4)
        # max_score_index = bird.find_max_value_index(p_C_record_list, max_candidate_pep_num, 4)
        # if 11080 == max_candidate_pep_num:
        #     print(max_score_index, p_C_record_list[max_candidate_pep_num * 4 + max_score_index])
        # # # step_result.append([max_score_index, p_C_record_list[max_candidate_pep_num * 4 + max_score_index]])
        # # if 9436 < max_candidate_pep_num:
        # #     step_result.append([9436, p_C_record_list[max_candidate_pep_num * 4 + 9436]])
        # # 跨交联的得分
        # print("[MAIN] bird.get_ion_index_match_score")
        bird.get_ion_index_match_score(p_C_spe_linker_mz,
                                       len(linker_mz),
                                       p_C_spe_linker_inten,
                                       self.p_C_ion_index_start,
                                       self.dp.myCFG.max_mz_index_size,
                                       max_candidate_pep_num,
                                       self.is_ppm_pre,
                                       self.fragment_ppm,
                                       p_C_record_list,
                                       3,
                                       5)
        # max_score_index = bird.find_max_value_index(p_C_record_list, max_candidate_pep_num, 5)
        # if 11080 == max_candidate_pep_num:
        #     print(max_score_index, p_C_record_list[max_candidate_pep_num * 5 + max_score_index])
        # # # step_result.append([max_score_index, p_C_record_list[max_candidate_pep_num * 5 + max_score_index]])
        # # if 9436 < max_candidate_pep_num:
        # #     step_result.append([9436, p_C_record_list[max_candidate_pep_num * 5 + 9436]])
        # # 最终得分
        #
        # print("[MAIN] bird.get_sum_score")
        bird.get_sum_score(p_C_record_list, max_candidate_pep_num, 4, 5, 6, max_candidate_score,
                           max_candidate_score_index)
        # max_score_index = bird.find_max_value_index(p_C_record_list, max_candidate_pep_num, 6)
        # if 11080 == max_candidate_pep_num:
        #     print(max_score_index, p_C_record_list[max_candidate_pep_num * 6 + max_score_index])
        # # step_result.append([max_score_index, p_C_record_list[max_candidate_pep_num * 6 + max_score_index]])
        # if 9436 < max_candidate_pep_num:
        #     step_result.append([9436, p_C_record_list[max_candidate_pep_num * 6 + 9436]])
        # formatted_data = ["{},{}".format(pair[0], pair[1]) for pair in step_result]
        # print("\t".join(formatted_data))
        # 找到最高分
        max_score_index = bird.find_max_value_index(p_C_record_list, max_candidate_pep_num, 6)
        # score_list = p_C_record_list[max_candidate_pep_num * 6: max_candidate_pep_num * 7]
        # match_no_link_num = p_C_record_list[: max_candidate_pep_num * 1]
        # match_link_num = p_C_record_list[max_candidate_pep_num * 1: max_candidate_pep_num * 2]
        return max_score_index

    # def get_ion_index_match_num(self, mz_list, intensity_list, ):

    # def read_ion_index_bin(self, read_path, input_mz_list_len):

    def __captain_get_candidate_peptide(self,
                                        candidate_pep_index,
                                        candidate_pep_BM25_score,
                                        candidate_pep_search_left_mass,
                                        candidate_pep_search_right_mass,
                                        mass_index,
                                        index_end_candidate_peptide,
                                        list_BM25_score,
                                        max_score_index,
                                        candidate_crosslink_index_list,
                                        candidate_crosslink_score_list,
                                        candidate_crosslink_mass_list,
                                        candidate_crosslink_separate_score_list,
                                        res_idx=None):
        for pep_1_index, pep_1_mass_index in enumerate(candidate_pep_index):
            # if pep_1_mass_index == 1588:
            #     print(1)
            pep_1_score = candidate_pep_BM25_score[pep_1_index]
            pep_1_link_site_type = self.list_link_site[pep_1_mass_index][1]

            pep_2_left_mass = candidate_pep_search_left_mass[pep_1_index]
            pep_2_right_mass = candidate_pep_search_right_mass[pep_1_index]

            if pep_2_left_mass > self.dp.myCFG.D7_MASS_PEP_UP:
                continue
            if pep_2_right_mass < self.dp.myCFG.D6_MASS_PEP_LOW:
                continue
            if pep_2_left_mass < self.dp.myCFG.D6_MASS_PEP_LOW:
                pep_2_left_mass = self.dp.myCFG.D6_MASS_PEP_LOW
            if pep_2_right_mass > self.dp.myCFG.D7_MASS_PEP_UP:
                pep_2_left_mass = self.dp.myCFG.D7_MASS_PEP_UP
            pep_2_left_mass_index = mass_index[int(pep_2_left_mass) - 1]
            pep_2_right_mass_index = mass_index[int(pep_2_right_mass)]
            if pep_2_left_mass_index >= pep_2_right_mass_index:
                continue

            # 查找候选的互补肽段所在整个肽段列表中的位置
            pep_2_mass_index_list_start_index = tool_binary_search_index(
                self.list_pep_mass[pep_2_left_mass_index:pep_2_right_mass_index], pep_2_left_mass)

            if pep_2_mass_index_list_start_index[0] == index_end_candidate_peptide:
                continue
            pep_2_mass_index = pep_2_mass_index_list_start_index[1] + pep_2_left_mass_index
            # 切片获取右边超出不会影响所以只做左边的-1处理
            if pep_2_mass_index == -1:
                pep_2_mass_index = 0
            while True:
                if pep_2_mass_index >= len(self.list_pep_mass):
                    break
                if self.list_pep_mass[pep_2_mass_index] > pep_2_right_mass:
                    break
                pep_2_score = list_BM25_score[pep_2_mass_index]
                pep_2_link_site_type = self.list_link_site[pep_2_mass_index][1]
                if pep_1_score + pep_2_score < list_BM25_score[max_score_index]:
                    pass
                else:
                    if pep_1_link_site_type + pep_2_link_site_type != 0:
                        pass
                    else:
                        # if res_idx == 130:
                        #     print("1: {}, 2: {}".format(pep_1_mass_index, pep_2_mass_index))
                        candidate_crosslink_index_list.append((pep_1_mass_index, pep_2_mass_index))
                        candidate_crosslink_score_list.append(pep_1_score + pep_2_score)
                        candidate_crosslink_separate_score_list.append((pep_1_score, pep_2_score))
                        candidate_crosslink_mass_list.append(
                            (self.list_pep_mass[pep_1_mass_index], self.list_pep_mass[pep_2_mass_index]))
                pep_2_mass_index += 1
