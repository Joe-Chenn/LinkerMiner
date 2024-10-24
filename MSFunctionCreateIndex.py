# -*- mode: python ; coding: utf-8 -*-..
import sys

from MSData import CModSite, CPeptidePKL
from MSDataC import ION_INDEX_MZ
from MSFunction import CFunctionPickle
from MSFunctionMutiProcess import CFunctionMultiProcessSortPKL
from MSSysterm import DLL_INDEX, DLL_FUNCTION
from MSSysterm import FOLDER_MASS, AA_NUM, VRA_PKL_FILE, TMP_FOLDER
from MSSysterm import FOLDER_INDEX, FILE_NAME_ION_INDEX, FILE_NAME_PEP_DATA, FILE_NAME_ION_NUM, FILE_NAME_PEP_MASS, FILE_NAME_PEP_PKL, FILE_NAME_LINK_SITE
from MSSysterm import PEP_PRO_INDEX_FILE
from MSTool import tool_create_aa_list
from MSTool_code import tool_get_peptide_info_code, tool_get_peptide_info, tool_get_fix_mod_site_list, tool_get_var_mod_data
from MSLogging import logGetError
from MSFunctionCreatePeptide import CFunctionAddSite

from MSOperator import op_get_fix_mod_site_list

import psutil
import time
import math
import os
import pickle
from ctypes import *

class CFunctionCreatePepIndex:

    def __init__(self, inputDP):

        self.dp = inputDP

    def generate_peptide_index_file(self, mutiprocess_pkl_file_list, del_old=True):
        # # 建立质量区间pkl文件
        if del_old:
            self.__captain_delete_old_pkl()  # 删除旧pkl文件
        self.__captain_get_mass_pkl(mutiprocess_pkl_file_list)  # 生成新pk文件

    def __captain_delete_old_pkl(self):
        if os.path.exists(self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_MASS + '\\'):
            for p in os.listdir(self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_MASS + '\\'):
                if p.endswith(".pkl"):
                    os.remove(os.path.join(self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_MASS + '\\', p))  # 删除该路径下存在的pkl文件
        else:
            os.makedirs(self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_MASS + '\\')

    def __captain_get_mass_pkl(self, mutiprocess_pkl_file_list):

        # 默认按100Da区间建pkl文件
        for i in range(len(self.dp.LIST_MASS_INDEX_FILE)):
            pkl_file = self.dp.LIST_MASS_INDEX_FILE[i]
            ind_file = self.dp.LIST_MASS_IND_FILE[i]
            f = open(pkl_file, 'wb')
            f.close()
            f = open(ind_file, 'wb')
            f.close()
            mutiprocess_pkl_file_list.append((pkl_file, ind_file))

    def create_peptide_index(self, protein_list, peptide_matrix, peptide_ind_matrix):

        time_1 = time.time()

        # 现对固定修饰去重，再加可变修饰
        memory_data = psutil.virtual_memory()
        before_used_memory = float(memory_data.used) / 1024 / 1024 / 1024

        peptide_dataset = []
        functionAddSite = CFunctionAddSite(self.dp)
        functionAddSite.create_fix_peptide_list(protein_list, peptide_dataset)

        set_peptide_dataset = []
        peptide2protein_index = []  # 记录肽段来自的蛋白，因为可变修饰和固定修饰肽段来源会有相同的，去冗余存储
        # 去重肽段索引
        time_2 = time.time()
        self.__captain_sort_fix_peptide_list_data(protein_list, peptide_dataset, set_peptide_dataset, peptide2protein_index)
        self.__captain_output_peptide_protein_index(tuple(peptide2protein_index))

        time_3 = time.time()
        print("[Info] Sort fix mod peptides time : %.4f s . "%(time_3 - time_2) )
        print("before sort : %d"%(len(peptide_dataset)))
        print("after sort : %d" % (len(set_peptide_dataset)))

        memory_data = psutil.virtual_memory()
        after_used_memory = float(memory_data.used) / 1024 / 1024 / 1024
        print("after_used_memory : %.4f GB" % (after_used_memory))
        print("used memory : %.4f GB"% (after_used_memory - before_used_memory))

        del peptide_dataset
        del peptide2protein_index

        self.__captain_init_peptide_matrix(peptide_matrix, False)
        self.__captain_init_peptide_matrix(peptide_ind_matrix, True)

        var_mod_peptide = []
        functionAddSite.generate_var_peptide_list(set_peptide_dataset, protein_list, var_mod_peptide)

        '''
        out_data = []
        for one_peptide in var_mod_peptide:
            is_target, protein_index, peptide_start, peptide_length = tool_get_peptide_info(one_peptide[1])
            mod_site_list = []
            tool_get_fix_mod_site_list(one_peptide[2], mod_site_list, peptide_length)
            pep_sq = protein_list.sq[protein_index][peptide_start: peptide_start + peptide_length]
            if len(one_peptide) == 6:
                tool_get_var_mod_data(one_peptide[5], mod_site_list)
            gdm = self.op_get_peptide_gdm(mod_site_list, pep_sq, peptide_length, 1, self.dp.G_matrix, AA_NUM)
            out_data.append([pep_sq, mod_site_list, one_peptide[0], gdm])
        out_data.sort(key=lambda x: (x[2], x[3]), reverse=False)

        with open("E:\\wqsh\\code_test\\eLink\\pLink_data\\eLink\\test_human\\new_peptides.txt", 'w') as f:
            for item in out_data:
                f.write(item[0])
                f.write('\t')
                # f.write(str(protein_list.ac[protein_index]))
                # f.write('\t')
                f.write(str(item[1]))
                f.write('\t')
                f.write(str(item[2]))
                f.write('\t')
                f.write(str(item[3]))
                f.write('\n')
        exit(0)
        '''

        self.__captain_sort_var_peptide_list_data(var_mod_peptide)
        self.__captain_separate_peptide(var_mod_peptide, peptide_matrix)
        self.__captain_merge_peptide_index(peptide_matrix, peptide_ind_matrix)

        # for one_var_mod_peptide in var_mod_peptide:
        #     self.dp.myCFG.

        # var_num = self.__soldier_sort_peptide_list_data(protein_list, var_mod_peptide, final_var_mod_peptide, set_data=False)
        # print("var_num : " + str(var_num))
        # 31860035 可变修饰肽段总数，总肽段数39470157
        # functionAddSite.create_var_peptide_new(set_peptide_dataset, peptide_matrix, protein_list, tmp)
        # functionAddSite.create_var_peptide(set_peptide_dataset, peptide_matrix, protein_list)
        # self.__captain_sort_peptide_matrix(protein_list, peptide_matrix, peptide_ind_matrix, set_data=False)

        '''
        var_mod_peptide.sort(key=lambda x: (x[0], x[4]), reverse=False)
        with open(self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\new_peptide.txt', 'w') as f:
            for index in range(len(var_mod_peptide)):
                one_peptide = var_mod_peptide[index]
                target, pro_index, pep_start_pos, pep_length = tool_get_peptide_info(one_peptide[1])
                peptide_sq = protein_list.sq[pro_index][pep_start_pos: pep_start_pos + pep_length]
                mod_site_list = []
                tool_get_fix_mod_site_list(one_peptide[2], mod_site_list, pep_length)
                if len(one_peptide) < 6:
                    pass
                else:
                    tool_get_var_mod_data(one_peptide[6], mod_site_list)
                f.write(peptide_sq)
                f.write('\t')
                f.write('%.6f'%one_peptide[0])
                f.write('\t')
                f.write('%.6f'%one_peptide[4])
                f.write('\t')
                f.write(str(mod_site_list))
                f.write('\n')
        '''
        # for index in range(len(peptide_matrix)):
        #     for one_peptide in peptide_matrix[index]:
        #         target, pro_index, pep_start_pos, pep_length = tool_get_peptide_info(one_peptide.pep_code)
        #         peptide_sq = protein_list.sq[pro_index][pep_start_pos: pep_start_pos + pep_length]
        #         if peptide_sq == 'GGGSK':
        #             record_protein_data = [protein_list.ac[pro_index]]
        #             for peptide_code in one_peptide.pep_code_list:
        #                 target, pro_index, pep_start_pos, pep_length = tool_get_peptide_info(peptide_code)
        #                 record_protein_data.append(protein_list.ac[pro_index])
        #             print('aa')
        #
        # with open(self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\mid_peptide.txt', 'w') as f:
        #     for index in range(len(peptide_matrix)):
        #         for one_peptide in peptide_matrix[index]:
        #             target, pro_index, pep_start_pos, pep_length = tool_get_peptide_info(one_peptide.pep_code)
        #             peptide_sq = protein_list.sq[pro_index][pep_start_pos: pep_start_pos + pep_length]
        #             mod_site_list = []
        #             tool_get_fix_mod_site_list(one_peptide.fix_mod_code, mod_site_list, pep_length)
        #             tool_get_var_mod_data(one_peptide.var_mod_code, mod_site_list)
        #             f.write(peptide_sq)
        #             f.write('\t')
        #             f.write('%.6f'%one_peptide.mass)
        #             f.write('\t')
        #             f.write('%.6f'%one_peptide.gdm)
        #             f.write('\t')
        #             f.write(str(mod_site_list))
        #             f.write('\n')

        memory_data = psutil.virtual_memory()
        after_used_memory = float(memory_data.used) / 1024 / 1024 / 1024
        print("before_used_memory : %.4f GB" % (before_used_memory))
        print("after_used_memory : %.4f GB" % (after_used_memory))
        print("used memory : %.4f GB"% (after_used_memory - before_used_memory))

        time_2 = time.time()
        print("[Info] Generate candidate peptides time : %.4f s."%(time_2 - time_1))

    def __debug_get_peptide_gdm(self, mod_site_list, sq, pep_len, fix_mod_num, G_matrix, aa_num):
        # 肽段序列，肽段修饰list，存着这种情况的修饰，一行是一个CModSite类,类的内容包括mod_name, site, mass
        # get the gdm value of peptide sequence with variable modifications
        # The same peptide sequence with same modifications and sites have the same gdm values
        gdm_value = 0.0
        for aa_index, mod_index in enumerate(mod_site_list):
            if aa_index == 0:
                if mod_site_list[aa_index] >= fix_mod_num:
                    gdm_value += G_matrix[aa_num + mod_site_list[aa_index] - fix_mod_num][aa_index]
                else:
                    gdm_value += 0
            elif aa_index == pep_len + 1:
                if mod_site_list[aa_index] >= fix_mod_num:
                    gdm_value += G_matrix[aa_num + mod_site_list[aa_index] - fix_mod_num][aa_index]
                else:
                    gdm_value += 0
            else:
                if mod_site_list[aa_index] >= fix_mod_num:
                    gdm_value += G_matrix[aa_num + mod_site_list[aa_index] - fix_mod_num][aa_index]
                else:
                    gdm_value += G_matrix[ord(sq[aa_index - 1]) - 65][aa_index]
        return gdm_value  # 返回肽段对应的哥德尔编码值

    def __captain_init_peptide_matrix(self, peptide_matrix, list_none=False):
        for i in range(len(self.dp.LIST_MASS_INDEX_FILE)):
            if list_none:
                peptide_matrix.append(None)
            else:
                peptide_matrix.append([])

    def __captain_sort_peptide_matrix(self, protein_list, peptide_matrix, peptide_ind_matrix, set_data=True):
        peptide_num = 0
        for peptide_mass_index in range(len(peptide_matrix)):
            time_1 = time.time()
            original_peptide_list = peptide_matrix[peptide_mass_index]
            sorted_peptide_list = []
            peptide_num += self.__soldier_sort_peptide_data(protein_list, original_peptide_list, sorted_peptide_list, set_data)
            # del original_peptide_list
            peptide_matrix[peptide_mass_index] = sorted_peptide_list
            start_mass = self.dp.myCFG.start_mass + peptide_mass_index * 100
            end_mass = self.dp.myCFG.start_mass + (peptide_mass_index + 1) * 100
            peptide_ind_matrix[peptide_mass_index] = self.__soldier_create_range_index(sorted_peptide_list, start_mass, end_mass, 1)
            time_2 = time.time()
            print("[Info] Sort peptides in " +
                  self.dp.LIST_MASS_INDEX_FILE[peptide_mass_index].split('\\')[-1].split('.')[0] + ' using time : %.4f s'%(time_2 - time_1))
        print("[Info] Sorted peptides num : %d"% peptide_num)

    def __captain_sort_fix_peptide_list_data(self, protein_list, peptide_list_data, new_all_peptides, peptide2protein_index):
        peptide_num = 0
        peptide_list_data.sort(key=lambda x: (x[0]), reverse=False)
        # 按gdm值排序去冗余
        start_i = 0  # 索引all_peptides这个list，相当于索引肽段
        while start_i < len(peptide_list_data):
            is_target = False
            is_pro_N = False
            is_pro_C = False

            new_all_peptides.append(peptide_list_data[start_i])
            new_all_peptides[-1].append(peptide_num)
            # 一个新的肽段信息
            cur_p = peptide_list_data[start_i]
            target, pro_index, pep_start_pos, pep_length = tool_get_peptide_info(peptide_list_data[start_i][1])
            cur_gdm = cur_p[0]  # 一个肽段哥德尔编码值
            if pep_start_pos == 0:
                is_pro_N = True
            if pep_length + pep_start_pos == len(protein_list.sq[pro_index]):
                is_pro_C = True
            # 因为排好序，下一个和这个不同那么后面都不会有与这个相同的了
            if target == 1:
                is_target = True
            # 下一个肽段索引
            record_same_sq = [peptide_list_data[start_i][1]]
            start_j = start_i + 1
            while start_j < len(peptide_list_data) and (peptide_list_data[start_j][0] == cur_gdm):
                target, pro_index, pep_start_pos, pep_length = tool_get_peptide_info(peptide_list_data[start_j][1])
                record_same_sq.append(peptide_list_data[start_j][1])
                # 下一个肽段的质量以及哥德尔编码值与上一个肽段值一样
                if is_target is True:
                    pass
                else:
                    if peptide_list_data[start_j][1] >> 63 == 1:
                        new_all_peptides[-1][1] = peptide_list_data[start_j][1]
                        is_target = True
                if not is_pro_N:
                    if pep_start_pos == 0 and target:
                        is_pro_N = True
                        new_all_peptides[-1][1] = peptide_list_data[start_j][1]
                if not is_pro_C:
                    if pep_length + pep_start_pos == len(protein_list.sq[pro_index]) and target:
                        is_pro_C = True
                        new_all_peptides[-1][1] = peptide_list_data[start_j][1]
                start_j += 1  # 再下一个肽段
            new_all_peptides[-1] = tuple(new_all_peptides[-1])
            peptide2protein_index.append(tuple(record_same_sq))
            del record_same_sq
            start_i = start_j
            peptide_num += 1
        # 整理好去重后的list
        return peptide_num

    def __captain_sort_var_peptide_list_data(self, peptide_list_data):

        peptide_list_data.sort(key=lambda x: (x[0]), reverse=False)
        peptide_num = len(peptide_list_data)
        return peptide_num

    def __captain_separate_peptide(self, peptide_list_data, peptide_matrix):
        start_mass = 0 + self.dp.myCFG.start_mass
        end_mass = start_mass + 100
        mass_index = 0
        record_peptide_list_index = 0
        for peptide_list_index, one_peptide_data in enumerate(peptide_list_data):
            peptide_mass = one_peptide_data[0]
            if peptide_mass > end_mass:
                end_mass += 100
                peptide_matrix[mass_index] = peptide_list_data[record_peptide_list_index: peptide_list_index]
                record_peptide_list_index = peptide_list_index
                mass_index += 1
        peptide_matrix[mass_index] = peptide_list_data[record_peptide_list_index:]

    def __captain_merge_peptide_index(self, peptide_matrix, peptide_ind_matrix, mul=1):
        for peptide_mass_index, peptide_list in enumerate(peptide_matrix):
            start_mass = self.dp.myCFG.start_mass + peptide_mass_index * 100
            end_mass = self.dp.myCFG.start_mass + (peptide_mass_index + 1) * 100
            num = int((end_mass - start_mass) * mul) + 1  # the number of list 区间按成0.01Da质量差进行分的份数
            index_list = [-1 for i in range(num)]  # index_list[X] is the min index in all_peptides whose mass is >= X
            for i, peptide_data in enumerate(peptide_list):  # 遍历一个质量区间的所有肽段，i表示第几个肽段
                index = int((peptide_data[0] - start_mass) * mul)  # 计算属于哪个区间的质量
                if index >= num:
                    index = num - 1
                if index_list[index] == -1:
                    index_list[index] = i
            # 这样得到的index_list表示这质量区间内，再分的小区间，其中对应的值是这个区间所在的质量在all_peptides中的右区间边界
            end_val = len(peptide_list)
            for i in range(num)[::-1]:
                if index_list[i] == -1:
                    index_list[i] = end_val  # index_list每一行表示到这个小区间已经有多少个肽段
                else:
                    end_val = index_list[i]
            peptide_ind_matrix[peptide_mass_index] = tuple(index_list)

    def __captain_output_peptide_protein_index(self, peptide_protein_index):
        path_peptide2protein_pkl_file = os.path.join(self.dp.myCFG.A4_PATH_FASTA_EXPORT, PEP_PRO_INDEX_FILE)
        functionPickle = CFunctionPickle()
        functionPickle.write_data_to_pkl(peptide_protein_index, path_peptide2protein_pkl_file)

    def __soldier_output_peptide_data(self, peptide_data, path_peptide_file, path_ind_file):

        functionPickle = CFunctionPickle()
        functionPickle.write_data_to_pkl(peptide_data, path_peptide_file)
        functionPickle.write_data_to_pkl(peptide_data, path_ind_file)

    def __soldier_sort_peptide_data(self, protein_list, peptide_list_data, new_all_peptides, set_data=True):
        peptide_num = 0
        peptide_list_data.sort(key=lambda x: (x.mass, x.gdm), reverse=False)
        # pkl_file = self.dp.LIST_MASS_INDEX_FILE[peptide_mass_index]
        # list 一行位置存一个类，表示一种修饰的肽段，这个已经按质量和哥德尔编码排好序了
        # 所属蛋白索引，肽段在蛋白中的起始位置，肽段在蛋白中的结束位置，肽段质量，修饰list修饰(一行是一个位点pkl类：修饰信息位置索引，位点# mod_index, site），哥德尔编码值
        # pro_index, start_pos, end_pos, mass, mods, gdm
        # mods : mod_index, site
        # all_peptides is the list of peptides sorted by masses and gdm values
        # return the list of peptides after removing the peptides with the same gdm values and update the protein indexes
        if set_data:
            start_i = 0  # 索引all_peptides这个list，相当于索引肽段
            while start_i < len(peptide_list_data):
                is_target = False
                is_pro_N = False
                is_pro_C = False
                new_all_peptides.append(peptide_list_data[start_i])
                peptide_num += 1
                cur_p = peptide_list_data[start_i]  # 一个肽段的信息
                cur_mass = cur_p.mass  # 一个肽段质量
                cur_gdm = cur_p.gdm  # 一个肽段哥德尔编码值
                start_j = start_i + 1  # 下一个肽段索引
                target, pro_index, pep_start_pos, pep_length = tool_get_peptide_info(peptide_list_data[start_i].pep_code)
                # pep_code_list = [peptide_list_data[start_i].pep_code]
                if pep_start_pos == 0:
                    is_pro_N = True
                if pep_length + pep_start_pos == len(protein_list.sq[pro_index]):
                    is_pro_C = True
                # 因为排好序，下一个和这个不同那么后面都不会有与这个相同的了
                if target == 1:
                    is_target = True
                while start_j < len(peptide_list_data) and (
                        peptide_list_data[start_j].mass == cur_mass) and (
                        peptide_list_data[start_j].gdm == cur_gdm):
                    target, pro_index, pep_start_pos, pep_length = tool_get_peptide_info(peptide_list_data[start_j].pep_code)
                    # 下一个肽段的质量以及哥德尔编码值与上一个肽段值一样
                    if is_target is True:
                        pass
                    else:
                        if peptide_list_data[start_j].pep_code >> 63 == 1:
                            new_all_peptides[-1].pep_code = peptide_list_data[start_j].pep_code
                            is_target = True
                    if not is_pro_N:
                        if pep_start_pos == 0 and target:
                            is_pro_N = True
                            new_all_peptides[-1].pep_code = peptide_list_data[start_j].pep_code
                    if not is_pro_C:
                        if pep_length + pep_start_pos == len(protein_list.sq[pro_index]) and target:
                            is_pro_C = True
                            new_all_peptides[-1].pep_code = peptide_list_data[start_j].pep_code
                    # pep_code_list.append(peptide_list_data[start_j].pep_code)  # 记录一样的肽段以及修饰一样的肽段来组不同的蛋白的索引号
                    start_j += 1  # 再下一个肽段
                # new_all_peptides[-1].pep_code_list = pep_code_list
                start_i = start_j
            # 整理好去重后的list
        else:
            for item in peptide_list_data:
                new_all_peptides.append(item)
                peptide_num += 1
        return peptide_num

    def __soldier_create_range_index(self, peptide_list, start_value, end_value, mul=1):
        # 去重后的的list，一行是一个肽段的类：所属蛋白索引，肽段在蛋白中的起始位置，肽段在蛋白中的结束位置，肽段质量，修饰list修饰(一行是一个位点pkl类：修饰信息位置索引，位点# mod_index, site），哥德尔编码值
        # 起始质量，结束质量，
        # 也按质量排好序了
        # all_peptides is the list of peptides sorted by masses and gdm values
        # return the list of mass index to the all_peptides
        num = int((end_value - start_value) * mul) + 1  # the number of list 区间按成0.01Da质量差进行分的份数
        index_list = [-1 for i in range(num)]  # index_list[X] is the min index in all_peptides whose mass is >= X
        for i in range(len(peptide_list)):  # 遍历一个质量区间的所有肽段，i表示第几个肽段
            p = peptide_list[i]  # 这个肽段的质量相对于区间起始质量的差*mul
            index = int((p.mass - start_value) * mul)  # 计算属于哪个区间的质量
            if index >= num:
                index = num - 1
            if index_list[index] == -1:
                index_list[index] = i
        # 这样得到的index_list表示这质量区间内，再分的小区间，其中对应的值是这个区间所在的质量在all_peptides中的右区间边界
        end_val = len(peptide_list)
        for i in range(num)[::-1]:
            if index_list[i] == -1:
                index_list[i] = end_val  # index_list每一行表示到这个小区间已经有多少个肽段
            else:
                end_val = index_list[i]
        return index_list

    def output_peptide_data(self, peptide_matrix, peptide_ind_matrix):

        functionPickle = CFunctionPickle()

        for peptide_matrix_index in range(len(peptide_matrix)):
            time_1 = time.time()
            peptide_index_name = self.dp.LIST_MASS_INDEX_FILE[peptide_matrix_index]
            peptide_ind_name = self.dp.LIST_MASS_IND_FILE[peptide_matrix_index]
            functionPickle.write_data_to_pkl(peptide_matrix[peptide_matrix_index], peptide_index_name)
            functionPickle.write_data_to_pkl(peptide_ind_matrix[peptide_matrix_index], peptide_ind_name)
            time_2 = time.time()
            print("[Info] Output peptide index : " + peptide_index_name + ' using time %.4f s'%(time_2 - time_1))

class CFunctionCreateOnePeptideIonIndex:

    peptide_type = 1
    link_type = 1

    def __init__(self, inputDP):

        self.dp = inputDP

        self.pep_num = 0
        self.used_pep_num_add_site = 0
        self.list_pep_data_pkl = [-1] * self.dp.myCFG.B9_MAX_PEPNUM  # 这个是记录满足交联条件的肽段索引，index->肽段编码
        self.list_ion_num = [-1] * self.dp.myCFG.B9_MAX_PEPNUM  # 这个是记录肽段的碎片离子数目，index->离子数
        self.list_pep_mass = [-1] * self.dp.myCFG.B9_MAX_PEPNUM  # 这个是记录肽段质量的，不包括交联剂质量
        self.list_link_site = [-1] * self.dp.myCFG.B9_MAX_PEPNUM

        self.link_1_pro_N = self.dp.myLINK.linker_1_data.same_site[-2] or self.dp.myLINK.linker_1_data.alpha_site[-2] or self.dp.myLINK.linker_1_data.beta_site[-2]
        self.link_1_pep_N = self.dp.myLINK.linker_1_data.same_site[-1] or self.dp.myLINK.linker_1_data.alpha_site[-1] or self.dp.myLINK.linker_1_data.beta_site[-1]
        self.link_1_pep_C = self.dp.myLINK.linker_1_data.same_site[1] or self.dp.myLINK.linker_1_data.alpha_site[1] or self.dp.myLINK.linker_1_data.beta_site[1]
        self.link_1_pro_C = self.dp.myLINK.linker_1_data.same_site[2] or self.dp.myLINK.linker_1_data.alpha_site[2] or self.dp.myLINK.linker_1_data.beta_site[2]

    # =================================================

    def import_DLL(self, bird_index, bird_function):

        # 导入动态链接库

        bird_index.malloc_ion_index.argtype = (c_int)
        bird_index.malloc_ion_index.restype = (POINTER(ION_INDEX_MZ))

        bird_index.free_ion_index.argtype = (c_int, POINTER(ION_INDEX_MZ))
        bird_index.free_ion_index.restype = (c_int)

        bird_index.put_peptide_to_ion_index.argtype = (ION_INDEX_MZ, c_double, c_int, c_int, c_long)
        bird_index.put_peptide_to_ion_index.restype = (c_int)

        bird_index.write_ion_index_bin.argtype = (c_char_p, POINTER(ION_INDEX_MZ), c_int)
        bird_index.write_ion_index_bin.restype = (c_int)

        # 导入动态链接库

        bird_function.tool_C_get_ion_code.argtype = (c_longlong, c_longlong, c_longlong, c_longlong)
        bird_function.tool_C_get_ion_code.restype = (c_longlong)

        bird_function.tool_C_new_get_ion_data.argtype = (c_longlong)
        bird_function.tool_C_new_get_ion_data.restype = POINTER(c_long)

        bird_function.tool_C_put_ion_to_ion_index_singlepeptide.argtype = (POINTER(ION_INDEX_MZ), c_char_p, c_int, POINTER(c_int), c_double, c_double, POINTER(c_double), POINTER(c_double), c_int, c_int, c_long)
        bird_function.tool_C_put_ion_to_ion_index_singlepeptide.restype = (c_int)

        bird_function.tool_C_put_ion_to_ion_index_crosslink.argtype = (POINTER(ION_INDEX_MZ), c_char_p, c_int, POINTER(c_int), c_double, c_double, POINTER(c_double), POINTER(c_double), c_int, c_int, c_int, c_long)
        bird_function.tool_C_put_ion_to_ion_index_crosslink.restype = (c_int)

        bird_function.tool_C_put_ion_to_ion_index_looplink.argtype = (POINTER(ION_INDEX_MZ), c_char_p, c_int, POINTER(c_int), c_double, c_double, POINTER(c_double), POINTER(c_double), c_int, c_int, c_double, c_int, c_int, c_long)
        bird_function.tool_C_put_ion_to_ion_index_looplink.restype = (c_int)

        bird_function.tool_C_put_ion_to_ion_index_crosscross.argtype = (POINTER(ION_INDEX_MZ), c_char_p, c_int, POINTER(c_int), c_double, c_double, POINTER(c_double), POINTER(c_double), c_int, c_int, c_int, c_int, c_long)
        bird_function.tool_C_put_ion_to_ion_index_crosscross.restype = (c_int)

        bird_function.tool_C_put_ion_to_ion_index_crossloop.argtype = (POINTER(ION_INDEX_MZ), c_char_p, c_int, POINTER(c_int), c_double, c_double, POINTER(c_double), POINTER(c_double), c_int, c_int, c_double, c_int, c_int, c_int, c_long)
        bird_function.tool_C_put_ion_to_ion_index_crossloop.restype = (c_int)

    def create_ion_index(self, protein_list, peptide_matrix=None, peptide_ind_matrix=None):

        # 导入动态链接库
        bird_index = CDLL(os.getcwd() + "\\" + DLL_INDEX)
        bird_function = CDLL(os.getcwd() + "\\" + DLL_FUNCTION)

        self.import_DLL(bird_index, bird_function)
        functionPickle = CFunctionPickle()
        p_mod_mass = self.__captain_C_get_p_mod_mass()
        p_aa_mass = self.__captain_C_get_p_aa_mass()
        p_C_ion_index = bird_index.malloc_ion_index(self.dp.myCFG.max_mz_index_size)
        time_1 = time.time()
        if self.dp.myCFG.D15_TYPE_LINK in [1]:
            # 每个质量区间遍历
            for peptide_list_index in range(len(self.dp.LIST_MASS_INDEX_FILE)):
                time_1 = time.time()
                if peptide_matrix is None:
                    peptide_list = functionPickle.load_pkl_to_data(self.dp.LIST_MASS_INDEX_FILE[peptide_list_index])
                    peptide_ind = functionPickle.load_pkl_to_data(self.dp.LIST_MASS_IND_FILE[peptide_list_index])
                else:
                    peptide_list = peptide_matrix[peptide_list_index]
                    peptide_ind = peptide_ind_matrix[peptide_list_index]
                self.create_single_peptide_ion_index(protein_list, peptide_list_index, peptide_list, peptide_ind, bird_function, p_mod_mass, p_aa_mass, p_C_ion_index)
                if peptide_matrix is not None:
                    peptide_matrix[peptide_list_index] = None
                    peptide_ind_matrix[peptide_list_index] = None
                del peptide_list
                del peptide_ind
                time_2 = time.time()
                print("[Info] Finishing creating single peptide ion index " + str(self.dp.LIST_MASS_INDEX_FILE[peptide_list_index].split('\\')[-1].split('.')[0]) + " using time %.4f s" % (time_2 - time_1))
            print("[info] Candidate single peptides : ", self.pep_num)

            # 输出离子索引
            if not os.path.exists(self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[1][1] + '\\'):
                os.makedirs(self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[1][1] + '\\')

            out_path = c_char_p(bytes(self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[1][1] + '\\' + FILE_NAME_ION_INDEX[1][1] + '.dc', 'utf-8'))
            bird_index.write_ion_index_bin(out_path, p_C_ion_index, self.dp.myCFG.max_mz_index_size)
            bird_index.free_ion_index(self.dp.myCFG.max_mz_index_size, p_C_ion_index)

        elif self.dp.myCFG.D15_TYPE_LINK in [2]:

            for peptide_list_index in range(len(self.dp.LIST_MASS_INDEX_FILE)):
                time_1 = time.time()
                if peptide_matrix is None:
                    peptide_list = functionPickle.load_pkl_to_data(self.dp.LIST_MASS_INDEX_FILE[peptide_list_index])
                    peptide_ind = functionPickle.load_pkl_to_data(self.dp.LIST_MASS_IND_FILE[peptide_list_index])
                else:
                    peptide_list = peptide_matrix[peptide_list_index]
                    peptide_ind = peptide_ind_matrix[peptide_list_index]
                self.create_looplink_ion_index(protein_list, peptide_list_index, peptide_list, peptide_ind, bird_function, p_mod_mass, p_aa_mass, p_C_ion_index)
                if peptide_matrix is not None:
                    peptide_matrix[peptide_list_index] = None
                    peptide_ind_matrix[peptide_list_index] = None
                del peptide_list
                del peptide_ind
                time_2 = time.time()
                print("[Info] Finishing creating loop-link peptide ion index " + str(self.dp.LIST_MASS_INDEX_FILE[peptide_list_index].split('\\')[-1].split('.')[0]) +" using time %.4f s" % (time_2 - time_1))
            print("[info] Candidate loop-link peptides : ", self.pep_num)

            # 输出离子索引
            if not os.path.exists(self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[1][2] + '\\'):
                os.makedirs(self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[1][2] + '\\')

            out_path = c_char_p(bytes(self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[1][2] + '\\' + FILE_NAME_ION_INDEX[1][2] + '.dc', 'utf-8'))
            bird_index.write_ion_index_bin(out_path, p_C_ion_index, self.dp.myCFG.max_mz_index_size)
            bird_index.free_ion_index(self.dp.myCFG.max_mz_index_size, p_C_ion_index)

        else:

            logGetError("[Error] The value of TYPE_PEPTIDE or TYPE_LINK may wrong !")
        time_2 = time.time()
        print("[Info] Using time : %.4f s."%(time_2 - time_1))
        self.__captain_out_put_ion_index()

    def __captain_C_get_p_mod_mass(self):

        # 先对这些数据结构转成C的做处理
        mod_num = len(self.dp.myMOD.fix_mod_list) + len(self.dp.myMOD.var_mod_list)
        p_mod_mass = (c_double * mod_num)()

        return p_mod_mass

    def __captain_C_get_p_aa_mass(self):

        p_aa_mass = (c_double * AA_NUM)()
        for aa_index, aa in enumerate(self.dp.myINI.DIC_AA.keys()):
            p_aa_mass[aa_index] = self.dp.myINI.DIC_AA[aa]

        return p_aa_mass

    # =================================================

    def create_single_peptide_ion_index(self, protein_list, peptide_list_index, peptide_list, peptide_ind, bird_function, p_mod_mass, p_aa_mass, p_C_ion_index):

        # 先对这些数据结构转成C的做处理
        old_ind_data = 0
        for pkl_file_ind_index, new_ind_data in enumerate(peptide_ind):
            if new_ind_data == old_ind_data:
                continue
            for one_ind_index in range(new_ind_data - old_ind_data):
                # 一条肽段
                one_pkl_peptide = peptide_list[old_ind_data + one_ind_index]
                is_target, pro_index, pep_start_pos, pep_length = tool_get_peptide_info(one_pkl_peptide[1])
                pep_end_pos = pep_start_pos + pep_length
                pep_sq = protein_list.sq[pro_index][pep_start_pos: pep_end_pos]

                p_sq = c_char_p(bytes(pep_sq, 'utf-8'))
                p_mod_site_list = (c_int * (pep_length + 2))()
                for i in range(pep_length + 2):
                    p_mod_site_list[i] = -1
                # 单肽没有link site的编码
                c_pep_code = bird_function.tool_C_get_ion_code(peptide_list_index, pkl_file_ind_index, one_ind_index, 0)

                # 产生解码修饰
                if len(one_pkl_peptide) == 6:
                    self.__captain_get_single_mod_site_list(pep_sq, pep_length, one_pkl_peptide[2], one_pkl_peptide[5], p_mod_site_list)
                else:
                    self.__captain_get_single_mod_site_list(pep_sq, pep_length, one_pkl_peptide[2], [], p_mod_site_list)

                self.__captain_put_single_peptide_in_ion_index(p_C_ion_index, c_pep_code, one_pkl_peptide[0], pep_length, p_sq, p_mod_site_list, p_mod_mass, p_aa_mass, bird_function)
            old_ind_data = new_ind_data

    def __captain_put_single_peptide_in_ion_index(self, p_C_ion_index, c_pep_code, pep_mass, pep_length, p_sq, p_mod_site_list, p_mod_mass, p_aa_mass, bird_function):
        bird_function.tool_C_put_ion_to_ion_index_singlepeptide(p_C_ion_index,
                                                                p_sq,
                                                                pep_length,
                                                                p_mod_site_list,
                                                                c_double(0.0),
                                                                c_double(self.dp.myINI.MASS_H2O),
                                                                p_mod_mass,
                                                                p_aa_mass,
                                                                self.dp.myCFG.D12_MULTI_MASS,
                                                                self.dp.myCFG.max_mz_index_size,
                                                                self.used_pep_num_add_site)
        self.list_pep_mass[self.used_pep_num_add_site] = pep_mass
        self.list_pep_data_pkl[self.used_pep_num_add_site] = c_pep_code
        self.list_ion_num[self.used_pep_num_add_site] = (pep_length - 1) * 2
        self.used_pep_num_add_site += 1

    def __captain_get_single_mod_site_list(self, pep_sq, pep_length, fix_mod_code, var_mod_code_list, p_mod_site_list):
        # 这里为了节约时间先走可变修饰再走固定修饰，因为固定修饰需要遍历全部位点,记录可以加交联位点的位置，这样就不用遍历整个肽段了
        for var_mod_code in var_mod_code_list:
            for i in range(5):
                one_var_mod_code = var_mod_code >> (12 * i) & 0xfff
                if one_var_mod_code == 0:
                    break
                var_mod_index = one_var_mod_code & 0x3f
                var_mod_site = one_var_mod_code >> 6 & 0x3f
                p_mod_site_list[var_mod_site] = var_mod_index
        for site in range(pep_length + 2):
            # site表示在修饰list中的位置，-1为第几个氨基酸，因为0是N端，索引氨基酸就是site-2
            site_data = fix_mod_code & (1 << site)
            # 这里优先处理修饰列表的值
            if site_data != 0:
                p_mod_site_list[site] = site_data
                if site == 0:
                    if pep_sq[0] in self.dp.myMOD.fix_mod_dic_pro_N.keys():
                        p_mod_site_list[site] = self.dp.myMOD.fix_mod_dic_pro_N[pep_sq[0]]
                    else:
                        if pep_sq[0] in self.dp.myMOD.fix_mod_dic_pep_N.keys():
                            p_mod_site_list[site] = self.dp.myMOD.fix_mod_dic_pep_N[pep_sq[0]]
                elif site == pep_length + 1:
                    if pep_sq[-1] in self.dp.myMOD.fix_mod_dic_pro_C.keys():
                        p_mod_site_list[site] = self.dp.myMOD.fix_mod_dic_pro_C[pep_sq[-1]]
                    else:
                        if pep_sq[-1] in self.dp.myMOD.fix_mod_dic_pep_C.keys():
                            p_mod_site_list[site] = self.dp.myMOD.fix_mod_dic_pep_C[pep_sq[-1]]
                else:
                    p_mod_site_list[site] = self.dp.myMOD.fix_mod_dic[pep_sq[site - 1]]

    # =================================================

    def create_looplink_ion_index(self, protein_list, peptide_list_index, peptide_list, peptide_ind, bird_function, p_mod_mass, p_aa_mass, p_C_ion_index):

        # 先对这些数据结构转成C的做处理
        old_ind_data = 0
        for pkl_file_ind_index, new_ind_data in enumerate(peptide_ind):
            if new_ind_data == old_ind_data:
                continue
            for one_ind_index in range(new_ind_data - old_ind_data):
                # 一条肽段
                one_pkl_peptide = peptide_list[old_ind_data + one_ind_index]
                is_target, pro_index, pep_start_pos, pep_length = tool_get_peptide_info(one_pkl_peptide[1])
                pep_end_pos = pep_start_pos + pep_length
                pep_sq = protein_list.sq[pro_index][pep_start_pos: pep_end_pos]

                p_sq = c_char_p(bytes(pep_sq, 'utf-8'))
                p_mod_site_list = (c_int * (pep_length + 2))()
                for i in range(pep_length + 2):
                    p_mod_site_list[i] = -1

                pro_length = len(protein_list.sq[pro_index])
                c_pep_code = bird_function.tool_C_get_ion_code(peptide_list_index, pkl_file_ind_index, one_ind_index, 0)

                candidate_site = []
                if len(one_pkl_peptide) == 6:
                    self.__captain_get_mod_loop_link_site_list(pep_sq, pep_length, pro_length, pep_start_pos, one_pkl_peptide[2], one_pkl_peptide[5], p_mod_site_list, candidate_site)
                else:
                    self.__captain_get_mod_loop_link_site_list(pep_sq, pep_length, pro_length, pep_start_pos, one_pkl_peptide[2], [], p_mod_site_list, candidate_site)

                if len(candidate_site) < 1:
                    continue
                else:
                    candidate_loop_link_site_combine = []
                    self.__captain_get_loop_link_site_combine(candidate_site, candidate_loop_link_site_combine)
                for one_loop_link in candidate_loop_link_site_combine:
                    if one_loop_link[0][0] == 1 and one_loop_link[1][0] == pep_length:
                        continue
                    self.__captain_put_looplink_peptide_in_ion_index(p_C_ion_index, c_pep_code, one_pkl_peptide[0],
                                                                     pep_length, p_sq, p_mod_site_list,
                                                                     p_mod_mass, p_aa_mass, bird_function,
                                                                     one_loop_link)

            old_ind_data = new_ind_data

    def __captain_get_loop_link_site(self, aa, site, candidate_site):

        if self.dp.myLINK.linker_1_data.same_site[aa] != 0:
            # 两端都支持的
            candidate_site.append((site, 0))
            exist = 1

        elif self.dp.myLINK.linker_1_data.alpha_site[aa] != 0:
            # 一端都支持的
            candidate_site.append((site, 1))
            exist = 1

        elif self.dp.myLINK.linker_1_data.beta_site[aa] != 0:
            # 一端都支持的
            candidate_site.append((site, -1))
            exist = 1
        else:
            exist = 0

        return exist

    def __captain_get_loop_link_site_combine(self, candidate_site, out_loop_link_site_combine):

        for link_site_1 in candidate_site:
            for link_site_2 in candidate_site:
                if link_site_1[0] >= link_site_2[0]:
                    continue
                else:
                    if link_site_1[1] + link_site_2[1] != 0:
                        continue
                    else:
                        out_loop_link_site_combine.append((link_site_1, link_site_2))

    def __captain_put_looplink_peptide_in_ion_index(self, p_C_ion_index, c_pep_code, pep_mass, pep_length, p_sq, p_mod_site_list, p_mod_mass, p_aa_mass, bird_function, link_site):
        # 在检索looplink的时候减掉linker的质量了，所以这里直接用肽段质量即可
        link_site_1 = link_site[0][0]
        link_site_2 = link_site[1][0]

        bird_function.tool_C_put_ion_to_ion_index_looplink(p_C_ion_index, p_sq, pep_length,
                                                  p_mod_site_list, c_double(0.0), c_double(self.dp.myINI.MASS_H2O),
                                                  p_mod_mass, p_aa_mass, link_site_1, link_site_2,
                                                  c_double(self.dp.myLINK.linker_1_data.loop_mass),
                                                  self.dp.myCFG.D12_MULTI_MASS,
                                                  self.dp.myCFG.max_mz_index_size,
                                                  self.used_pep_num_add_site)
        self.list_pep_mass[self.used_pep_num_add_site] = pep_mass
        self.list_pep_data_pkl[self.used_pep_num_add_site] = c_pep_code
        self.list_ion_num[self.used_pep_num_add_site] = (pep_length - 1 - abs(link_site_1 - link_site_2)) * 2
        self.list_link_site[self.used_pep_num_add_site] = link_site
        self.used_pep_num_add_site += 1

        self.pep_num += 1
        if self.used_pep_num_add_site > self.dp.myCFG.B9_MAX_PEPNUM:
            logGetError("The MAX_PEP_NUM may small, but your computer can finished this work. please change the MAX_PEP_NUM !")

    def __captain_get_mod_loop_link_site_list(self, pep_sq, pep_length, pro_length, pep_start_pos,
                                              fix_mod_code, var_mod_code_list,
                                              p_mod_site_list, link_site_list):
        # 这里为了节约时间先走可变修饰再走固定修饰，因为固定修饰需要遍历全部位点,记录可以加交联位点的位置，这样就不用遍历整个肽段了
        pep_end_pos = pep_start_pos + pep_length
        for var_mod_code in var_mod_code_list:
            for i in range(5):
                one_var_mod_code = var_mod_code >> (12 * i) & 0xfff
                if one_var_mod_code == 0:
                    break
                var_mod_index = one_var_mod_code & 0x3f
                var_mod_site = one_var_mod_code >> 6 & 0x3f
                p_mod_site_list[var_mod_site] = var_mod_index
        for site in range(pep_length + 2):
            # site表示在修饰list中的位置，-1为第几个氨基酸，因为0是N端，索引氨基酸就是site-2
            site_data = fix_mod_code & (1 << site)
            # 这里优先处理修饰列表的值
            if site_data != 0:
                p_mod_site_list[site] = site_data
                if site == 0:
                    if pep_sq[0] in self.dp.myMOD.fix_mod_dic_pro_N.keys():
                        p_mod_site_list[site] = self.dp.myMOD.fix_mod_dic_pro_N[pep_sq[0]]
                    else:
                        if pep_sq[0] in self.dp.myMOD.fix_mod_dic_pep_N.keys():
                            p_mod_site_list[site] = self.dp.myMOD.fix_mod_dic_pep_N[pep_sq[0]]
                elif site == pep_length + 1:
                    if pep_sq[-1] in self.dp.myMOD.fix_mod_dic_pro_C.keys():
                        p_mod_site_list[site] = self.dp.myMOD.fix_mod_dic_pro_C[pep_sq[-1]]
                    else:
                        if pep_sq[-1] in self.dp.myMOD.fix_mod_dic_pep_C.keys():
                            p_mod_site_list[site] = self.dp.myMOD.fix_mod_dic_pep_C[pep_sq[-1]]
                else:
                    p_mod_site_list[site] = self.dp.myMOD.fix_mod_dic[pep_sq[site - 1]]
            # 经过上面判断后再确定交联位点
            if p_mod_site_list[site] == -1:
                exist_link_site = False
                if site == 0:
                    if pep_start_pos == 0 and self.link_1_pro_N:
                        # 这个是蛋白N端交联
                        self.__soldier_check_loop_link_site(-2, 1, link_site_list)
                        exist_link_site = True
                    if not exist_link_site and self.link_1_pep_N:
                        self.__soldier_check_loop_link_site(-1, 1, link_site_list)
                elif site == pep_length + 1:
                    if (pro_length == pep_end_pos) and self.link_1_pro_C:
                        # 蛋白C端
                        self.__soldier_check_loop_link_site(2, pep_length, link_site_list)
                        exist_link_site = True
                    if not exist_link_site and self.link_1_pep_C:
                        self.__soldier_check_loop_link_site(1, pep_length, link_site_list)
                        exist_link_site = True
                else:
                    aa = pep_sq[site - 1]
                    if self.dp.myLINK.linker_1_data.dic_site[aa]:
                        # 肽段C端需要再处理一下，只有蛋白C端正好是交联位点氨基酸才可以作为交联肽段
                        if site != pep_length:
                            self.__soldier_check_loop_link_site(aa, site, link_site_list)
                        else:
                            if self.dp.myENZ.enzyme_C and pro_length == pep_end_pos:
                                self.__soldier_check_loop_link_site(aa, site, link_site_list)
                                exist_link_site = True

    def __soldier_check_loop_link_site(self, aa, site, candidate_site):
        if self.dp.myLINK.linker_1_data.same_site[aa] != 0:
            # 两端都支持的
            candidate_site.append((site, 0))

        elif self.dp.myLINK.linker_1_data.alpha_site[aa] != 0:
            # 一端都支持的
            candidate_site.append((site, 1))

        elif self.dp.myLINK.linker_1_data.beta_site[aa] != 0:
            # 一端都支持的
            candidate_site.append((site, -1))
        else:
            pass

    # =================================================

    def __captain_out_put_ion_index(self):

        functionPickle = CFunctionPickle()

        path = self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[self.dp.myCFG.D14_TYPE_PEPTIDE][self.dp.myCFG.D15_TYPE_LINK] + '\\' + \
               FILE_NAME_ION_NUM[self.dp.myCFG.D14_TYPE_PEPTIDE][self.dp.myCFG.D15_TYPE_LINK] + '.pkl'
        functionPickle.write_data_to_pkl(tuple(self.list_ion_num[: self.used_pep_num_add_site]), path)

        path = self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[self.dp.myCFG.D14_TYPE_PEPTIDE][self.dp.myCFG.D15_TYPE_LINK] + '\\' + \
               FILE_NAME_PEP_MASS[self.dp.myCFG.D14_TYPE_PEPTIDE][self.dp.myCFG.D15_TYPE_LINK] + '.pkl'
        functionPickle.write_data_to_pkl(tuple(self.list_pep_mass[: self.used_pep_num_add_site]), path)

        path = self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[self.dp.myCFG.D14_TYPE_PEPTIDE][self.dp.myCFG.D15_TYPE_LINK] + '\\' + \
               FILE_NAME_PEP_PKL[self.dp.myCFG.D14_TYPE_PEPTIDE][self.dp.myCFG.D15_TYPE_LINK] + '.pkl'
        functionPickle.write_data_to_pkl(tuple(self.list_pep_data_pkl[: self.used_pep_num_add_site]), path)

        if self.link_type not in FILE_NAME_LINK_SITE[self.dp.myCFG.D15_TYPE_LINK].keys():
            pass
        else:
            path = self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[self.dp.myCFG.D14_TYPE_PEPTIDE][self.dp.myCFG.D15_TYPE_LINK] + '\\' + \
                   FILE_NAME_LINK_SITE[self.dp.myCFG.D14_TYPE_PEPTIDE][self.dp.myCFG.D15_TYPE_LINK] + '.pkl'
            functionPickle.write_data_to_pkl(tuple(self.list_link_site[: self.used_pep_num_add_site]), path)

class CFunctionCreateTwoPeptideIonIndex:

    def __init__(self, inputDP):

        self.dp = inputDP

        self.pep_num = 0
        self.used_pep_num_add_site = 0
        self.list_pep_data_pkl = [None] * self.dp.myCFG.B9_MAX_PEPNUM  # 这个是记录满足交联条件的肽段索引，index->肽段编码
        self.list_ion_num = [None] * self.dp.myCFG.B9_MAX_PEPNUM  # 这个是记录肽段的碎片离子数目，index->离子数
        self.list_pep_mass = [None] * self.dp.myCFG.B9_MAX_PEPNUM  # 这个是记录肽段质量的，不包括交联剂质量
        self.list_link_site = [None] * self.dp.myCFG.B9_MAX_PEPNUM

        self.link_1_pro_N = self.dp.myLINK.linker_1_data.same_site[-2] or self.dp.myLINK.linker_1_data.alpha_site[-2] or self.dp.myLINK.linker_1_data.beta_site[-2]
        self.link_1_pep_N = self.dp.myLINK.linker_1_data.same_site[-1] or self.dp.myLINK.linker_1_data.alpha_site[-1] or self.dp.myLINK.linker_1_data.beta_site[-1]
        self.link_1_pep_C = self.dp.myLINK.linker_1_data.same_site[1] or self.dp.myLINK.linker_1_data.alpha_site[1] or self.dp.myLINK.linker_1_data.beta_site[1]
        self.link_1_pro_C = self.dp.myLINK.linker_1_data.same_site[2] or self.dp.myLINK.linker_1_data.alpha_site[2] or self.dp.myLINK.linker_1_data.beta_site[2]

    def import_DLL(self, bird_index, bird_function):

        # 导入动态链接库

        bird_index.malloc_ion_index.argtype = (c_int)
        bird_index.malloc_ion_index.restype = (POINTER(ION_INDEX_MZ))

        bird_index.free_ion_index.argtype = (c_int, POINTER(ION_INDEX_MZ))
        bird_index.free_ion_index.restype = (c_int)

        bird_index.put_peptide_to_ion_index.argtype = (ION_INDEX_MZ, c_double, c_int, c_int, c_long)
        bird_index.put_peptide_to_ion_index.restype = (c_int)

        bird_index.write_ion_index_bin.argtype = (c_char_p, POINTER(ION_INDEX_MZ), c_int)
        bird_index.write_ion_index_bin.restype = (c_int)

        # 导入动态链接库

        bird_function.tool_C_get_ion_code.argtype = (c_longlong, c_longlong, c_longlong, c_longlong)
        bird_function.tool_C_get_ion_code.restype = (c_longlong)

        bird_function.tool_C_new_get_ion_data.argtype = (c_longlong)
        bird_function.tool_C_new_get_ion_data.restype = POINTER(c_long)

        bird_function.tool_C_put_ion_to_ion_index_singlepeptide.argtype = (POINTER(ION_INDEX_MZ), c_char_p, c_int, POINTER(c_int), c_double, c_double, POINTER(c_double), POINTER(c_double), c_int, c_int, c_long)
        bird_function.tool_C_put_ion_to_ion_index_singlepeptide.restype = (c_int)

        bird_function.tool_C_put_ion_to_ion_index_crosslink.argtype = (POINTER(ION_INDEX_MZ), c_char_p, c_int, POINTER(c_int), c_double, c_double, POINTER(c_double), POINTER(c_double), c_int, c_int, c_int, c_long)
        bird_function.tool_C_put_ion_to_ion_index_crosslink.restype = (c_int)

        bird_function.tool_C_put_ion_to_ion_index_looplink.argtype = (POINTER(ION_INDEX_MZ), c_char_p, c_int, POINTER(c_int), c_double, c_double, POINTER(c_double), POINTER(c_double), c_int, c_int, c_double, c_int, c_int, c_long)
        bird_function.tool_C_put_ion_to_ion_index_looplink.restype = (c_int)

        bird_function.tool_C_put_ion_to_ion_index_crosscross.argtype = (POINTER(ION_INDEX_MZ), c_char_p, c_int, POINTER(c_int), c_double, c_double, POINTER(c_double), POINTER(c_double), c_int, c_int, c_int, c_int, c_long)
        bird_function.tool_C_put_ion_to_ion_index_crosscross.restype = (c_int)

        bird_function.tool_C_put_ion_to_ion_index_crossloop.argtype = (POINTER(ION_INDEX_MZ), c_char_p, c_int, POINTER(c_int), c_double, c_double, POINTER(c_double), POINTER(c_double), c_int, c_int, c_double, c_int, c_int, c_int, c_long)
        bird_function.tool_C_put_ion_to_ion_index_crossloop.restype = (c_int)

    def create_ion_index(self, protein_list, peptide_matrix=None, peptide_ind_matrix=None):

        # 导入动态链接库
        bird_index = CDLL(os.getcwd() + "\\" + DLL_INDEX)
        bird_function = CDLL(os.getcwd() + "\\" + DLL_FUNCTION)

        self.import_DLL(bird_index, bird_function)

        p_mod_mass = self.__captain_C_get_p_mod_mass()
        p_aa_mass = self.__captain_C_get_p_aa_mass()
        p_C_ion_index = bird_index.malloc_ion_index(self.dp.myCFG.max_mz_index_size)

        functionPickle = CFunctionPickle()
        time_start = time.time()
        if self.dp.myCFG.D15_TYPE_LINK in [1]:
            # cross
            # 先判断存索引的文件夹是否存在
            if not os.path.exists(self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[2][1] + '\\'):
                os.makedirs(self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[2][1] + '\\')
            # 每个质量区间遍历
            for peptide_list_index in range(len(self.dp.LIST_MASS_INDEX_FILE)):
                time_1 = time.time()
                if peptide_matrix is None:
                    peptide_list = functionPickle.load_pkl_to_data(self.dp.LIST_MASS_INDEX_FILE[peptide_list_index])
                    peptide_ind = functionPickle.load_pkl_to_data(self.dp.LIST_MASS_IND_FILE[peptide_list_index])
                else:
                    peptide_list = peptide_matrix[peptide_list_index]
                    peptide_ind = peptide_ind_matrix[peptide_list_index]
                last_num = self.pep_num
                self.create_only_cross_ion_index(protein_list, peptide_list_index, peptide_list, peptide_ind, bird_function, p_mod_mass, p_aa_mass, p_C_ion_index)

                if peptide_matrix is not None:
                    peptide_matrix[peptide_list_index] = None
                    peptide_ind_matrix[peptide_list_index] = None
                del peptide_list
                del peptide_ind
                time_2 = time.time()
                print("[Info] Finishing creating crosslink peptide ion index " +
                      str(self.dp.LIST_MASS_INDEX_FILE[peptide_list_index].split('\\')[-1].split('.')[0]) +
                      " using time %.4f s" % (time_2 - time_1) + " total peptides : %d" % (self.pep_num - last_num))
            print("[info] Candidate crosslink peptides : ", self.pep_num)
            # 输出离子索引
            out_path = self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[2][1] + '\\' + FILE_NAME_ION_INDEX[2][1] + '.dc'
            p_out_path = c_char_p(bytes(out_path, 'utf-8'))
            bird_index.write_ion_index_bin(p_out_path, p_C_ion_index, self.dp.myCFG.max_mz_index_size)
            bird_index.free_ion_index(self.dp.myCFG.max_mz_index_size, p_C_ion_index)
            # 其他的索引文件
            self.__captain_out_put_ion_index_data()

        elif self.dp.myCFG.D15_TYPE_LINK in [2]:
            pass
            # cross+cross形式
            # self.create_cross_cross_ion_index(protein_list, peptide_list, bird_index, bird_function)
        elif self.dp.myCFG.D15_TYPE_LINK in [3]:
            pass
            # cross+loop形式
            # self.create_cross_loop_ion_index(protein_list, peptide_list, bird_index, bird_function)
        else:
            logGetError("[Error] The value of TYPE_PEPTIDE or TYPE_LINK may wrong !")

        time_end = time.time()
        print("[Info] Using time : %.4f s."%(time_end - time_start))

    def create_only_cross_ion_index(self, protein_list, peptide_list_index, peptide_list, peptide_ind,
                                    bird_function, p_mod_mass, p_aa_mass, p_C_ion_index):
        # 这个是C版本的建索引
        old_ind_data = 0
        for pkl_file_ind_index, new_ind_data in enumerate(peptide_ind):
            if new_ind_data == old_ind_data:
                continue
            for one_ind_index in range(new_ind_data - old_ind_data):
                # 一条肽段
                one_pkl_peptide = peptide_list[old_ind_data + one_ind_index]
                is_target, pro_index, pep_start_pos, pep_length = tool_get_peptide_info(one_pkl_peptide[1])
                pro_length = len(protein_list.sq[pro_index])
                pep_end_pos = pep_start_pos + pep_length
                pep_sq = protein_list.sq[pro_index][pep_start_pos: pep_end_pos]

                p_sq = c_char_p(bytes(pep_sq, 'utf-8'))
                p_mod_site_list = (c_int * (pep_length + 2))()
                for i in range(pep_length + 2):
                    p_mod_site_list[i] = -1
                link_site_list = []
                # 这个里面-4 -3 -2 -1分别对应蛋白N端、肽段N端、肽段C端、蛋白C端，其余数字为氨基酸索引
                if len(one_pkl_peptide) == 6:
                    self.__captain_get_mod_link_site_list(pep_sq, pep_length, pro_length, pep_start_pos,
                                                          one_pkl_peptide[2], one_pkl_peptide[5],
                                                          p_mod_site_list, link_site_list)
                else:
                    self.__captain_get_mod_link_site_list(pep_sq, pep_length, pro_length, pep_start_pos,
                                                          one_pkl_peptide[2], [],
                                                          p_mod_site_list, link_site_list)
                link_site_num = 0
                # 不带上交联位点的编码
                c_pep_without_site_code = bird_function.tool_C_get_ion_code(peptide_list_index, pkl_file_ind_index, one_ind_index, link_site_num)
                is_link_peptide = False
                for link_site in link_site_list:
                    if link_site == -4:
                        exist = self.__captain_put_crosslink_peptide_in_ion_index(p_C_ion_index,
                                                                                  c_pep_without_site_code,
                                                                                  one_pkl_peptide[0], pep_length,
                                                                                  p_sq, p_mod_site_list, p_mod_mass,
                                                                                  p_aa_mass,
                                                                                  -2, 1,
                                                                                  bird_function,
                                                                                  )
                    elif link_site == -3:
                        exist = self.__captain_put_crosslink_peptide_in_ion_index(p_C_ion_index,
                                                                                  c_pep_without_site_code,
                                                                                  one_pkl_peptide[0], pep_length,
                                                                                  p_sq, p_mod_site_list, p_mod_mass,
                                                                                  p_aa_mass,
                                                                                  -1, 1,
                                                                                  bird_function,
                                                                                  )
                    elif link_site == -2:
                        exist = self.__captain_put_crosslink_peptide_in_ion_index(p_C_ion_index,
                                                                                  c_pep_without_site_code,
                                                                                  one_pkl_peptide[0], pep_length,
                                                                                  p_sq, p_mod_site_list, p_mod_mass,
                                                                                  p_aa_mass,
                                                                                  2, pep_length,
                                                                                  bird_function,
                                                                                  )
                    elif link_site == -1:
                        exist = self.__captain_put_crosslink_peptide_in_ion_index(p_C_ion_index,
                                                                                  c_pep_without_site_code,
                                                                                  one_pkl_peptide[0], pep_length,
                                                                                  p_sq, p_mod_site_list, p_mod_mass,
                                                                                  p_aa_mass,
                                                                                  1, pep_length,
                                                                                  bird_function,
                                                                                  )
                    else:
                        aa = pep_sq[link_site]
                        exist = self.__captain_put_crosslink_peptide_in_ion_index(p_C_ion_index,
                                                                                  c_pep_without_site_code,
                                                                                  one_pkl_peptide[0], pep_length,
                                                                                  p_sq, p_mod_site_list, p_mod_mass,
                                                                                  p_aa_mass,
                                                                                  aa, link_site + 1,
                                                                                  bird_function
                                                                                  )
                    if exist:
                        is_link_peptide = True
                if is_link_peptide:
                    self.pep_num += 1

            old_ind_data = new_ind_data

    def __captain_get_mod_link_site_list(self, pep_sq, pep_length, pro_length, pep_start_pos, fix_mod_code, var_mod_code_list, p_mod_site_list, link_site_list):
        # 这里为了节约时间先走可变修饰再走固定修饰，因为固定修饰需要遍历全部位点,记录可以加交联位点的位置，这样就不用遍历整个肽段了
        pep_end_pos = pep_start_pos + pep_length
        for var_mod_code in var_mod_code_list:
            for i in range(5):
                one_var_mod_code = var_mod_code >> (12 * i) & 0xfff
                if one_var_mod_code == 0:
                    break
                var_mod_index = one_var_mod_code & 0x3f
                var_mod_site = one_var_mod_code >> 6 & 0x3f
                p_mod_site_list[var_mod_site] = var_mod_index
        for site in range(pep_length + 2):
            # site表示在修饰list中的位置，-1为第几个氨基酸，因为0是N端，索引氨基酸就是site-2
            site_data = fix_mod_code & (1 << site)
            # 这里优先处理修饰列表的值
            if site_data != 0:
                p_mod_site_list[site] = site_data
                if site == 0:
                    if pep_sq[0] in self.dp.myMOD.fix_mod_dic_pro_N.keys():
                        p_mod_site_list[site] = self.dp.myMOD.fix_mod_dic_pro_N[pep_sq[0]]
                    else:
                        if pep_sq[0] in self.dp.myMOD.fix_mod_dic_pep_N.keys():
                            p_mod_site_list[site] = self.dp.myMOD.fix_mod_dic_pep_N[pep_sq[0]]
                elif site == pep_length + 1:
                    if pep_sq[-1] in self.dp.myMOD.fix_mod_dic_pro_C.keys():
                        p_mod_site_list[site] = self.dp.myMOD.fix_mod_dic_pro_C[pep_sq[-1]]
                    else:
                        if pep_sq[-1] in self.dp.myMOD.fix_mod_dic_pep_C.keys():
                            p_mod_site_list[site] = self.dp.myMOD.fix_mod_dic_pep_C[pep_sq[-1]]
                else:
                    p_mod_site_list[site] = self.dp.myMOD.fix_mod_dic[pep_sq[site - 1]]
            # 经过上面判断后再确定交联位点
            if p_mod_site_list[site] == -1:
                exist_link_site = False
                if site == 0:
                    if pep_start_pos == 0 and self.link_1_pro_N:
                        # 这个是蛋白N端交联
                        link_site_list.append(-4)
                        exist_link_site = True
                    if not exist_link_site and self.link_1_pep_N:
                        link_site_list.append(-3)
                elif site == pep_length + 1:
                    if (pro_length == pep_end_pos) and self.link_1_pro_C:
                        # 蛋白C端
                        link_site_list.append(-2)
                        exist_link_site = True
                    if not exist_link_site and self.link_1_pep_C:
                        link_site_list.append(-1)
                        exist_link_site = True
                else:
                    aa = pep_sq[site - 1]
                    if self.dp.myLINK.linker_1_data.dic_site[aa]:
                        # 肽段C端需要再处理一下，只有蛋白C端正好是交联位点氨基酸才可以作为交联肽段
                        if site != pep_length:
                            link_site_list.append(site - 1)
                        else:
                            if self.dp.myENZ.enzyme_C and pro_length == pep_end_pos:
                                link_site_list.append(site - 1)
                                exist_link_site = True

    def __captain_C_get_p_mod_mass(self):

        # 先对这些数据结构转成C的做处理
        mod_num = len(self.dp.myMOD.fix_mod_list) + len(self.dp.myMOD.var_mod_list)
        p_mod_mass = (c_double * mod_num)()

        return p_mod_mass

    def __captain_C_get_p_aa_mass(self):

        p_aa_mass = (c_double * AA_NUM)()
        for aa_index, aa in enumerate(self.dp.myINI.DIC_AA.keys()):
            p_aa_mass[aa_index] = self.dp.myINI.DIC_AA[aa]

        return p_aa_mass

    def __captain_put_crosslink_peptide_in_ion_index(self, p_C_ion_index, c_pep_without_site_code, pep_mass, pep_length, p_sq, p_mod_site_list, p_mod_mass, p_aa_mass, aa, link_site, bird_function):
        exist = False

        if self.dp.myLINK.linker_1_data.same_site[aa] != 0:
            # 两端都支持的
            bird_function.tool_C_put_ion_to_ion_index_crosslink(p_C_ion_index, p_sq, pep_length,
                                                      p_mod_site_list, c_double(0.0), c_double(self.dp.myINI.MASS_H2O),
                                                      p_mod_mass, p_aa_mass, link_site,
                                                      self.dp.myCFG.D12_MULTI_MASS,
                                                      self.dp.myCFG.max_mz_index_size,
                                                      self.used_pep_num_add_site)
            site_flag = 0
            exist = True

        elif self.dp.myLINK.linker_1_data.alpha_site[aa] != 0:
            # 一端都支持的
            bird_function.tool_C_put_ion_to_ion_index_crosslink(p_C_ion_index, p_sq, pep_length,
                                                      p_mod_site_list, c_double(0.0), c_double(self.dp.myINI.MASS_H2O),
                                                      p_mod_mass, p_aa_mass, link_site,
                                                      self.dp.myCFG.D12_MULTI_MASS,
                                                      self.dp.myCFG.max_mz_index_size,
                                                      self.used_pep_num_add_site)
            site_flag = -1
            exist = True

        elif self.dp.myLINK.linker_1_data.beta_site[aa] != 0:
            # 一端都支持的
            bird_function.tool_C_put_ion_to_ion_index_crosslink(p_C_ion_index, p_sq, pep_length,
                                                      p_mod_site_list, c_double(0.0), c_double(self.dp.myINI.MASS_H2O),
                                                      p_mod_mass, p_aa_mass, link_site,
                                                      self.dp.myCFG.D12_MULTI_MASS,
                                                      self.dp.myCFG.max_mz_index_size,
                                                      self.used_pep_num_add_site)
            site_flag = 1
            exist = True

        else:
            exist = 0
            site_flag = None
        
        if exist == True:
            self.list_pep_mass[self.used_pep_num_add_site] = pep_mass
            self.list_pep_data_pkl[self.used_pep_num_add_site] = c_pep_without_site_code
            self.list_ion_num[self.used_pep_num_add_site] = pep_length
            self.list_link_site[self.used_pep_num_add_site] = (link_site, site_flag)
            self.used_pep_num_add_site += 1

        return exist

    # =================================================
    def __captain_out_put_ion_index_data(self):

        functionPickle = CFunctionPickle()

        path = self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[self.dp.myCFG.D14_TYPE_PEPTIDE][self.dp.myCFG.D15_TYPE_LINK] + '\\' + FILE_NAME_ION_NUM[self.dp.myCFG.D14_TYPE_PEPTIDE][self.dp.myCFG.D15_TYPE_LINK] + '.pkl'
        functionPickle.write_data_to_pkl(tuple(self.list_ion_num[: self.used_pep_num_add_site]), path)

        del self.list_ion_num

        path = self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[self.dp.myCFG.D14_TYPE_PEPTIDE][self.dp.myCFG.D15_TYPE_LINK] + '\\' + FILE_NAME_PEP_MASS[self.dp.myCFG.D14_TYPE_PEPTIDE][self.dp.myCFG.D15_TYPE_LINK] + '.pkl'
        functionPickle.write_data_to_pkl(tuple(self.list_pep_mass[: self.used_pep_num_add_site]), path)

        del self.list_pep_mass

        path = self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[self.dp.myCFG.D14_TYPE_PEPTIDE][self.dp.myCFG.D15_TYPE_LINK] + '\\' + FILE_NAME_PEP_PKL[self.dp.myCFG.D14_TYPE_PEPTIDE][self.dp.myCFG.D15_TYPE_LINK] + '.pkl'
        functionPickle.write_data_to_pkl(tuple(self.list_pep_data_pkl[: self.used_pep_num_add_site]), path)

        del self.list_pep_data_pkl

        if self.dp.myCFG.D15_TYPE_LINK not in FILE_NAME_LINK_SITE[self.dp.myCFG.D14_TYPE_PEPTIDE].keys():
            pass
        else:
            path = self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[self.dp.myCFG.D14_TYPE_PEPTIDE][self.dp.myCFG.D15_TYPE_LINK] + '\\' + FILE_NAME_LINK_SITE[self.dp.myCFG.D14_TYPE_PEPTIDE][self.dp.myCFG.D15_TYPE_LINK] + '.pkl'
            functionPickle.write_data_to_pkl(tuple(self.list_link_site[: self.used_pep_num_add_site]), path)

            del self.list_link_site