# -*- mode: python ; coding: utf-8 -*-..
from MSFunction import CFunctionEnzyme, CFunctionPickle
from MSOperator import op_restore_fix_mod_site_list, op_copy_CPeptide
from MSOperator import op_get_fix_mod_code, op_get_fix_mod_site_list
from MSData import CPeptidePKL
from MSTool import tool_get_N_combine
from MSTool_code import tool_get_peptide_info_code, tool_get_peptide_info, tool_update_var_mod_code_by_site_data
from MSTool_code import tool_get_fix_mod_code, tool_get_fix_mod_site_list, tool_get_var_mod_code, tool_get_var_mod_data
from MSSysterm import AA_NUM, DLL_FUNCTION, TMP_FOLDER
from MSLogging import logGetError

from ctypes import *

import os
import time
import copy

class CFunctionAddSite():

    def __init__(self, inputDP):

        self.dp = inputDP
        self.pep_num = 0

        self.p_DIC_AA = None
        self.p_G_matrix = None

        self.record_peptide_sq = []

    def import_DLL_function(self):

        # 导入动态链接库

        bird = CDLL(os.getcwd() + '\\' + DLL_FUNCTION)

        self.p_DIC_AA = (c_double * len(self.dp.myINI.DIC_AA))()

        for aa_index, aa in enumerate(self.dp.myINI.DIC_AA.keys()):
            self.p_DIC_AA[aa_index] = self.dp.myINI.DIC_AA[aa]

        G_line_num = len(self.dp.G_matrix)
        G_row_num = len(self.dp.G_matrix[0])
        self.p_G_matrix = (c_double * (G_line_num * G_row_num))()

        for line_index, G_line in enumerate(self.dp.G_matrix):
            for row_index, G_value in enumerate(G_line):
                self.p_G_matrix[line_index * G_row_num + row_index] = G_value

        bird.calculate_peptide_mass.argtype = (c_char_p, c_int, c_double, POINTER(c_double))
        bird.calculate_peptide_mass.restype = (c_double)

        bird.calculate_gdm.argtype = (c_char_p, c_int, POINTER(c_int), c_int, POINTER(c_double), c_int, c_int)
        bird.calculate_gdm.restype = (c_double)

        return bird

    def create_fix_peptide_list(self, protein_list, peptide_set):

        # 导入动态链接库
        bird = self.import_DLL_function()
        time_1 = time.time()
        functionEnzyme = CFunctionEnzyme(self.dp)
        for pro_index, protein_sq in enumerate(protein_list.sq):  # 遍历蛋白质list,一个个蛋白来，一行是名称ac，描述de，序列sq
            # print("[Info] protein name: " + protein_list.ac[pro_index])
            if protein_list.ac[pro_index].startswith("REV_"):
                is_target = 0
            else:
                is_target = 1
            split_sites = []
            functionEnzyme.find_split_site(protein_sq, split_sites)
            # specific cut
            # For example: protein sequence is "ACDEFKGGGRPSQED"
            # peptide_seq_list, peptide_pos_list: "ACDEFK", "GGGR", "PSQED"; [0,6],[6,10],[10,15]
            # print("protein_num : " + str(len(protein_list.sq)) + " now_num : " + str(pro_index))
            for i in range(len(split_sites) - 1):
                start_site = split_sites[i] + 1  # 开始加1是因为切片获取的时候，C端切就不包括位点那个位置
                for j in range(1, self.dp.myCFG.B3_NUMBER_MAX_MISS_CLV + 2):
                    if i + j >= len(split_sites):
                        break
                    end_site = split_sites[i + j] + 1  # 结尾加1是因为切片获取右边位置需要对应数只能取到左边一个
                    if not self.check_length(start_site, end_site):
                        continue
                    pep_sq = protein_sq[start_site: end_site]  # 满足条件的序列
                    pep_length = end_site - start_site
                    self.__captain_generate_peptide_dataset_list(peptide_set, pep_sq, is_target, pro_index, start_site, end_site, pep_length, bird)
                    # split the protein when the first amino acid is "M"，还要补一个M开头剔除的序列
                    if start_site == 0 and pep_sq[0] == "M":
                        if not self.check_length(1, end_site):
                            continue
                        pep_sq = protein_sq[1:end_site]
                        pep_length = end_site - 1
                        self.__captain_generate_peptide_dataset_list(peptide_set, pep_sq, is_target, pro_index, 1, end_site, pep_length, bird)
            if pro_index % 5000 == 0:
                print("[Info] Enzyme protein : {} ".format("%.4f %% "%(pro_index * 100 / len(protein_list.sq)) + str(pro_index) + " / " + str(len(protein_list.sq))))
        time_2 = time.time()
        print("[Info] Enzyme protein : 100%")
        print("[Info] Enzyme protein get peptides : ", self.pep_num)
        print("[Info] Enzyme protein time : %.4f s ."%(time_2 - time_1))

    def __captain_generate_peptide_dataset_list(self, peptide_set, pep_sq, is_target, pro_index, pep_pos_start, pep_pos_end, pep_length, bird):

        fix_mod_code = 0  # 这个是用来记录固定修饰的
        p_mod_site = (c_int * (pep_length + 2))()  # 这个是用来算固定修饰gdm编码的
        for aa_index, aa in enumerate(pep_sq):
            if aa in self.dp.myMOD.fix_mod_dic:
                p_mod_site[aa_index + 1] = self.dp.myMOD.fix_mod_dic[aa]
                fix_mod_code = fix_mod_code | (1 << (aa_index + 1))
            else:
                p_mod_site[aa_index + 1] = -1

        if pep_pos_start == 0:
            if pep_sq[0] in self.dp.myMOD.fix_mod_dic_pro_N:
                one_fix_mod_index = self.dp.myMOD.fix_mod_dic_pro_N[pep_sq[0]]
                p_mod_site[0] = one_fix_mod_index
                fix_mod_code = fix_mod_code | (1 << 0)

        if p_mod_site[0] == -1:
            if pep_sq[0] in self.dp.myMOD.fix_mod_dic_pep_N:
                one_fix_mod_index = self.dp.myMOD.fix_mod_dic_pep_N[pep_sq[0]]
                p_mod_site[0] = one_fix_mod_index
                fix_mod_code = fix_mod_code | (1 << 0)

        if pep_pos_end == pep_length:
            if pep_sq[pep_length - 1] in self.dp.myMOD.fix_mod_dic_pro_C:
                fix_mod_code = fix_mod_code | (1 << (pep_length + 1))
                p_mod_site[pep_length + 1] = self.dp.myMOD.fix_mod_dic_pro_C[pep_sq[pep_length - 1]]

        if p_mod_site[pep_length + 1] == -1:
            if pep_sq[pep_length - 1] in self.dp.myMOD.fix_mod_dic_pep_C:
                fix_mod_code = fix_mod_code | (1 << (pep_length + 1))
                p_mod_site[pep_length + 1] = self.dp.myMOD.fix_mod_dic_pep_C[pep_sq[pep_length - 1]]

        ''' 先不算肽段质量，等去冗余了再算'''
        pep_code = tool_get_peptide_info_code(is_target, pro_index, pep_pos_start, pep_length)
        p_sq = c_char_p(bytes(pep_sq, 'utf-8'))
        gdm = self.__soldier_update_peptide_gdm(bird, p_sq, self.p_G_matrix, pep_length, p_mod_site)
        peptide_set.append([gdm, pep_code, fix_mod_code])
        self.pep_num += 1

    def __soldier_update_peptide_gdm(self, bird, p_sq, p_G_matrix, pep_length, p_mod_site):
        # 获取哥德尔编码值
        return bird.calculate_gdm(p_sq, pep_length, p_mod_site, self.dp.myMOD.fix_mod_num, p_G_matrix, AA_NUM, self.dp.max_sq_len)

    def generate_var_peptide_list(self, set_peptide_dataset, protein_list, peptide_buffer):
        peptide_num = 0
        peptide_buffer_num = 0
        bird = self.import_DLL_function()

        for peptide_index, one_CPeptide in enumerate(set_peptide_dataset):
            one_CPeptide = list(one_CPeptide)
            if (peptide_index + 1) % 20000 == 0:
                print("[Info] Generating peptides : %.4f%%" % (peptide_index / len(set_peptide_dataset) * 100))
            record_var_dic = {}
            is_target, pro_index, pep_start_pos, pep_length = tool_get_peptide_info(one_CPeptide[1])
            pep_end_pos = pep_start_pos + pep_length
            pro_length = len(protein_list.sq[pro_index])
            pep_sq = protein_list.sq[pro_index][pep_start_pos: pep_start_pos + pep_length]
            p_sq = c_char_p(bytes(pep_sq, 'utf-8'))
            # 肽段质量
            peptide_mass = bird.calculate_peptide_mass(p_sq, pep_length, c_double(self.dp.myINI.MASS_H2O), self.p_DIC_AA)
            fix_mod_mass = self.__captain_get_fix_peptide_mass(pep_sq, pep_length, one_CPeptide[2])
            fix_peptide_mass = peptide_mass + fix_mod_mass
            one_CPeptide[0] = fix_peptide_mass
            # 0:gdm/质量 1:肽段获取索引 2:固定修饰编码 3:蛋白索引
            # 可变修饰肽段才有 4:可变修饰数量 5:可变修饰编码
            if self.check_mass(fix_peptide_mass):
                # 这里是先将固定修饰加到候选列表里，不用校验修饰位点了
                self.pep_num += 1
                peptide_buffer_num += 1
                peptide_num += 1
                peptide_buffer.append(tuple(one_CPeptide))
            # 修饰位点：[可变修饰名索引]
            self.__captain_get_var_mod_site(pep_sq, pep_start_pos, pep_end_pos, pep_length, one_CPeptide[2], record_var_dic)
            # 存在可变修饰位点
            if len(record_var_dic.keys()) > 0:
                var_mod_site_combination_list = []
                # 现在有N个要取m个的递归求解
                tool_get_N_combine(list(record_var_dic.keys()), var_mod_site_combination_list, self.dp.myCFG.B6_NUMBER_MAX_MOD)

                for one_var_mod_site_combination in var_mod_site_combination_list:
                    # 与one_var_mod_site_combination对应
                    # var_mod_list的索引号放到one_var_mod_site_combination里得到位点位
                    # var_mod_list的索引号放到var_mod_lis里得到修饰索引
                    candidate_var_mod_list = []
                    self.__captain_get_var_mod_list(one_var_mod_site_combination, record_var_dic, 0, [], candidate_var_mod_list)

                    for var_mod_list in candidate_var_mod_list:
                        new_var_pep = copy.deepcopy(one_CPeptide)

                        var_mod_add_mass = 0.0
                        cur_pep_var_mod_num = 0
                        var_mod_code = 0
                        var_mod_code_list = []
                        pep_N_aa_mod = None
                        pep_C_aa_mod = None
                        for var_mod_list_index, var_mod_index in enumerate(var_mod_list):
                            one_site = one_var_mod_site_combination[var_mod_list_index]  # 修饰在的位点
                            var_mod_index_total = var_mod_index + self.dp.myMOD.fix_mod_num  # 修饰在所有添加修饰你中的索引
                            var_mod_add_mass += self.dp.myINI.DIC_MOD[self.dp.myMOD.var_mod_list[var_mod_index]].mass
                            cur_pep_var_mod_num += 1
                            var_mod_code = tool_update_var_mod_code_by_site_data(var_mod_code, self.dp.myMOD.fix_mod_num, cur_pep_var_mod_num, one_site, var_mod_index_total)
                            # 修饰是肽段两端的修饰没有影响，主要是检查在两端氨基酸上的修饰
                            if one_site == 1:
                                pep_N_aa_mod = var_mod_index_total
                            if one_site == pep_length:
                                pep_C_aa_mod = var_mod_index_total
                            if cur_pep_var_mod_num % 5 == 0:
                                var_mod_code_list.append(var_mod_code)
                                var_mod_code = 0
                        if var_mod_code != 0:
                            var_mod_code_list.append(var_mod_code)
                        new_var_pep[0] = new_var_pep[0] + var_mod_add_mass
                        new_var_pep.append(cur_pep_var_mod_num)
                        new_var_pep.append(tuple(var_mod_code_list))
                        if self.check_mass(new_var_pep[0]) and self.check_mod_site(pep_sq, pep_start_pos, pep_end_pos, pro_length, pep_N_aa_mod, pep_C_aa_mod):
                            peptide_buffer.append(tuple(new_var_pep))
                            peptide_buffer_num += 1
                            peptide_num += 1
                        del var_mod_code
                        del new_var_pep
        print("[Info] Generating peptides : 100%")
        print("[Info] Adding variable modification num : %d" % peptide_num)

    def __captain_get_fix_peptide_mass(self, pep_sq, pep_length, fix_mod_code):
        fix_mod_mass = 0.0
        for site in range(pep_length + 2):
            site_data = fix_mod_code & 1
            if site_data == 0:
                pass
            else:
                if site == 0:
                    if pep_sq[0] in self.dp.myMOD.fix_mod_dic_pro_N.keys():
                        fix_mod_mass += self.dp.myINI.DIC_MOD[self.dp.myMOD.fix_mod_list[self.dp.myMOD.fix_mod_dic_pro_N[pep_sq[0]]]].mass
                    else:
                        if pep_sq[0] in self.dp.myMOD.fix_mod_dic_pep_N.keys():
                            fix_mod_mass += self.dp.myINI.DIC_MOD[self.dp.myMOD.fix_mod_list[self.dp.myMOD.fix_mod_dic_pep_N[pep_sq[0]]]].mass
                elif site == pep_length + 1:
                    if pep_sq[-1] in self.dp.myMOD.fix_mod_dic_pro_C.keys():
                        fix_mod_mass += self.dp.myINI.DIC_MOD[self.dp.myMOD.fix_mod_list[self.dp.myMOD.fix_mod_dic_pro_C[pep_sq[-1]]]].mass
                    else:
                        if pep_sq[-1] in self.dp.myMOD.fix_mod_dic_pep_C.keys():
                            fix_mod_mass += self.dp.myINI.DIC_MOD[self.dp.myMOD.fix_mod_list[self.dp.myMOD.fix_mod_dic_pep_C[pep_sq[-1]]]].mass
                        else:
                            pass
                else:
                    fix_mod_mass += self.dp.myINI.DIC_MOD[self.dp.myMOD.fix_mod_list[self.dp.myMOD.fix_mod_dic[pep_sq[site - 1]]]].mass
            fix_mod_code = fix_mod_code >> 1
        return fix_mod_mass

    def __captain_get_var_mod_site(self, pep_sq, pep_start_pos, pep_end_pos, pep_length, fix_mod_code, record_var_dic):
        pep_C_fix_mod = fix_mod_code & 1
        pep_N_fix_mod = (fix_mod_code >> pep_length) & 1
        # 记录可变修饰两端
        if pep_C_fix_mod == 0:
            if pep_start_pos == 0:
                if pep_sq[0] in self.dp.myMOD.var_mod_dic_pro_N:
                    record_var_dic[0] = self.dp.myMOD.var_mod_dic_pro_N[pep_sq[0]]

            if pep_sq[0] in self.dp.myMOD.var_mod_dic_pep_N:
                if 0 in record_var_dic.keys():
                    record_var_dic[0] += self.dp.myMOD.var_mod_dic_pep_N[pep_sq[0]]
                else:
                    record_var_dic[0] = self.dp.myMOD.var_mod_dic_pep_N[pep_sq[0]]
        if pep_N_fix_mod == 0:
            if pep_end_pos == pep_length:
                if pep_sq[pep_length - 1] in self.dp.myMOD.var_mod_dic_pro_C:
                    record_var_dic[pep_length + 1] = self.dp.myMOD.var_mod_dic_pro_C[pep_sq[pep_length - 1]]

            if pep_sq[pep_length - 1] in self.dp.myMOD.var_mod_dic_pep_C:
                if pep_length + 1 in record_var_dic.keys():
                    record_var_dic[pep_length + 1] += self.dp.myMOD.var_mod_dic_pro_C[pep_sq[pep_length - 1]]
                else:
                    record_var_dic[pep_length + 1] = self.dp.myMOD.var_mod_dic_pro_C[pep_sq[pep_length - 1]]

        for aa_index, aa in enumerate(pep_sq):
            if aa in self.dp.myMOD.var_mod_dic:
                if (fix_mod_code >> (aa_index + 1)) & 1 == 0:
                    record_var_dic[aa_index + 1] = self.dp.myMOD.var_mod_dic[aa]

    def __captain_get_var_mod_list(self, var_site_combination, record_var_dic, site_index, cur_record, out_record):
        '''
        :param var_site_combination: # 这条肽段存在的修饰位置[site1, site2, site3]
        :param record_var_dic: 表示这条肽段存在的可变修饰位点字典 {site1: [var_mod_index], site2: [var_mod_index]}
        :param site_index: var_site_combination的起始索引
        :param cur_record:递归记录的list
        :param out_record:最后返回的list
        :return:
        '''
        if site_index == len(var_site_combination):
            out_record.append(cur_record)
        else:
            site_list = record_var_dic[var_site_combination[site_index]]

            for one_site_mod in site_list:
                tmp_record = copy.copy(cur_record)
                tmp_record.append(one_site_mod)
                self.__captain_get_var_mod_list(var_site_combination, record_var_dic, site_index + 1, tmp_record, out_record)

    def check_mass(self, mass):
        # check the mass of peptide sequence is in [min_mass, max_mass]
        return mass >= self.dp.myCFG.D6_MASS_PEP_LOW and mass <= self.dp.myCFG.D7_MASS_PEP_UP

    def check_length(self, start_site, end_site):
        # check the length of peptide sequence is in [min_len, max_len]
        len_peptide_seq = end_site - start_site
        return len_peptide_seq >= self.dp.myCFG.D8_LEN_PEP_LOW and len_peptide_seq <= self.dp.myCFG.D9_LEN_PEP_UP

    def check_mod_site(self, pep_sq, pep_start_pos, pep_end_pos, pro_length, pep_N_aa_mod, pep_C_aa_mod):
        # 这个是检查修饰是否可以产生在肽段两端的氨基酸上的
        return_value = True
        # 因为正常的修饰在肽段产生的时候就做好处理了，这里主要是过滤mono-link的修饰
        # 修饰在肽段两端氨基酸需要专门判断一下
        # pep_N_mod = mod_index, pep_C_mod = mod_index
        if self.dp.myENZ.aa_enzyme_C[ord(pep_sq[-1]) - ord('A')]:
            # C端酶切
            if pep_end_pos == pro_length:
                # 蛋白质C端发生什么修饰都不影响酶切，对上氨基酸即可
                pass
            else:
                # 这个肽段C端的氨基酸是C端酶切的位点
                if pep_C_aa_mod is not None:
                    # 肽段C端有修饰不是C端的修饰那么就不能出这个肽段
                    if pep_C_aa_mod >= self.dp.myMOD.fix_mod_num + self.dp.myMOD.var_mod_num:
                        # 这是mono-link修饰
                        if pep_C_aa_mod - self.dp.myMOD.fix_mod_num - self.dp.myMOD.var_mod_num == 0:
                            if 1 not in self.dp.myLINK.linker_1_data.linker_site:
                                return_value = False
                        elif pep_C_aa_mod - self.dp.myMOD.fix_mod_num - self.dp.myMOD.var_mod_num == 1:
                            if 1 not in self.dp.myLINK.linker_2_data.linker_site:
                                return_value = False
                        else:
                            return_value = False
        if self.dp.myENZ.aa_enzyme_N[ord(pep_sq[0]) - ord('A')]:
            if pep_start_pos == 0:
                # 蛋白质N端发生什么修饰都不影响酶切，对上氨基酸即可
                pass
            else:
                # 这个肽段N端的氨基酸是N端酶切的位点,那就必须得是肽段N端修饰才可以
                if pep_N_aa_mod is not None:
                    if pep_N_aa_mod >= self.dp.myMOD.fix_mod_num + self.dp.myMOD.var_mod_num:
                        # 这是mono-link修饰
                        if pep_N_aa_mod - self.dp.myMOD.fix_mod_num - self.dp.myMOD.var_mod_num == 0:
                            if -1 not in self.dp.myLINK.linker_1_data.linker_site:
                                return_value = False
                            elif pep_N_aa_mod - self.dp.myMOD.fix_mod_num - self.dp.myMOD.var_mod_num == 1:
                                if -1 not in self.dp.myLINK.linker_2_data.linker_site:
                                    return_value = False
                            else:
                                return_value = False
        return return_value
