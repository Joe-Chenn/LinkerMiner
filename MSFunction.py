# -*- mode: python ; coding: utf-8 -*-..
from MSOperator import op_INIT_CFILE_MS2, op_modify_timsReaderconfig
from MSTool import tool_compute_mass_from_element, tool_get_peptide_mass, tool_compute_mass_from_element_for_isotope
from MSData import CModData, CProtein, CLinkerData, CTimsTOFReaderConfig
from MSSysterm import FOLDER_INDEX, FOLDER_MASS, AA_NUM, PRO_PKL_FILE, FILE_NAME_PSM_RESULT, SITE_DIC
from MSSysterm import VALUE_MAX_SCAN, TimsTOFReadConfigName
from MSSysterm import PPARSE_EXE, TIMSTOFREADER_EXE
from MSLogging import logGetError, logGetWarning
# from MSSysterm import ATOM_MASS_P, MASS_PKL

from MSData import CFileMS2

import numpy as np
import os
import pickle
import platform
import math
import copy
import operator
import time
import sys

class CFunctionConfigIO:

    def config2file(self, path, config):

        with open(path, 'w') as f:
            f.write("#[data]\n")
            f.write("TYPE_MS2=" + config.A1_TYPE_MS2 + "\n")
            f.write("PATH_MS2=" + config.A2_PATH_MS2 + "\n")
            f.write("PATH_FASTA=" + config.A3_PATH_FASTA + "\n")
            f.write("PATH_FASTA_EXPORT=" + config.A4_PATH_FASTA_EXPORT + "\n")
            f.write("PATH_RESULT_EXPORT=" + config.A5_PATH_RESULT_EXPORT + "\n")
            f.write("\n")

            f.write("#[biology]\n")
            f.write("NAME_ENZYME=" + config.B1_NAME_ENZYME + "\n")
            f.write("TYPE_DIGEST=" + config.B2_TYPE_DIGEST + "\n")
            f.write("NUMBER_MAX_MISS_CLV=" + config.B3_NUMBER_MAX_MISS_CLV + "\n")
            f.write("NAME_MOD_FIX=" + config.B4_NAME_MOD_FIX + "\n")
            f.write("NAME_MOD_VAR=" + config.B5_NAME_MOD_VAR + "\n")
            f.write("NUMBER_MAX_MOD=" + config.B6_NUMBER_MAX_MOD + "\n")
            f.write("LINKER_1_NAME=" + config.B7_LINKER_1_NAME + "\n")
            f.write("LINKER_2_NAME=" + config.B8_LINKER_2_NAME + "\n")
            f.write("MAX_PEPNUM=" + config.B9_MAX_PEPNUM + '\n')
            f.write("\n")

            f.write("#[mass spectrometry]\n")
            f.write("PPM_TOL_PRECURSOR=" + config.C1_PPM_TOL_PRECURSOR + "\n")
            f.write("PPM_TOL_FRAGMENT=" + config.C2_PPM_TOL_FRAGMENT + "\n")
            f.write("TYPE_ACTIVATION=" + config.C3_activation_type + "\n")
            f.write("\n")

            f.write("#[performance]\n")
            f.write("NUMBER_THREAD=" + config.D1_NUMBER_THREAD + "\n")
            f.write("TYPE_THREAD=" + config.D2_TYPE_THREAD + "\n")
            f.write("NUMBER_SELECT_PEAK=" + config.D3_NUMBER_SELECT_PEAK + "\n")
            f.write("NUMBER_SPECTRUM=" + config.D4_NUMBER_SPECTRUM + "\n")
            f.write("LEN_MAX_PROTEIN=" + config.D5_LEN_MAX_PROTEIN + "\n")
            f.write("MASS_PEP_LOW=" + config.D6_MASS_PEP_LOW + "\n")
            f.write("MASS_PEP_UP=" + config.D7_MASS_PEP_UP + "\n")
            f.write("LEN_PEP_LOW=" + config.D8_LEN_PEP_LOW + "\n")
            f.write("LEN_PEP_UP=" + config.D9_LEN_PEP_UP + "\n")
            f.write("INDEX_SPLIT_MASS=" + config.D10_INDEX_SPLIT_MASS + "\n")
            f.write("NUMBER_TOP_RESULT=" + config.D11_NUMBER_TOP_RESULT + "\n")
            f.write("\n")
            f.write("MULTI_MASS=" + config.D12_MULTI_MASS + "\n")
            f.write("TYPE_FLOW=" + config.D13_TYPE_FLOW + "\n")
            f.write("TYPE_PEPTIDE=" + config.D14_TYPE_PEPTIDE + '\n')
            f.write("TYPE_LINK=" + config.D15_TYPE_LINK + '\n')
            f.write("\n")
            f.write("#[filter]\n")
            f.write("FDR_PSM=" + config.E1_FDR_PSM + "\n")
            f.write("FDR_SEPARATE=" + config.E2_FDR_SEPARATE + "\n")
            f.write("MOBILITY_RERANK=" + config.E3_MOBILITY_RERANK + "\n")
            f.write("\n")
            f.write("#[ini]\n")
            f.write("PATH_INI_ELEMENT=" + config.F1_PATH_INI_ELEMENT + "\n")
            f.write("PATH_INI_AA=" + config.F2_PATH_INI_AA + "\n")
            f.write("PATH_INI_MOD=" + config.F3_PATH_INI_MOD + "\n")
            f.write("PATH_INI_XLINK=" + config.F4_PATH_INI_XLINK + "\n")

    def file2config(self, file_path, config):

        sys_flag = platform.system()
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.find('=') < 0:
                    continue
                if line.find('#') >= 0:
                    line = line[:line.find('#')].strip()
                tmp = line.split('=')
                if len(tmp) < 2:
                    continue
                name = tmp[0]
                value = tmp[1]
                print(name, ":", value)
                try:
                    if name == "PATH_INI_ELEMENT":
                        config.F1_PATH_INI_ELEMENT = value
                    elif name == "PATH_INI_AA":
                        try:
                            config.F2_PATH_INI_AA = value
                        except Exception as err:
                            print("[Error] The element_file must be in front of AA_file.", err)
                            info = "[Error] The element_file must be in front of AA_file."
                            logGetError(info)
                    elif name == "PATH_INI_MOD":
                        config.F3_PATH_INI_MOD = value
                    elif name == "PATH_INI_XLINK":
                        config.F4_PATH_INI_XLINK = value

                    # [data]
                    elif name == "TYPE_MS2":
                        config.A1_TYPE_MS2 = value
                    elif name == "PATH_MS2":
                        if sys_flag == "Linux":
                            value = value.replace('\\', '/')
                        else:
                            value = value.replace('/', '\\')
                        config.A2_PATH_MS2 = value
                    elif name == "PATH_FASTA":
                        if sys_flag == "Linux":
                            value = value.replace('\\', '/')
                        else:
                            value = value.replace('/', '\\')
                        config.A3_PATH_FASTA = value
                    elif name == "PATH_FASTA_EXPORT":
                        if sys_flag == "Linux":
                            value = value.replace('\\', '/')
                        else:
                            value = value.replace('/', '\\')
                        config.A4_PATH_FASTA_EXPORT = value
                        if not os.path.exists(config.A4_PATH_FASTA_EXPORT):
                            os.makedirs(config.A4_PATH_FASTA_EXPORT)

                    elif name == "PATH_RESULT_EXPORT":
                        if sys_flag == "Linux":
                            value = value.replace('\\', '/')
                        else:
                            value = value.replace('/', '\\')
                        config.A5_PATH_RESULT_EXPORT = value
                        if not os.path.exists(config.A5_PATH_RESULT_EXPORT):
                            os.makedirs(config.A5_PATH_RESULT_EXPORT)

                    # [biology]
                    elif name == "NAME_ENZYME":
                        config.B1_NAME_ENZYME = value
                    elif name == "TYPE_DIGEST":
                        config.B2_TYPE_DIGEST = int(value)
                    elif name == "NUMBER_MAX_MISS_CLV":
                        config.B3_NUMBER_MAX_MISS_CLV = int(value)
                    elif name == "NAME_MOD_FIX":
                        config.B4_NAME_MOD_FIX = value
                    elif name == "NAME_MOD_VAR":
                        config.B5_NAME_MOD_VAR = value
                    elif name == "NUMBER_MAX_MOD":
                        config.B6_NUMBER_MAX_MOD = int(value)
                    elif name == 'LINKER_1_NAME':
                        config.B7_LINKER_1_NAME = value
                    elif name == 'LINKER_2_NAME':
                        config.B8_LINKER_2_NAME = value
                    elif name == 'MAX_PEPNUM':
                        config.B9_MAX_PEPNUM = int(value)

                    # [mass spectrometry]
                    elif name == "PPM_TOL_PRECURSOR":
                        config.is_ppm_pre = True
                        if value.endswith("ppm") or value.endswith("PPM"):
                            config.C1_PPM_TOL_PRECURSOR = float(value[:-3]) * 1e-6
                        elif value.endswith("da") or value.endswith("DA") or value.endswith("Da"):
                            config.is_ppm_pre = False
                            config.C1_PPM_TOL_PRECURSOR = float(value[:-2])
                        else:
                            config.C1_PPM_TOL_PRECURSOR = 20e-6  # default is 20ppm
                    elif name == "PPM_TOL_FRAGMENT":
                        config.is_ppm_fra = True
                        if value.endswith("ppm") or value.endswith("PPM"):
                            config.C2_PPM_TOL_FRAGMENT = float(value[:-3]) * 1e-6
                        elif value.endswith("da") or value.endswith("DA") or value.endswith("Da"):
                            config.is_ppm_fra = False
                            config.C2_PPM_TOL_FRAGMENT = float(value[:-2])
                        else:
                            config.C2_PPM_TOL_FRAGMENT = 20e-6  # default is 20ppm

                    # [performance]
                    elif name == "NUMBER_THREAD":
                        config.D1_NUMBER_THREAD = int(value)
                    elif name == "TYPE_THREAD":
                        config.D2_TYPE_THREAD = int(value)
                    elif name == "NUMBER_SELECT_PEAK":
                        config.D3_NUMBER_SELECT_PEAK = int(value)
                    elif name == "NUMBER_SPECTRUM":
                        config.D4_NUMBER_SPECTRUM = int(value)
                    elif name == 'LEN_MAX_PROTEIN':
                        config.D5_LEN_MAX_PROTEIN = int(value)
                    elif name == "MASS_PEP_LOW":
                        config.D6_MASS_PEP_LOW = float(value)
                        config.start_mass = math.ceil(float(value) / 100) * 100
                    elif name == "MASS_PEP_UP":
                        config.D7_MASS_PEP_UP = float(value)
                        config.end_mass = int(float(value) / 100) * 100
                    elif name == "LEN_PEP_LOW":
                        config.D8_LEN_PEP_LOW = int(value)
                    elif name == "LEN_PEP_UP":
                        config.D9_LEN_PEP_UP = int(value)
                    elif name == "INDEX_SPLIT_MASS":
                        config.D10_INDEX_SPLIT_MASS = float(value)
                    elif name == "NUMBER_TOP_RESULT":
                        config.D11_NUMBER_TOP_RESULT = int(value)
                    elif name == "MULTI_MASS":
                        config.D12_MULTI_MASS = int(value)
                    elif name == "TYPE_FLOW":
                        config.D13_TYPE_FLOW = int(value)
                    elif name == 'TYPE_PEPTIDE':
                        config.D14_TYPE_PEPTIDE = int(value)
                    elif name == 'TYPE_LINK':
                        config.D15_TYPE_LINK = int(value)
                    elif name == 'FILTER_LOOP':
                        config.D16_FILTER_LOOP = int(value)
                    elif name == "FDR_PSM":
                        config.E1_FDR_PSM = float(value)
                    elif name == "FDR_SEPARATE":
                        if value == '1':
                            config.E2_FDR_SEPARATE = True
                        else:
                            config.E2_FDR_SEPARATE = False
                    elif name == "MOBILITY_RERANK":
                        if value == '1':
                            config.E3_MOBILITY_RERANK = True
                        else:
                            config.E3_MOBILITY_RERANK = False
                    elif name == "THRESHOLD_SCORE":
                        config.score_t = float(value)
                    elif name == "NUMBER_CAND":
                        config.cand_num = int(value)

                except Exception as err:
                    print("[Error] occured in configure.txt", err)
                    # info = "[Error] occured in configure.txt"
                    # logGetError(info)

class CFunctionPickle:
    def __init__(self):
        pass

    def write_data_to_pkl(self, inputData, path):
        f = open(path, 'wb')
        pickle.dump(inputData, f)
        f.close()

    def load_pkl_to_data(self, path):
        f = open(path, 'rb')
        load_data = pickle.load(f)
        f.close()
        return load_data

    def get_start_end_mass(self, path):
        # From pkl file path to get the start and end mass
        if path.find('/') >= 0:
            filename = path[path.rfind('/')+1:]
        else:
            filename = path[path.rfind('\\')+1:]
        filename = filename[:filename.find('.')]
        index = filename.find('_')
        start_mass = float(filename[:index])
        end_mass = float(filename[index+1:])
        return start_mass, end_mass

    def get_all_pkl_files(self, folder):
        pkl_file_list, ind_file_list = [], []
        for p in os.listdir(folder):
            if not p.endswith('.pkl'): continue
            if p.endswith("_ind.pkl"):
                pass
            else:
                pkl_file_list.append(os.path.join(folder, p))
                ind_file_list.append(os.path.join(folder, p[:-4] + "_ind.pkl"))
        assert(len(pkl_file_list) == len(ind_file_list))
        return pkl_file_list, ind_file_list

class CFunctionLoadINIFile:

    def __init__(self, inputDP):
        self.dp = inputDP

    def get_INI_data(self):
        self.__captainLoad_element_file()
        self.__captainLoad_aa_file()
        self.__captainLoad_mod_file()
        self.__captainLoad_xlink_file()

    def get_protein_data(self):

        functionPickle = CFunctionPickle()

        protein_list = CProtein()
        ac, de, sq = "", "", ""
        list_ac = []
        list_de = []
        list_sq = []
        list_sq_original = []
        with open(self.dp.myCFG.A3_PATH_FASTA, 'r') as f:
            for line in f:
                line = line.strip()
                if line == "":
                    continue
                if line[0] == '>':
                    if ac != "":
                        sq_replace = sq.replace('I', 'L')
                        list_ac.append(ac)
                        list_de.append(de)
                        list_sq.append(sq_replace)
                        list_sq_original.append(sq)
                        if self.dp.myCFG.using_decoy:
                            decoy_pro = self.__captain_generate_decoy_protein(ac, sq)
                            list_ac.append(decoy_pro[0])
                            list_de.append(decoy_pro[1])
                            list_sq.append(decoy_pro[2])
                            list_sq_original.append(decoy_pro[3])
                    line = line[1:]
                    ind = line.find(' ')
                    if ind >= 0:
                        ac = line[:ind]
                    else:
                        ac = line
                    de = line[ind + 1:]
                    sq = ""
                else:
                    sq += line
        if ac != "":
            sq_replace = sq.replace('I', 'L')
            list_ac.append(ac)
            list_de.append(de)
            list_sq.append(sq_replace)
            list_sq_original.append(sq)
            if self.dp.myCFG.using_decoy:
                decoy_pro = self.__captain_generate_decoy_protein(ac, sq)
                list_ac.append(decoy_pro[0])
                list_de.append(decoy_pro[1])
                list_sq.append(decoy_pro[2])
                list_sq_original.append(decoy_pro[3])
        protein_list.ac = protein_list.ac + list_ac
        protein_list.de = protein_list.de + list_de
        protein_list.sq = protein_list.sq + list_sq
        protein_list.sq_original = protein_list.sq_original + list_sq_original

        functionPickle.write_data_to_pkl(protein_list, self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + PRO_PKL_FILE)

    def __captainLoad_xlink_file(self):
        with open(self.dp.myCFG.F4_PATH_INI_XLINK, 'r') as f:
            for line in f:
                if line.startswith('['):
                    continue
                elif line.startswith('total'):
                    continue
                elif line.startswith('name'):
                    continue
                elif line.strip() == '':
                    continue
                else:
                    link_name = line.split('=')[0]
                    link_data = line.split('=')[1]
                    link_data.strip()

                    list_link_data = link_data.split(' ')

                    alpha_site_str = list_link_data[0]
                    beta_site_str = list_link_data[1]

                    link_site = []
                    alpha_site_list = []
                    beta_site_list = []
                    same_site_list = []

                    self.__soldier_set_link_site_data(alpha_site_str, beta_site_str, link_site, alpha_site_list,
                                                      beta_site_list, same_site_list)

                    alpha_site = copy.deepcopy(SITE_DIC)
                    self.__soldier_get_site_dic(alpha_site_list, alpha_site)
                    beta_site = copy.deepcopy(SITE_DIC)
                    self.__soldier_get_site_dic(beta_site_list, beta_site)
                    same_site = copy.deepcopy(SITE_DIC)
                    self.__soldier_get_site_dic(same_site_list, same_site)

                    dic_site = copy.deepcopy(SITE_DIC)
                    self.__soldier_get_site_dic(alpha_site_list, dic_site)
                    self.__soldier_get_site_dic(beta_site_list, dic_site)
                    self.__soldier_get_site_dic(same_site_list, dic_site)

                    cross_mass = self.__soldier_check_mass(link_name, list_link_data[6], float(list_link_data[2]))
                    loop_mass = self.__soldier_check_mass(link_name, list_link_data[6], float(list_link_data[3]))
                    mono_mass = self.__soldier_check_mass(link_name, list_link_data[7], float(list_link_data[4]))

                    self.dp.myINI.DIC_LINKER[link_name] = CLinkerData(link_name, link_site, dic_site, alpha_site, beta_site,
                                                                      same_site, cross_mass, loop_mass, mono_mass, list_link_data[6])

    def __soldier_check_mass(self, linker_name, intput_molecular_formula, input_mass):

        # 这里是做一下检验，怕分子式和填写的数值不一致
        cal_mass = tool_compute_mass_from_element_for_isotope(intput_molecular_formula, self.dp.myINI.DIC_ELEMENT_MASS)
        if abs(cal_mass - input_mass) < 0.001:
            out_mass = input_mass
        else:
            out_mass = cal_mass
            # info = linker_name + " mass is different from the calculated mass !"
            # logGetWarning(info)
        return out_mass

    def __soldier_set_link_site_data(self, alpha_site_str, beta_site_str, all_stie_list, alpha_site_list, beta_site_list, same_site_list):
        for item in alpha_site_str:
            if item in beta_site_str:
                if item == '[':
                    find_key = -2
                elif item == '{':
                    find_key = -1
                elif item == ']':
                    find_key = 2
                elif item == '}':
                    find_key = 1
                else:
                    find_key = item
                if find_key in same_site_list:
                    pass
                else:
                    same_site_list.append(find_key)
            else:
                if item == '[':
                    find_key = -2
                elif item == '{':
                    find_key = -1
                elif item == ']':
                    find_key = 2
                elif item == '}':
                    find_key = 1
                else:
                    find_key = item
                if find_key in alpha_site_list:
                    pass
                else:
                    alpha_site_list.append(find_key)
            if find_key in all_stie_list:
                if find_key in beta_site_list:
                    del beta_site_list[beta_site_list.index(find_key)]
                if find_key in alpha_site_list:
                    del alpha_site_list[alpha_site_list.index(find_key)]
            else:
                all_stie_list.append(find_key)

        for item in beta_site_str:
            if item == '[':
                find_key = -2
            elif item == '{':
                find_key = -1
            elif item == ']':
                find_key = 2
            elif item == '}':
                find_key = 1
            else:
                find_key = item
            if find_key in same_site_list:
                pass
            else:
                beta_site_list.append(find_key)
            if find_key in all_stie_list:
                if find_key in beta_site_list:
                    del beta_site_list[beta_site_list.index(find_key)]
                if find_key in alpha_site_list:
                    del alpha_site_list[alpha_site_list.index(find_key)]
            else:
                all_stie_list.append(find_key)

    def __soldier_get_site_dic(self, input_site_list, out_site_dic):

        for item in input_site_list:
            if item not in out_site_dic:
                info = "A linker site may not exist!"
                logGetWarning(info)
            else:
                out_site_dic[item] = 1

    def __captainLoad_element_file(self):
        with open(self.dp.myCFG.F1_PATH_INI_ELEMENT, 'r') as f:
            for line in f:
                if line.find('=') < 0:
                    continue
                if line.strip() == "" or line[0] == '@':
                    continue
                line = line.strip()
                line = line[line.find('=')+1:]
                tmp = line.split('|')
                element_name = tmp[0]
                element_mass = [float(v) for v in tmp[1][:-1].split(',')]
                element_percentage = [float(v) for v in tmp[2][:-1].split(',')]
                assert(len(element_mass) == len(element_percentage))
                index = element_percentage.index(max(element_percentage))
                self.dp.myINI.DIC_ELEMENT_MASS[element_name] = element_mass[index]
                self.dp.myINI.DIC_ELEMENT_ALL_MASS[element_name] = element_mass
                self.dp.myINI.DIC_ELEMENT_MASS_PERCENTAGE[element_name] = element_percentage


    def __captainLoad_aa_file(self):
        with open(self.dp.myCFG.F2_PATH_INI_AA, 'r') as f:
            for line in f:
                if line.find('=') < 0:
                    continue
                if line.strip() == "" or line[0] == '@':
                    continue
                name, value = line.strip().split('=')
                aa, elem_str = value[:-1].split('|')
                self.dp.myINI.DIC_AA[aa] = tool_compute_mass_from_element(elem_str, self.dp.myINI.DIC_ELEMENT_MASS)
                self.dp.myINI.DIC_AA_COMPOSITION[aa] = elem_str

    def __captainLoad_mod_file(self):
        with open(self.dp.myCFG.F3_PATH_INI_MOD, 'r') as f:
            for line in f:
                if line.find('=') < 0: continue
                if line.strip() == "" or line[0] == '@': continue
                if line.startswith("name"): continue
                name, value = line.strip().split('=')
                tmp = value.split(' ')
                mass = float(tmp[2])
                aa = tmp[0]
                composition = ""
                if len(tmp) > 6:
                    composition = tmp[7]
                else:
                    composition = tmp[5]
                if tmp[1] == "PEP_N":
                    type = -1
                elif tmp[1] == "PEP_C":
                    type = 1
                elif tmp[1] == "PRO_N":
                    type = -2
                elif tmp[1] == "PRO_C":
                    type = 2
                elif tmp[1] == "NORMAL":
                    type = 0
                else:
                    return -1
                self.dp.myINI.DIC_MOD[name] = CModData(name, mass, type, aa, composition)

    def __captain_generate_decoy_protein(self, ac, sq):
        new_ac = "REV_" + ac
        new_sq = sq[::-1]
        new_sq_replace = new_sq.replace('I', 'L')
        decoy_protein = (new_ac, "", new_sq_replace, new_sq)  # 蛋白名 描述 序列
        return decoy_protein

class CFunctionFillConfigData:

    def __init__(self, inputDP):

        self.dp = inputDP
        self.dp.myCFG.max_mz_index_size = int(self.dp.myCFG.max_mz * self.dp.myCFG.D10_INDEX_SPLIT_MASS)

    def generate_mod_config(self):

        ret = self.__captain_init_mod()
        if ret != 0:
            print("[Error] occured in _init_mod")
            info = "[Error] occured in configure.txt"
            logGetError(info)

    def generate_link_config(self):
        self.__captain_init_link()

    def generate_ms2_config(self):
        self.__captain_init_ms2()

    def generate_G_matrix(self):

        var_mod = self.dp.myCFG.B5_NAME_MOD_VAR.split(';')
        prime_num = self.dp.myCFG.D9_LEN_PEP_UP + self.dp.myCFG.B6_NUMBER_MAX_MOD + 1
        aa_num = AA_NUM + len(var_mod) + 1

        for i in range(aa_num):
            tmp = [0.0 for j in range(prime_num)]
            self.dp.G_matrix.append(tmp)
        prime_list = [0 for j in range(prime_num)]
        prime_list[0] = 2
        for i in range(1, prime_num):
            start_v = prime_list[i - 1] + 1
            while True:
                is_prime = True
                v = int(math.sqrt(start_v))
                if v == 1: v = 2
                for k in range(v, start_v):
                    if start_v % k == 0:
                        is_prime = False
                        break
                if is_prime: break
                start_v += 1
            prime_list[i] = start_v

        for i in range(aa_num):
            for j in range(prime_num):
                self.dp.G_matrix[i][j] = (i + 1) * 1.0 * math.log(prime_list[j])
        self.dp.aa_num = aa_num
        self.dp.max_sq_len = prime_num

    def generate_enzyme_config(self):
        value_list = self.dp.myCFG.B1_NAME_ENZYME.split(';')
        for v in value_list:
            if len(v.split(' ')) < 3:
                info = "[Error] Enzyme data may be wrong !"
                logGetError(info)
            tn, ta, tf = v.split(' ')
            self.dp.myENZ.enzyme_name.append(tn)
            self.dp.myENZ.enzyme_aa.append(ta)
            self.dp.myENZ.enzyme_flag.append(tf)
        self.__captain_init_aa_enzyme_CN()

    def __captain_init_aa_enzyme_CN(self):

        self.dp.myENZ.aa_enzyme_C = [False for i in range(AA_NUM)]
        self.dp.myENZ.aa_enzyme_N = [False for i in range(AA_NUM)]

        for i, aa in enumerate(self.dp.myENZ.enzyme_aa):
            for one_aa in aa:
                flag = self.dp.myENZ.enzyme_flag[i]
                if flag == "C":
                    self.dp.myENZ.aa_enzyme_C[ord(one_aa) - ord('A')] = True
                    self.dp.myENZ.enzyme_C = True
                elif flag == "N":
                    self.dp.myENZ.aa_enzyme_N[ord(one_aa) - ord('A')] = True
                    self.dp.myENZ.enzyme_N = True

    def generate_massindex_config(self):
        pkl_num = int((self.dp.myCFG.end_mass - self.dp.myCFG.start_mass) / 100) + 1
        # 默认按100Da区间建pkl文件
        for i in range(pkl_num):
            start_mass = self.dp.myCFG.start_mass + i * 100
            end_mass = self.dp.myCFG.start_mass + (i + 1) * 100
            if end_mass > self.dp.myCFG.end_mass:
                end_mass = self.dp.myCFG.end_mass
            if start_mass == end_mass:
                continue
            pkl_file = os.path.join(self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_MASS + '\\', "{0}_{1}.pkl".format(int(start_mass), int(end_mass)))
            ind_file = pkl_file[:pkl_file.rfind('.')] + "_ind.pkl"
            self.dp.LIST_MASS_INDEX_FILE.append(pkl_file)
            self.dp.LIST_MASS_IND_FILE.append(ind_file)

    def add_mono_var_mod(self, linker_name, complex_search=False):

        self.__soldier_init_mono_mod(linker_name, complex_search)

    def __captain_init_link(self):

        if self.dp.myCFG.B7_LINKER_1_NAME in self.dp.myINI.DIC_LINKER.keys():
            self.dp.myLINK.linker_1_data = self.dp.myINI.DIC_LINKER[self.dp.myCFG.B7_LINKER_1_NAME]
        else:
            logGetError('The LINKER_1_NAME may be wrong !')
        if len(self.dp.myCFG.B8_LINKER_2_NAME) == 0:
            pass
        else:
            if self.dp.myCFG.B8_LINKER_2_NAME in self.dp.myINI.DIC_LINKER.keys():
                self.dp.myLINK.linker_2_data = self.dp.myINI.DIC_LINKER[self.dp.myCFG.B8_LINKER_2_NAME]
            else:
                logGetError('The LINKER_2_NAME may be wrong !')
        if self.dp.myCFG.B7_LINKER_1_NAME == self.dp.myCFG.B8_LINKER_2_NAME:
            self.dp.myLINK.same_linker = 1

    def __captain_init_ms2(self):
        # 这里是把那些文件名都给弄好
        if self.dp.myCFG.A1_TYPE_MS2 in ['mgf']:

            if not os.path.isdir(self.dp.myCFG.A2_PATH_MS2) and "|" not in self.dp.myCFG.A2_PATH_MS2:
                self.dp.myMS2.ms2_num = 1
                self.dp.myMS2.original_file_name.append(self.dp.myCFG.A2_PATH_MS2)
                self.dp.myMS2.ms2_file.append(self.dp.myCFG.A2_PATH_MS2)
                file_name = self.dp.myCFG.A2_PATH_MS2.split('.')[0]
                self.dp.myMS2.original_file_name.append(self.dp.myCFG.A2_PATH_MS2)
                self.dp.myMS2.ms2_file.append(file_name + '.mgf')
                self.dp.myMS2.trans_ms2_file.append(file_name + '.pkl')
                self.dp.myMS2.trans_ms2_file_linker_ion.append(file_name + '_trans_linker.pkl')
                self.dp.myMS2.trans_ms2_file_no_linker_ion.append(file_name + '_trans_no_linker.pkl')
            else:
                ms2_file_list = self.dp.myCFG.A2_PATH_MS2.split('|')
                for one_ms2_file in ms2_file_list:
                    if not os.path.exists(one_ms2_file):
                        Info = "Path " + one_ms2_file + ' is not exit !'
                        logGetWarning(Info)
                    self.dp.myMS2.ms2_num += 1
                    file_name = one_ms2_file.split('.')[0]
                    self.dp.myMS2.original_file_name.append(one_ms2_file)
                    self.dp.myMS2.ms2_file.append(file_name + '.mgf')
                    self.dp.myMS2.trans_ms2_file.append(file_name + '.pkl')
                    self.dp.myMS2.trans_ms2_file_linker_ion.append(file_name + '_trans_linker.pkl')
                    self.dp.myMS2.trans_ms2_file_no_linker_ion.append(file_name + '_trans_no_linker.pkl')

        elif self.dp.myCFG.A1_TYPE_MS2 in ['d']:

            if not os.path.isdir(self.dp.myCFG.A2_PATH_MS2):
                self.dp.myMS2.ms2_num = 1
                file_name = os.path.dirname(self.dp.myCFG.A2_PATH_MS2) + '\\' + self.dp.myCFG.A2_PATH_MS2.split('\\')[-1].split('.')[0]
                self.dp.myMS2.original_file_name.append(self.dp.myCFG.A2_PATH_MS2)
                self.dp.myMS2.ms2_file.append(file_name + '.mgf')
                self.dp.myMS2.trans_ms2_file.appene(file_name + '.pkl')
                self.dp.myMS2.trans_ms2_file_linker_ion.append(file_name + '_trans_linker.pkl')
                self.dp.myMS2.trans_ms2_file_no_linker_ion.append(file_name + '_trans_no_linker.pkl')
            else:
                ms2_file_list = self.dp.myCFG.A2_PATH_MS2.split('|')
                for one_ms2_file in ms2_file_list:
                    if not os.path.exists(one_ms2_file):
                        Info = "Path " + one_ms2_file + ' is not exit !'
                        logGetWarning(Info)
                        continue
                    self.dp.myMS2.ms2_num += 1
                    file_name = os.path.dirname(one_ms2_file) + '\\' + one_ms2_file.split('\\')[-1].split('.')[0]
                    self.dp.myMS2.original_file_name.append(one_ms2_file)
                    self.dp.myMS2.ms2_file.append(file_name + '.mgf')
                    self.dp.myMS2.trans_ms2_file.append(file_name + '.pkl')
                    self.dp.myMS2.trans_ms2_file_linker_ion.append(file_name + '_trans_linker.pkl')
                    self.dp.myMS2.trans_ms2_file_no_linker_ion.append(file_name + '_trans_no_linker.pkl')

        elif self.dp.myCFG.A1_TYPE_MS2 in ['raw']:

            if not os.path.isdir(self.dp.myCFG.A2_PATH_MS2) and "|" not in self.dp.myCFG.A2_PATH_MS2:
                self.dp.myMS2.ms2_num = 1
                file_name = self.dp.myCFG.A2_PATH_MS2.split('.')[0] + '_HCDFT'
                self.dp.myMS2.original_file_name.append(self.dp.myCFG.A2_PATH_MS2)
                self.dp.myMS2.ms2_file.append(file_name + '.mgf')
                self.dp.myMS2.trans_ms2_file.append(file_name + '.pkl')
                self.dp.myMS2.trans_ms2_file_linker_ion.append(file_name + '_trans_linker.pkl')
                self.dp.myMS2.trans_ms2_file_no_linker_ion.append(file_name + '_trans_no_linker.pkl')
            else:
                ms2_file_list = self.dp.myCFG.A2_PATH_MS2.split('|')
                for one_ms2_file in ms2_file_list:
                    if not os.path.exists(one_ms2_file):
                        Info = "Path " + one_ms2_file + ' is not exit !'
                        logGetWarning(Info)
                        continue
                    self.dp.myMS2.ms2_num += 1
                    file_name = one_ms2_file.split('.')[0] + '_HCDFT'
                    self.dp.myMS2.original_file_name.append(one_ms2_file)
                    self.dp.myMS2.ms2_file.append(file_name + '.mgf')
                    self.dp.myMS2.trans_ms2_file.append(file_name + '.pkl')
                    self.dp.myMS2.trans_ms2_file_linker_ion.append(file_name + '_trans_linker.pkl')
                    self.dp.myMS2.trans_ms2_file_no_linker_ion.append(file_name + '_trans_no_linker.pkl')

    def __captain_init_mod(self):

        # set fixed modifications and variable modifications
        if self.dp.myCFG.B4_NAME_MOD_FIX.strip() == "":
            fix_mod = []
        else:
            fix_mod = self.dp.myCFG.B4_NAME_MOD_FIX.split(';')
        if self.dp.myCFG.B5_NAME_MOD_VAR.strip() == "":
            var_mod = []
        else:
            var_mod = self.dp.myCFG.B5_NAME_MOD_VAR.split(';')

        for mod_name in fix_mod:
            if mod_name in self.dp.myMOD.fixmodname2index:
                return -1
            self.dp.myMOD.fixmodname2index[mod_name] = len(self.dp.myMOD.fix_mod_list)
            self.dp.myMOD.fix_mod_list.append(mod_name)
        for mod_name in var_mod:
            if mod_name in self.dp.myMOD.fixmodname2index:
                return -1
            self.dp.myMOD.varmodname2index[mod_name] = len(self.dp.myMOD.var_mod_list)
            self.dp.myMOD.var_mod_list.append(mod_name)
        # For self.fix_**: the key is the amino acid, value is the Mod object
        # For self.var_**: the key is the amino acid, value is the list of Mod objects
        self.__soldier_init_fix_mod(fix_mod)
        self.__soldier_init_var_mod(var_mod)
        # initial the modification index

        return 0

    def __soldier_init_fix_mod(self, fix_mod):

        for mod_name in fix_mod:
            if mod_name not in self.dp.myINI.DIC_MOD:
                print("[Error] {} not in modification.ini".format(mod_name))
                return -1
            mod_index = self.dp.myMOD.fixmodname2index[mod_name]
            mod = self.dp.myINI.DIC_MOD[mod_name]
            self.dp.myMOD.fix_mod_num += 1
            for a in mod.aa:
                if mod.type == 0:
                    if a not in self.dp.myMOD.fix_mod_dic.keys():
                        self.dp.myMOD.fix_mod_dic[a] = mod_index
                elif mod.type == - 1:
                    if a not in self.dp.myMOD.fix_mod_dic_pep_N.keys():
                        self.dp.myMOD.fix_mod_dic_pep_N[a] = mod_index
                elif mod.type == 1:
                    if a not in self.dp.myMOD.fix_mod_dic_pep_C.keys():
                        self.dp.myMOD.fix_mod_dic_pep_C[a].append(mod_index)
                elif mod.type == - 2:
                    if a not in self.dp.myMOD.fix_mod_dic_pro_N.keys():
                        self.dp.myMOD.fix_mod_dic_pro_N[a].append(mod_index)
                elif mod.type == 2:
                    if a not in self.dp.myMOD.fix_mod_dic_pro_C.keys():
                        self.dp.myMOD.fix_mod_dic_pro_C[a].append(mod_index)

    def __soldier_init_var_mod(self, var_mod):

        for mod_name in var_mod:
            if mod_name not in self.dp.myINI.DIC_MOD:
                print("[Error] {} not in modification.ini".format(mod_name))
                return -1
            mod_index = self.dp.myMOD.varmodname2index[mod_name]
            mod = self.dp.myINI.DIC_MOD[mod_name]
            self.dp.myMOD.var_mod_num += 1
            for a in mod.aa:
                if mod.type == 0:
                    if a not in self.dp.myMOD.var_mod_dic.keys():
                        self.dp.myMOD.var_mod_dic[a] = [mod_index]
                    else:
                        self.dp.myMOD.var_mod_dic[a].append(mod_index)
                elif mod.type == - 1:
                    if a not in self.dp.myMOD.var_mod_dic_pep_N.keys():
                        self.dp.myMOD.var_mod_dic_pep_N[a] = [mod_index]
                    else:
                        self.dp.myMOD.var_mod_dic_pep_N[a].append(mod_index)
                elif mod.type == 1:
                    if a not in self.dp.myMOD.var_mod_dic_pep_C.keys():
                        self.dp.myMOD.var_mod_dic_pep_C[a] = [mod_index]
                    else:
                        self.dp.myMOD.var_mod_dic_pep_C[a].append(mod_index)
                elif mod.type == - 2:
                    if a not in self.dp.myMOD.var_mod_dic_pro_N.keys():
                        self.dp.myMOD.var_mod_dic_pro_N[a] = [mod_index]
                    else:
                        self.dp.myMOD.var_mod_dic_pro_N[a].append(mod_index)
                elif mod.type == 2:
                    if a not in self.dp.myMOD.var_mod_dic_pro_C.keys():
                        self.dp.myMOD.var_mod_dic_pro_C[a] = [mod_index]
                    else:
                        self.dp.myMOD.var_mod_dic_pro_C[a].append(mod_index)

    def __soldier_init_mono_mod(self, linker_name, complex_search):

        mono_name = "mono_" + linker_name

        if self.dp.myINI.DIC_LINKER[linker_name].mono_mass == 0.0:
            pass
        else:
            # 加可变修饰处理
            if complex_search:
                self.dp.myINI.DIC_MOD[mono_name] = CModData(mono_name, self.dp.myINI.DIC_LINKER[linker_name].mono_mass, 3, '')
            else:
                if mono_name in self.dp.myINI.DIC_MOD.keys():
                    logGetError("[Error] Mono-link may same with one modification, name : " + mono_name)
                else:
                    self.dp.myINI.DIC_MOD[mono_name] = CModData(mono_name, self.dp.myINI.DIC_LINKER[linker_name].mono_mass, 3, '')

            self.dp.myMOD.var_mod_list.append(mono_name)
            self.dp.myMOD.varmodname2index[mono_name] = len(self.dp.myMOD.var_mod_list) - 1

            mod_index = self.dp.myMOD.varmodname2index[mono_name]

            for aa in self.dp.myINI.DIC_LINKER[linker_name].same_site.keys():
                if self.dp.myINI.DIC_LINKER[linker_name].same_site[aa] == 1:
                    if aa == -2:
                        if aa not in self.dp.myMOD.var_mod_dic_pro_N.keys():
                            self.dp.myMOD.var_mod_dic_pro_N[aa] = [mod_index]
                        else:
                            self.dp.myMOD.var_mod_dic_pro_N[aa].append(mod_index)
                    elif aa == -1:
                        if aa not in self.dp.myMOD.var_mod_dic_pep_N.keys():
                            self.dp.myMOD.var_mod_dic_pep_N[aa] = [mod_index]
                        else:
                            self.dp.myMOD.var_mod_dic_pep_N[aa].append(mod_index)
                    elif aa == 1:
                        if aa not in self.dp.myMOD.var_mod_dic_pep_C.keys():
                            self.dp.myMOD.var_mod_dic_pep_C[aa] = [mod_index]
                        else:
                            self.dp.myMOD.var_mod_dic_pep_C[aa].append(mod_index)
                    elif aa == 2:
                        if aa not in self.dp.myMOD.var_mod_dic_pro_C.keys():
                            self.dp.myMOD.var_mod_dic_pro_C[aa] = [mod_index]
                        else:
                            self.dp.myMOD.var_mod_dic_pro_C[aa].append(mod_index)
                    else:
                        if aa not in self.dp.myMOD.var_mod_dic.keys():
                            self.dp.myMOD.var_mod_dic[aa] = [mod_index]
                        else:
                            self.dp.myMOD.var_mod_dic[aa].append(mod_index)

            for aa in self.dp.myINI.DIC_LINKER[linker_name].alpha_site.keys():
                if self.dp.myINI.DIC_LINKER[linker_name].alpha_site[aa] == 1:
                    if aa == -2:
                        if aa not in self.dp.myMOD.var_mod_dic_pro_N.keys():
                            self.dp.myMOD.var_mod_dic_pro_N[aa] = [mod_index]
                        else:
                            self.dp.myMOD.var_mod_dic_pro_N[aa].append(mod_index)
                    elif aa == -1:
                        if aa not in self.dp.myMOD.var_mod_dic_pep_N.keys():
                            self.dp.myMOD.var_mod_dic_pep_N[aa] = [mod_index]
                        else:
                            self.dp.myMOD.var_mod_dic_pep_N[aa].append(mod_index)
                    elif aa == 1:
                        if aa not in self.dp.myMOD.var_mod_dic_pep_C.keys():
                            self.dp.myMOD.var_mod_dic_pep_C[aa] = [mod_index]
                        else:
                            self.dp.myMOD.var_mod_dic_pep_C[aa].append(mod_index)
                    elif aa == 2:
                        if aa not in self.dp.myMOD.var_mod_dic_pro_C.keys():
                            self.dp.myMOD.var_mod_dic_pro_C[aa] = [mod_index]
                        else:
                            self.dp.myMOD.var_mod_dic_pro_C[aa].append(mod_index)
                    else:
                        if aa not in self.dp.myMOD.var_mod_dic.keys():
                            self.dp.myMOD.var_mod_dic[aa] = [mod_index]
                        else:
                            self.dp.myMOD.var_mod_dic[aa].append(mod_index)

            for aa in self.dp.myINI.DIC_LINKER[linker_name].beta_site.keys():
                if self.dp.myINI.DIC_LINKER[linker_name].beta_site[aa] == 1:
                    if aa == -2:
                        if aa not in self.dp.myMOD.var_mod_dic_pro_N.keys():
                            self.dp.myMOD.var_mod_dic_pro_N[aa] = [mod_index]
                        else:
                            self.dp.myMOD.var_mod_dic_pro_N[aa].append(mod_index)
                    elif aa == -1:
                        if aa not in self.dp.myMOD.var_mod_dic_pep_N.keys():
                            self.dp.myMOD.var_mod_dic_pep_N[aa] = [mod_index]
                        else:
                            self.dp.myMOD.var_mod_dic_pep_N[aa].append(mod_index)
                    elif aa == 1:
                        if aa not in self.dp.myMOD.var_mod_dic_pep_C.keys():
                            self.dp.myMOD.var_mod_dic_pep_C[aa] = [mod_index]
                        else:
                            self.dp.myMOD.var_mod_dic_pep_C[aa].append(mod_index)
                    elif aa == 2:
                        if aa not in self.dp.myMOD.var_mod_dic_pro_C.keys():
                            self.dp.myMOD.var_mod_dic_pro_C[aa] = [mod_index]
                        else:
                            self.dp.myMOD.var_mod_dic_pro_C[aa].append(mod_index)
                    else:
                        if aa not in self.dp.myMOD.var_mod_dic.keys():
                            self.dp.myMOD.var_mod_dic[aa] = [mod_index]
                        else:
                            self.dp.myMOD.var_mod_dic[aa].append(mod_index)

            self.dp.myMOD.mono_mod_num += 1

class CFunctionEnzyme:

    def __init__(self, inputDP):
        self.dp = inputDP
        for flag in self.dp.myENZ.enzyme_flag:
            if flag != "N" and flag != "C":
                print("[Error] Enzyme split flag can only be N or C.")
                info = "[Error] Enzyme split flag can only be N or C."
                logGetError(info)
        if self.dp.myCFG.B2_TYPE_DIGEST < 0 or self.dp.myCFG.B2_TYPE_DIGEST > 2:
            print("[Error] Enzyme type can only be 0, 1, 2.")
            info = "[Error] Enzyme type can only be 0, 1, 2."
            logGetError(info)
        if self.dp.myCFG.B3_NUMBER_MAX_MISS_CLV < 0:
            print("[Error] Enzyme miss cleavage number must be >= 0.")
            info = "[Error] Enzyme miss cleavage number must be >= 0."
            logGetError(info)
        if self.dp.myCFG.D6_MASS_PEP_LOW < 0:
            print("[Error] Peptide mass minimum must be >= 0.")
            info = "[Error] Peptide mass minimum must be >= 0."
            logGetError(info)
        if self.dp.myCFG.D8_LEN_PEP_LOW < 0:
            print("[Error] Peptide length minimum must be >= 0.")
            info = "[Error] Peptide length minimum must be >= 0."
            logGetError(info)

    def split_protein(self, protein_seq, peptide_seq_list, peptide_pos_list, print_flag=False, pro_index=0):
        # split the protein sequence to list of CPeptide
        split_sites = []
        self.find_split_site(protein_seq, split_sites)
        if self.dp.myCFG.B2_TYPE_DIGEST == 0:
            self.__captain_specific_split(split_sites, protein_seq, peptide_seq_list, peptide_pos_list, pro_index)
        elif self.dp.myCFG.B2_TYPE_DIGEST == 1:
            self.__captain_semi_specific_split(split_sites, protein_seq, peptide_seq_list, peptide_pos_list)
        else:
            self.__captain_non_specific_split(split_sites, protein_seq, peptide_seq_list, peptide_pos_list)

        if print_flag:
            print("[Info]Protein sq", protein_seq)
            for p in peptide_seq_list:
                print("[Info]Peptide sq", p)

    def find_split_site(self, protein_seq, sites):
        site_dict = {}
        aa_split_type = {}
        for i, aa in enumerate(self.dp.myENZ.enzyme_aa):
            for one_aa in aa:
                flag = self.dp.myENZ.enzyme_flag[i]
                if flag == "C":
                    if one_aa not in aa_split_type:
                        aa_split_type[one_aa] = 0
                    aa_split_type[one_aa] += 1
                elif flag == "N":
                    if one_aa not in aa_split_type:
                        aa_split_type[one_aa] = 0
                    aa_split_type[one_aa] -= 1
        sites.append(-1)
        site_dict[-1] = 1
        for i in range(len(protein_seq)):
            p = protein_seq[i]
            if p not in aa_split_type:
                continue
            if aa_split_type[p] >= 0:
                if i != len(protein_seq) - 1 and i not in site_dict:
                    sites.append(i)
                    site_dict[i] = 1
            if aa_split_type[p] <= 0:
                if i != 0 and i - 1 not in site_dict:
                    sites.append(i - 1)
                    site_dict[i - 1] = 1

        if len(protein_seq) - 1 not in site_dict:
            sites.append(len(protein_seq) - 1)
            site_dict[len(protein_seq) - 1] = 1
        sites.sort()

    def __captain_specific_split(self, split_sites, protein_seq, peptide_seq_list, peptide_pos_list, pro_index=0):
        # specific cut
        # For example: protein sequence is "ACDEFKGGGRPSQED"
        # peptide_seq_list, peptide_pos_list: "ACDEFK", "GGGR", "PSQED"; [0,6],[6,10],[10,15]
        for i in range(len(split_sites) - 1):
            start_site = split_sites[i] + 1
            for j in range(1, self.dp.myCFG.B3_NUMBER_MAX_MISS_CLV + 2):
                if i + j >= len(split_sites):
                    break
                end_site = split_sites[i + j]
                if not self.__soldier_check_length(start_site, end_site + 1):
                    continue
                peptide_seq = protein_seq[start_site:end_site + 1]
                peptide_seq_list.append(peptide_seq)
                peptide_pos_list.append([start_site, end_site + 1])
        # split the protein when the first amino acid is "M"
        cur_len = len(peptide_pos_list)
        for i in range(cur_len):
            if peptide_pos_list[i][0] == 0 and peptide_seq_list[i].startswith("M"):
                end_site = peptide_pos_list[i][1]
                if not self.__soldier_check_length(1, end_site):
                    continue
                peptide_seq_list.append(peptide_seq_list[i][1:])
                peptide_pos_list.append([1, end_site])

    def __captain_semi_specific_split(self, split_sites, protein_seq, peptide_seq_list, peptide_pos_list):

        pass

    def __captain_non_specific_split(self, split_sites, protein_seq, peptide_seq_list, peptide_pos_list):

        pass

    def __soldier_check_length(self, start_site, end_site):
        # check the length of peptide sequence is in [min_len, max_len]
        len_peptide_seq = end_site - start_site
        return len_peptide_seq >= self.dp.myCFG.D8_LEN_PEP_LOW and len_peptide_seq <= self.dp.myCFG.D9_LEN_PEP_UP

class CFunctionSpectrum:

    def __init__(self, inputDP):

        self.dp = inputDP

    def trans_ms2(self):
        for n in range(self.dp.myMS2.ms2_num):
            self.__captain_trans_ms2TOpkl(self.dp.myMS2.ms2_file[n], self.dp.myMS2.trans_ms2_file[n], False)

    def __captain_trans_ms2TOpkl(self, one_mgf_file, one_trans_mgf_file, out_mgf=True):
        if not os.path.exists(one_trans_mgf_file):
            # init
            dataMS2 = CFileMS2()
            op_INIT_CFILE_MS2(dataMS2)
            # 原始谱图变成为pkl文件存储
            if one_mgf_file.split('.')[-1] in ['ms2', 'MS2']:
                self.__soldier_load_ms2(one_mgf_file, dataMS2)
            elif one_mgf_file.split('.')[-1] in ['mgf', 'MGF']:
                time_1 = time.time()
                self.__soldier_load_mgf(one_mgf_file, dataMS2)
                time_2 = time.time()
                print("Loading mgf time : %.4f"%(time_2 - time_1))
            else:
                info = "MS2 file may be wrong !"
                logGetError(info)
            self.__soldier_trans_ms2_mz(dataMS2)
            time_1 = time.time()
            self.__soldier_output_ms2_pkl(one_trans_mgf_file, dataMS2)
            time_2 = time.time()
            print("out original mgf pkl time : ", time_2 - time_1)
            if out_mgf:
                self.__soldier_output_ms2_mgf(one_mgf_file, dataMS2)
            del dataMS2

        else:

            pass

    def __soldier_load_ms2(self, path, spectrum_list):

        pass

    def __soldier_load_mgf(self, path, dataMS2):
        print('[Info] MSMS file path : ', path)
        peak_list = []
        with open(path, 'r') as f:
            first_time_scan = True
            buffer = f.readlines()
            len_buffer = len(buffer)
            line_num = 1
            spe_num = 0
            mobility = None
            print("[Info]Loading spectrum")
            for now_index, line in enumerate(buffer):
                line_num += 1
                line = line.strip()
                if line == "":
                    continue
                if line[0] == 'T':
                    title = line[6:]
                    scan = int(line.split(".")[-4])
                    # TITLE=xxx...x.SCAN.SCAN.CHARGE.RANK.dta
                    # 逆序主要是怕前面的文件名里有.字符
                    if scan in dataMS2.INDEX_SCAN:
                        # 已经进行了初始化，且记录了二级谱信息
                        dataMS2.MATRIX_FILE_NAME[scan].append(line[6:])
                        first_time_scan = False
                    else:
                        # 第一次记录信息
                        first_time_scan = True
                        dataMS2.INDEX_SCAN.append(scan)
                        dataMS2.LIST_PEAK_MOZ[scan] = []  # 初始化，不能省略
                        dataMS2.LIST_PEAK_INT[scan] = []

                        dataMS2.MATRIX_PEAK_MOZ_NO_LINKER[scan] = []  # 初始化，不能省略
                        dataMS2.MATRIX_PEAK_INT_NO_LINKER[scan] = []

                        dataMS2.MATRIX_PEAK_MOZ_LINKER[scan] = []  # 初始化，不能省略
                        dataMS2.MATRIX_PEAK_INT_LINKER[scan] = []

                        dataMS2.MATRIX_FILE_NAME[scan] = [title]
                        dataMS2.MATRIX_PRECURSOR_MOZ[scan] = []
                        dataMS2.MATRIX_PRECURSOR_CHARGE[scan] = []
                        dataMS2.MATRIX_MOBILITY[scan] = []

                elif line[0] == 'P':
                    precursor_moz = float(line.split('=')[1].split(' ')[0])

                    dataMS2.MATRIX_PEAK_MOZ_NO_LINKER[scan].append([])
                    dataMS2.MATRIX_PEAK_INT_NO_LINKER[scan].append([])

                    dataMS2.MATRIX_PEAK_MOZ_LINKER[scan].append([])
                    dataMS2.MATRIX_PEAK_INT_LINKER[scan].append([])

                elif line[0] == 'C':
                    if line[-1] == '+':
                        line = line[:-1]
                    try:
                        charge = int(line[7])
                    except:
                        charge = None
                elif line.startswith("MOBILITY"):
                    mobility = float(line.split('=')[1])
                elif line[0] == "R":
                    RT = float(line[12:])
                    dataMS2.LIST_RET_TIME[scan] = RT
                    # dataMS2.LIST_ION_INJECTION_TIME[tmpScan] = float(t)
                elif line[0] >= '1' and line[0] <= '9':
                    try:
                        ind = line.find('\t')
                        p_mz = float(line[:ind])
                    except:
                        ind = line.find(' ')
                        p_mz = float(line[:ind])
                    p_inten = float(line[ind + 1:])
                    peak_list.append((p_mz, p_inten))
                elif line[0] == 'E':
                    if charge is None:
                        continue
                    precursor_moz = precursor_moz * charge - (charge - 1) * self.dp.myINI.MASS_PROTON_MONO
                    dataMS2.MATRIX_PRECURSOR_CHARGE[scan].append(charge)
                    dataMS2.MATRIX_PRECURSOR_MOZ[scan].append(precursor_moz)
                    dataMS2.MATRIX_MOBILITY[scan].append(mobility)
                    mobility = None
                    if spe_num % 10000 == 0:
                        print('[Info]Loading spectrum num : ' + str(spe_num) + ' len_buffer : ' + str(len_buffer) + ' now :' + str(line_num))
                    spe_num += 1
                    precursor_mass = precursor_moz - self.dp.myINI.MASS_PROTON_MONO
                    if len(peak_list) == 0:
                        pass
                    else:

                        if len(peak_list) - 1 >= self.dp.myCFG.D3_NUMBER_SELECT_PEAK:
                            n = self.dp.myCFG.D3_NUMBER_SELECT_PEAK
                        else:
                            n = len(peak_list) - 1
                        # 按强度排序取top-n
                        peak_list.sort(key=lambda x: x[1], reverse=True)
                        coarse_peak_list = peak_list[:n + 1]
                        # 再按质荷比规整好
                        coarse_peak_list.sort(key=lambda x: x[0], reverse=False)
                        if first_time_scan:
                            # 如果是第一次扫描到这个谱图，那么添加在原始谱图中要做添加，否则跳过
                            # 记录完整的原始谱图,scan号一样母离子不一样也不用多记录，因为都是一个谱图
                            for one_mz_inten in coarse_peak_list:
                                dataMS2.LIST_PEAK_MOZ[scan].append(one_mz_inten[0])
                                dataMS2.LIST_PEAK_INT[scan].append(one_mz_inten[1])
                        no_liner_peak_list = []
                        for one_mz_inten in coarse_peak_list:
                            for charge_num in range(charge):
                                mz_charge = charge_num + 1
                                mz_mass = one_mz_inten[0] * mz_charge - mz_charge * self.dp.myINI.MASS_PROTON_MONO  # 这根峰不带电荷质量
                                int_mz_mass = int(round(mz_mass * self.dp.myCFG.D12_MULTI_MASS))
                                if int_mz_mass <= 0 or int_mz_mass > math.ceil((precursor_mass - self.dp.myCFG.D6_MASS_PEP_LOW) * self.dp.myCFG.D12_MULTI_MASS) or int_mz_mass > self.dp.myCFG.max_mz_index_size:
                                    pass
                                else:
                                    no_liner_peak_list.append((int_mz_mass, one_mz_inten[1]))
                        no_liner_peak_list.sort(key=lambda x: x[0], reverse=False)
                        for one_mz_inten in no_liner_peak_list:
                            dataMS2.MATRIX_PEAK_MOZ_NO_LINKER[scan][-1].append(one_mz_inten[0])
                            dataMS2.MATRIX_PEAK_INT_NO_LINKER[scan][-1].append(one_mz_inten[1])
                    peak_list = []
        dataMS2.MS2_NUM = spe_num
        print('\n[Info]Loading spectrum num : ' + str(dataMS2.MS2_NUM) + ' len_buffer : ' + str(len_buffer) + ' now :' + str(line_num))

    def __soldier_trans_ms2_mz(self, dataMS2):
        for spectrum_num, scan in enumerate(dataMS2.INDEX_SCAN):
            for one_precursor_moz_index, one_precursor_moz in enumerate(dataMS2.MATRIX_PRECURSOR_MOZ[scan]):
                precursor_mass = one_precursor_moz - self.dp.myINI.MASS_PROTON_MONO
                linker_peak_list = []
                for one_mz_index, one_mz in enumerate(dataMS2.LIST_PEAK_MOZ[scan]):
                    one_mz_inten = dataMS2.LIST_PEAK_INT[scan][one_mz_index]
                    for charge_num in range(self.dp.myCFG.max_fragment_charge):
                        # 假设这根峰是跨linker的b离子，把互补的y离子算出来即y离子不会跨linker
                        # 现在都假设只带单点荷和两电荷
                        # 不带电荷离子质量
                        mz_charge = charge_num + 1
                        mz_mass = one_mz * mz_charge - mz_charge * self.dp.myINI.MASS_PROTON_MONO  # 这根峰不带电荷质量
                        mz_complementation_mass = precursor_mass - mz_mass  # 这根峰假设为跨linker的峰不带电荷质量
                        int_mz_mass = int(round(mz_complementation_mass * self.dp.myCFG.D12_MULTI_MASS))
                        # 因为都换算成不带电质量，所以子离子超出母离子的部分直接就不算了
                        if mz_mass <= self.dp.myCFG.D6_MASS_PEP_LOW or int_mz_mass > math.ceil((precursor_mass - self.dp.myCFG.D6_MASS_PEP_LOW) * self.dp.myCFG.D12_MULTI_MASS) or int_mz_mass < 0 or int_mz_mass > self.dp.myCFG.max_mz_index_size:
                            pass
                        else:
                            # 不把质量相同的峰强度归到一个里面
                            linker_peak_list.append((int_mz_mass, one_mz_inten))
                linker_peak_list.sort(key=lambda x: x[0], reverse=False)
                one_spe_linker_mz = []
                one_spe_linker_inten = []
                for one_mz_inten in linker_peak_list:
                    one_spe_linker_mz.append(one_mz_inten[0])
                    one_spe_linker_inten.append(one_mz_inten[1])
                dataMS2.MATRIX_PEAK_MOZ_LINKER[scan][one_precursor_moz_index] = one_spe_linker_mz
                dataMS2.MATRIX_PEAK_INT_LINKER[scan][one_precursor_moz_index] = one_spe_linker_inten

    def __soldier_output_ms2_pkl(self, path, dataMS2):

        functionPickle = CFunctionPickle()
        functionPickle.write_data_to_pkl(dataMS2, path)

    def __soldier_output_ms2_mgf(self, path, dataMS2):

        with open(path.split('.')[0] + '_trans_no_linker.mgf', 'w') as f:
            for scan in dataMS2.INDEX_SCAN:
                for index, title in enumerate(dataMS2.MATRIX_FILE_NAME[scan]):
                    charge = dataMS2.MATRIX_PRECURSOR_CHARGE[scan][index]
                    precursor_moz = (dataMS2.MATRIX_PRECURSOR_MOZ[scan][index] + (charge - 1) * self.dp.myINI.MASS_PROTON_MONO) / charge
                    list_moz = dataMS2.MATRIX_PEAK_MOZ_NO_LINKER[scan][index]
                    list_inten = dataMS2.MATRIX_PEAK_INT_NO_LINKER[scan][index]
                    out_str = ''
                    f.write('BEGIN IONS\nTITLE=')
                    f.write(title)
                    f.write('\nCHARGE=')
                    f.write(str(charge) + '+')
                    f.write('\nPEPMASS=')
                    f.write('%.6f' %(precursor_moz))
                    f.write('\n')
                    for mz_index in range(len(list_moz)):
                        out_str = out_str + '%.0f' % (list_moz[mz_index]) + ' %.2f\n' % (list_inten[mz_index])
                    f.write(out_str)
                    f.write('END IONS\n')
        with open(path.split('.')[0] + '_trans_linker.mgf', 'w') as f:
            for scan in dataMS2.INDEX_SCAN:
                for index, title in enumerate(dataMS2.MATRIX_FILE_NAME[scan]):
                    charge = dataMS2.MATRIX_PRECURSOR_CHARGE[scan][index]
                    precursor_moz = (dataMS2.MATRIX_PRECURSOR_MOZ[scan][index] + (charge - 1) * self.dp.myINI.MASS_PROTON_MONO) / charge
                    list_moz = dataMS2.MATRIX_PEAK_MOZ_LINKER[scan][index]
                    list_inten = dataMS2.MATRIX_PEAK_INT_LINKER[scan][index]
                    out_str = ''
                    f.write('BEGIN IONS\nTITLE=')
                    f.write(title)
                    f.write('\nCHARGE=')
                    f.write(str(charge) + '+')
                    f.write('\nPEPMASS=')
                    f.write('%.6f' %(precursor_moz))
                    f.write('\n')
                    for mz_index in range(len(list_moz)):
                        out_str = out_str + '%.0f' % (list_moz[mz_index]) + ' %.2f\n' % (list_inten[mz_index])
                    f.write(out_str)
                    f.write('END IONS\n')

class CFunctionCheck:

    def __init__(self, inputDP):

        self.dp = inputDP

    def check_file(self):

        self.__check_open_result_file()

    def __check_open_result_file(self):

        #  输出文件有时候用excel打开了，需要检查下文件是否正常
        for i in range(len(FILE_NAME_PSM_RESULT)):
            for nameFile_key in FILE_NAME_PSM_RESULT[i].keys():
                path_1 = self.dp.myCFG.A5_PATH_RESULT_EXPORT + '\\' + FILE_NAME_PSM_RESULT[i][nameFile_key] + '.txt'
                path_2 = self.dp.myCFG.A5_PATH_RESULT_EXPORT + '\\res_' + FILE_NAME_PSM_RESULT[i][nameFile_key] + '.txt'
                self.__captain_try_open(path_1)
                self.__captain_try_open(path_2)

    def __captain_try_open(self, path):

        if os.path.exists(path):

            try:
                os.access(path, os.W_OK)

            except IOError:

                logGetError(path + ' is opened! Please close it and run the program again!')

class CFunctionFDR():

    def __init__(self, inputDP):
        self.dp = inputDP

    def calculate_cross_result_FDR(self, psm_result):
        TT = 0
        TD = 0
        DD = 0
        last_q_value = 0

        inter_TT = 0
        inter_TD = 0
        inter_DD = 0
        inter_last_q_value = 0

        intra_TT = 0
        intra_TD = 0
        intra_DD = 0
        intra_last_q_value = 0

        for index_all_match_result, one_match_result in enumerate(psm_result):
            if one_match_result.peptide_score == 0:
                break
            else:
                if one_match_result.protein_type_math == 1:
                    TT += 1
                elif one_match_result.protein_type_math == 0:
                    TD += 1
                elif one_match_result.protein_type_math == -1:
                    DD += 1

                if TT == 0:
                    now_q_value = 1
                else:
                    now_q_value = (TD - DD) / TT
                if last_q_value >= now_q_value:
                    pass
                else:
                    last_q_value = now_q_value
                if last_q_value >= 1:
                    last_q_value = 1
                one_match_result.FDR = last_q_value

                if one_match_result.inter_protein:
                    if one_match_result.protein_type_math == 1:
                        inter_TT += 1
                    elif one_match_result.protein_type_math == 0:
                        inter_TD += 1
                    elif one_match_result.protein_type_math == -1:
                        inter_DD += 1

                    if inter_TT == 0:
                        inter_now_q_value = 1
                    else:
                        inter_now_q_value = (inter_TD - inter_DD) / inter_TT
                    if inter_last_q_value >= inter_now_q_value:
                        pass
                    else:
                        inter_last_q_value = inter_now_q_value
                    if inter_last_q_value >= 1:
                        inter_last_q_value = 1
                    one_match_result.inter_FDR = inter_last_q_value

                else:
                    if one_match_result.protein_type_math == 1:
                        intra_TT += 1
                    elif one_match_result.protein_type_math == 0:
                        intra_TD += 1
                    elif one_match_result.protein_type_math == -1:
                        intra_DD += 1

                    if intra_TT == 0:
                        intra_now_q_value = 1
                    else:
                        intra_now_q_value = (intra_TD - intra_DD) / intra_TT
                    if intra_last_q_value >= intra_now_q_value:
                        pass
                    else:
                        intra_last_q_value = intra_now_q_value
                    if intra_last_q_value >= 1:
                        intra_last_q_value = 1
                    one_match_result.intra_FDR = intra_last_q_value

    def calculate_single_result_FDR(self, psm_result):
        T = 0
        D = 0
        last_q_value = 0
        for index_all_match_result, one_match_result in enumerate(psm_result):
            if one_match_result.peptide_score == 0:
                break
            else:
                if one_match_result.protein_type_math == 1:
                    T += 1
                elif one_match_result.protein_type_math == -1:
                    D += 1
            if T == 0:
                now_q_value = 1
            else:
                now_q_value = D / T
            if last_q_value >= now_q_value:
                FDR = last_q_value
            else:
                last_q_value = now_q_value
                FDR = now_q_value
            one_match_result.FDR = FDR

    def result_clean_scan(self, all_match_result):

        scan_record = [[] for i in range(VALUE_MAX_SCAN)]
        result_index = 0
        for one_match_result in all_match_result:
            if one_match_result.spectrum_title in scan_record[one_match_result.scan]:
                del all_match_result[result_index]
            else:
                scan_record[one_match_result.scan].append(one_match_result.spectrum_title)
            result_index += 1

    def result_filer_continue(self, all_match_result):
        out = []
        result_index = 0
        for one_match_result in all_match_result:
            if one_match_result.feature.continue_alpha_score > 1 and one_match_result.feature.continue_beta_score > 0:
                out.append(one_match_result)
            else:
                pass
                # del all_match_result[result_index]
            result_index += 1
        return out

class CFunctionExtractMS2():

    def __init__(self, inputDP):

        self.dp = inputDP

    def extract_raw(self, raw_index):

        path_raw = self.dp.myMS2.original_file_name[raw_index]
        print("[Info] Start pParse " + path_raw)
        cmd_pParse = PPARSE_EXE + " -D " + path_raw + " -t -0.5 -z 5 -i 1 -C 1 -a 0 -p 1 -S 1"
        os.system(cmd_pParse)
        print("[Info] Finish extract mgf file : " + path_raw)

    def extract_tims(self, tims_index):

        example_config_path = os.getcwd() + '\\' + TimsTOFReadConfigName
        config = CTimsTOFReaderConfig()
        self.timsTOFReader_file2config(example_config_path, config)
        path_d_file = self.dp.myMS2.original_file_name[tims_index]
        print("[Info] Start timsTOFReader " + path_d_file)
        folder_mgf_out = os.path.dirname(path_d_file)
        ms_name = self.dp.myMS2.original_file_name[tims_index].split('\\')[-1].split('.')[0]
        op_modify_timsReaderconfig(config, path_d_file, folder_mgf_out + "\\", ms_name)
        config_out_path = folder_mgf_out + '\\' + TimsTOFReadConfigName
        self.timsTOFReader_config2file(config_out_path, config)
        cmd_timsTOFReader = TIMSTOFREADER_EXE + " " + config_out_path
        os.system(cmd_timsTOFReader)
        print("[Info] Finish extract mgf file : " + path_d_file)

    def timsTOFReader_config2file(self, path, config):

        with open(path, 'w') as f:
            f.write("#[data]\n")
            f.write("FOLDER_TIMS_DATA=" + config.A1_FOLDER_TIMS_DATA + "\n")
            f.write("FOLDER_EXPORT=" + config.A2_FOLDER_EXPORT + "\n")
            f.write("MS_NAME=" + config.A3_MS_NAME + "\n")
            f.write("\n")
            f.write("#[performance]\n")
            f.write("EXTRACT_MS1=" + str(config.B1_EXTRACT_MS1) + "\n")
            f.write("EXTRACT_MS2=" + str(config.B2_EXTRACT_MS2) + "\n")
            f.write("EXTRACT_MGF=" + str(config.B3_EXTRACT_MGF) + "\n")
            f.write("MERGE_MS2=" + str(config.B4_MERGE_MS2) + "\n")
            f.write("MS1_INTENSITY_THRESHOLD=" + str(config.B5_MS1_INTENSITY_THRESHOLD) + "\n")
            f.write("MS2_INTENSITY_THRESHOLD=" + str(config.B6_MS2_INTENSITY_THRESHOLD) + "\n")
            f.write("MS1_FILTER=" + str(config.B7_MS1_FILTER) + "\n")
            f.write("MS2_FILTER=" + str(config.B8_MS2_FILTER) + "\n")
            f.write("\n")
            f.write("#[task]\n")
            f.write("TYPE_FLOW=" + str(config.C1_TYPE_FLOW) + "\n")
            f.write("\n")

    def timsTOFReader_file2config(self, file_path, config):

        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.find('=') < 0:
                    continue
                if line.find('#') >= 0:
                    line = line[:line.find('#')].strip()
                tmp = line.split('=')
                if len(tmp) < 2:
                    continue
                name = tmp[0]
                value = tmp[1]
                print(name, ":", value)
                try:
                    # [data]
                    if name == "FOLDER_TIMS_DATA":
                        config.A1_FOLDER_TIMS_DATA = value
                    elif name == "FOLDER_EXPORT":
                        config.A2_FOLDER_EXPORT = value
                        if not os.path.exists(config.A2_FOLDER_EXPORT):
                            os.makedirs(config.A2_FOLDER_EXPORT)
                    elif name == "MS_NAME":
                        config.A3_MS_NAME = value

                    # [performance]
                    elif name == "EXTRACT_MS1":
                        if int(value) == 1:
                            config.B1_EXTRACT_MS1 = 1
                        else:
                            config.B1_EXTRACT_MS1 = 0
                    elif name == "EXTRACT_MS2":
                        if int(value) == 1:
                            config.B2_EXTRACT_MS2 = 1
                        else:
                            config.B2_EXTRACT_MS2 = 0
                    elif name == "EXTRACT_MGF":
                        if int(value) == 1:
                            config.B3_EXTRACT_MGF = 1
                        else:
                            config.B3_EXTRACT_MGF = 0
                    elif name == "MERGE_MS2":
                        if int(value) == 1:
                            config.B4_MERGE_MS2 = 1
                        else:
                            config.B4_MERGE_MS2 = 0

                    elif name == "MS1_INTENSITY_THRESHOLD":
                        config.B5_MS1_INTENSITY_THRESHOLD = float(value)
                    elif name == "MS2_INTENSITY_THRESHOLD":
                        config.B6_MS2_INTENSITY_THRESHOLD = float(value)

                    elif name == "MS1_FILTER":
                        config.B7_MS1_FILTER = int(value)
                    elif name == "MS2_FILTER":
                        config.B8_MS2_FILTER = int(value)

                    # [task]
                    elif name == "TYPE_FLOW":
                        config.C1_TYPE_FLOW = int(value)

                except Exception as err:
                    print("[Error] occured in configure.txt", err)
