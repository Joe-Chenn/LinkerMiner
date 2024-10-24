# -*- mode: python ; coding: utf-8 -*-..
from MSDataResult import CRerankOnePeptideFeature, CWriteResultSinglePeptide
from MSData import CModSite, CPeak, CSpectrum
from MSDataResult import CSinglePeptideCoarseData, CSinglePeptideResult
from MSDataC import ION_INDEX_MZ_LOCK, ION_INDEX_PEP_LOCK
from MSSysterm import PRO_PKL_FILE, PEP_PRO_INDEX_FILE
from MSSysterm import DLL_SINGLE
from MSSysterm import FOLDER_INDEX, FILE_NAME_PSM_RESULT, FILE_NAME_LINK_SITE, FILE_NAME_ION_INDEX, FILE_NAME_PEP_DATA, FILE_NAME_ION_NUM, FILE_NAME_PEP_MASS, FILE_NAME_PEP_PKL
from MSFunction import CFunctionPickle, CFunctionFDR
from MSFunctionRerank import CFunctionOutOnePeptidePSM, CFunctionOutTwoPeptidePSM, CFunctionReadReasult
from MSOperator import op_create_peaks_index, op_fill_CSinglePeptideResult, op_fill_CWriteSinglePeptdidResult, op_fill_CSingleCoarseData
from MSOperatorScore import op_get_continue_score, op_get_peptide_by_continue_data, op_get_match_ion_score, op_get_peptide_continue_score
from MSTool_code import tool_new_get_ion_code, tool_new_get_ion_data
from MSTool import tool_create_aa_list, tool_generate_singlepepdite_by_ion, tool_generate_loop_type_by_ion, tool_binary_search_index, tool_set_two_list
from MSLogging import logGetError

from ctypes import *
import os
import math
import time
import operator
import numpy as np

class CFunctionOnePeptideSearch:

    link_type = 1  # 1单肽，2是loop-link

    def __init__(self, inputDP):

        self.dp = inputDP

    def update_class_variable(self, link_type):

        self.link_type = link_type

    def search(self, protein_list, id_psm_name=None, id_psm_scan=None):

        # 先把粗打分的参数初始化好
        functionSearchCoarse = CFunctionOnePeptideSearchCoarse(self.dp)
        pep_num = functionSearchCoarse.update_class_variable(self.link_type)
        # 细打分的参数初始化
        functionSearchFine = CFunctionOnePeptideSearchFine(self.dp)
        functionSearchFine.update_class_variable(self.link_type)

        functionPickle = CFunctionPickle()
        all_match_result = []
        # 一次粗打分和细打分调用时给出一个mgf的检索结果
        for mgf_index in range(self.dp.myMS2.ms2_num):
            dataMS2 = functionPickle.load_pkl_to_data(self.dp.myMS2.trans_ms2_file[mgf_index])
            # 把谱图弄成列表形式
            list_coarse_data = [None for i in range(dataMS2.MS2_NUM)]
            list_fine_data = [None for i in range(dataMS2.MS2_NUM)]
            self.trans_dataMS2_to_CSpectrum(dataMS2, list_coarse_data, list_fine_data, id_psm_name, id_psm_scan)
            print("[Info]#Spectra in {0} is {1}".format(self.dp.myMS2.trans_ms2_file[mgf_index], dataMS2.MS2_NUM))
            time_1 = time.time()
            coarse_peptide = [[] for i in range(dataMS2.MS2_NUM)]
            list_useful_result = [-1 for i in range(pep_num)]
            del dataMS2
            # 粗打分
            functionSearchCoarse.search(list_coarse_data, coarse_peptide, list_useful_result)
            time_2 = time.time()
            # 细打分
            s2pep_res = []
            functionSearchFine.search(protein_list, list_fine_data, coarse_peptide, list_useful_result, s2pep_res)
            all_match_result += s2pep_res
            time_3 = time.time()
            print("[Info] Coarse search using time : %.4f s ."%(time_2 - time_1))
            print("[Info] Fine search using time : %.4f s ."%(time_3 - time_2))
        all_match_result.sort(key=lambda all_match_result: all_match_result.peptide_score, reverse=True)
        # 写出结果
        if self.link_type in [1]:
            functionFDR = CFunctionFDR(self.dp)
            functionFDR.calculate_single_result_FDR(all_match_result)

            functionOutOnePeptidePSM = CFunctionOutOnePeptidePSM(self.dp)

            write_path = self.dp.myCFG.A5_PATH_RESULT_EXPORT
            write_pkl_path = self.dp.myCFG.A5_PATH_RESULT_EXPORT + "\\tmp\\"
            if not os.path.exists(write_pkl_path):
                os.makedirs(write_pkl_path)
            if not os.path.exists(write_path):
                os.makedirs(write_path)
            functionOutOnePeptidePSM.outputSinglePeptidePSM(all_match_result, write_pkl_path, write_path, False)
            functionOutOnePeptidePSM.update_class_variable(1, 0, 0)

        elif self.link_type in [2]:

            functionFDR = CFunctionFDR(self.dp)
            functionFDR.calculate_single_result_FDR(all_match_result)

            functionOutOnePeptidePSM = CFunctionOutOnePeptidePSM(self.dp)

            write_path = self.dp.myCFG.A5_PATH_RESULT_EXPORT
            write_pkl_path = self.dp.myCFG.A5_PATH_RESULT_EXPORT + "\\tmp\\"

            if not os.path.exists(write_pkl_path):
                os.makedirs(write_pkl_path)
            if not os.path.exists(write_path):
                os.makedirs(write_path)

            functionOutOnePeptidePSM.update_class_variable(1, 0, 0)
            functionOutOnePeptidePSM.outputLooplinkPSM(all_match_result, write_pkl_path, write_path)

    def trans_dataMS2_to_CSpectrum(self, dataMS2, list_coarse_data, list_fine_data, id_psm_name=None, id_psm_scan=None):
        '''
        :param dataMS2:
        :param split_coarse_MS2_data:
        :param id_psm_name: 这个是raw文件名的list，[raw1, raw2, raw3...]
        :param id_psm_scan: 这个是上面对应索引的位置按scan号大小存储，存储的spectrum tilte,
                            [
                             [
                              [], [raw1_spectrum.1.1.2.0.dta, raw1_spectrum.1.1.2.1.dta], [raw1_spectrum.2.2.3.0.dta], []...
                             ],
                             [
                              [], [], [raw2_spectrum.2.2.2.0.dta], [raw2_spectrum.3.3.3.0.dta, raw1_spectrum.3.3.3.1.dta]...
                             ],
                             [
                              [], [raw3_spectrum.1.1.1.0.dta], [raw3_spectrum.2.2.3.0.dta], [raw3_spectrum.3.3.3.0.dta]...
                             ]
                             ...
                            ]
        :return:
        '''
        # 需要把已有结果的谱图置空
        cur_dataMS2_name = dataMS2.MATRIX_FILE_NAME[dataMS2.INDEX_SCAN[0]][0].split('.')[0]
        if id_psm_name is None:
            id_index = None
            id_scan_list = []  # 这个是记录这个mgf文件的哪些scan号谱图被用过
        else:
            if cur_dataMS2_name in id_psm_name:
                id_index = id_psm_name.index(cur_dataMS2_name)
                id_scan_list = id_psm_scan[id_index]
            else:
                id_index = None
                id_scan_list = []

        spectrum_num = 0
        for ms2_scan_index, ms2_scan in enumerate(dataMS2.INDEX_SCAN):
            # 判断当前scan号的谱图是否已有结果
            if id_index is None:
                if len(dataMS2.LIST_PEAK_MOZ[ms2_scan]) < 1:
                    continue
                spe_max_int, spe_all_int, peaks_list, peaks_index = self.__captain_get_spectrum_data(dataMS2, ms2_scan)

                for precursor_moz_index, precursor_moz in enumerate(dataMS2.MATRIX_PRECURSOR_MOZ[ms2_scan]):

                    list_coarse_data[spectrum_num] = [
                                                    dataMS2.MATRIX_PRECURSOR_MOZ[ms2_scan][precursor_moz_index],
                                                    dataMS2.MATRIX_PRECURSOR_CHARGE[ms2_scan][precursor_moz_index],
                                                    dataMS2.MATRIX_PEAK_MOZ_NO_LINKER[ms2_scan][precursor_moz_index],
                                                    dataMS2.MATRIX_PEAK_INT_NO_LINKER[ms2_scan][precursor_moz_index],
                                                    dataMS2.MATRIX_FILE_NAME[ms2_scan][precursor_moz_index]
                                                   ]

                    list_fine_data[spectrum_num] = CSpectrum(
                                                             dataMS2.MATRIX_FILE_NAME[ms2_scan][precursor_moz_index],
                                                             precursor_moz,
                                                             dataMS2.MATRIX_PRECURSOR_CHARGE[ms2_scan][precursor_moz_index],
                                                             spe_max_int,
                                                             spe_all_int,
                                                             peaks_list,
                                                             peaks_index,
                                                             dataMS2.MATRIX_MOBILITY[ms2_scan][precursor_moz_index]
                                                            )
                    spectrum_num += 1
            else:
                if len(dataMS2.LIST_PEAK_MOZ[ms2_scan]) < 1:
                    continue
                spe_max_int, spe_all_int, peaks_list, peaks_index = self.__captain_get_spectrum_data(dataMS2, ms2_scan)
                for precursor_moz_index, precursor_moz in enumerate(dataMS2.MATRIX_PRECURSOR_MOZ[ms2_scan]):
                    if len(id_scan_list[ms2_scan]) < 1:
                        list_coarse_data[spectrum_num] = [
                                                        dataMS2.MATRIX_PRECURSOR_MOZ[ms2_scan][precursor_moz_index],
                                                        dataMS2.MATRIX_PRECURSOR_CHARGE[ms2_scan][precursor_moz_index],
                                                        dataMS2.MATRIX_PEAK_MOZ_NO_LINKER[ms2_scan][precursor_moz_index],
                                                        dataMS2.MATRIX_PEAK_INT_NO_LINKER[ms2_scan][precursor_moz_index],
                                                        dataMS2.MATRIX_FILE_NAME[ms2_scan][precursor_moz_index]
                                                       ]
                        list_fine_data[spectrum_num] = CSpectrum(
                                                                 dataMS2.MATRIX_FILE_NAME[ms2_scan][precursor_moz_index],
                                                                 precursor_moz,
                                                                 dataMS2.MATRIX_PRECURSOR_CHARGE[ms2_scan][precursor_moz_index],
                                                                 spe_max_int,
                                                                 spe_all_int,
                                                                 peaks_list,
                                                                 peaks_index,
                                                                 dataMS2.MATRIX_MOBILITY[ms2_scan][precursor_moz_index]
                                                                )
                    else:
                        if dataMS2.MATRIX_FILE_NAME[ms2_scan][precursor_moz_index] in id_scan_list[ms2_scan]:
                            pass
                        else:
                            list_coarse_data[spectrum_num] = [
                                                            dataMS2.MATRIX_PRECURSOR_MOZ[ms2_scan][precursor_moz_index],
                                                            dataMS2.MATRIX_PRECURSOR_CHARGE[ms2_scan][precursor_moz_index],
                                                            dataMS2.MATRIX_PEAK_MOZ_NO_LINKER[ms2_scan][precursor_moz_index],
                                                            dataMS2.MATRIX_PEAK_INT_NO_LINKER[ms2_scan][precursor_moz_index],
                                                            dataMS2.MATRIX_FILE_NAME[ms2_scan][precursor_moz_index]
                                                           ]
                            list_fine_data[spectrum_num] = CSpectrum(
                                                                     dataMS2.MATRIX_FILE_NAME[ms2_scan][precursor_moz_index],
                                                                     precursor_moz,
                                                                     dataMS2.MATRIX_PRECURSOR_CHARGE[ms2_scan][precursor_moz_index],
                                                                     spe_max_int,
                                                                     spe_all_int,
                                                                     peaks_list,
                                                                     peaks_index,
                                                                     dataMS2.MATRIX_MOBILITY[ms2_scan][precursor_moz_index]
                                                                    )
                    spectrum_num += 1

    def __captain_get_spectrum_data(self, dataMS2, scan):
        peak_list = []
        mz_list = dataMS2.LIST_PEAK_MOZ[scan]
        inten_list = dataMS2.LIST_PEAK_INT[scan]
        spe_all_int = 0.0
        for one_inten_index, one_inten in enumerate(inten_list):
            spe_all_int += one_inten
            peak_list.append(CPeak(mz_list[one_inten_index], one_inten))
        spe_max_int = max(inten_list)
        peaks_index = op_create_peaks_index(peak_list, 1)  # 这里是把谱图的碎片离子按质量大小排好，1mz作为一个区间
        return spe_max_int, spe_all_int, peak_list, peaks_index

class CFunctionOnePeptideSearchCoarse:

    # 类变量，即使多个它们只存一份就行，谁在使用之前就做一次更新
    mz_index_data = []  # 这个是离子索引，mz-> [肽段编码]
    list_pep_data_pkl = []  # 这个是记录满足交联条件的肽段索引，index->肽段编码
    list_ion_num = []  # 这个是记录肽段的碎片离子数目，index->离子数
    list_pep_mass = []  # 这个是记录肽段质量的，不包括交联剂质量
    list_link_site = []

    list_code_split_pkl_file = []

    link_type = 1  # 1是单肽，2是loop-link

    def __init__(self, inputDP):

        self.dp = inputDP

    def update_class_variable(self, link_type):
        # 调用下面的必须先弄这个

        self.link_type = link_type
        self.__captain_load_index()
        self.__captain_get_pkl_file_split_code()
        pep_num = len(self.list_pep_data_pkl)
        return pep_num

    def __captain_load_index(self):

        functionPickle = CFunctionPickle()

        path = self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[1][self.link_type] + '\\' + FILE_NAME_ION_NUM[1][self.link_type] + '.pkl'
        self.list_ion_num = functionPickle.load_pkl_to_data(path)

        path = self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[1][self.link_type] + '\\' + FILE_NAME_PEP_MASS[1][self.link_type] + '.pkl'
        self.list_pep_mass = functionPickle.load_pkl_to_data(path)

        path = self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[1][self.link_type] + '\\' + FILE_NAME_PEP_PKL[1][self.link_type] + '.pkl'
        self.list_pep_data_pkl = functionPickle.load_pkl_to_data(path)

        if self.link_type != 1:
            path = self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[1][self.link_type] + '\\' + FILE_NAME_LINK_SITE[1][self.link_type] + '.pkl'
            self.list_link_site = functionPickle.load_pkl_to_data(path)

    def __captain_get_pkl_file_split_code(self):

        for i in range(int((self.dp.myCFG.D7_MASS_PEP_UP - self.dp.myCFG.D6_MASS_PEP_LOW) / self.dp.myCFG.D10_INDEX_SPLIT_MASS)):
            my_code = tool_new_get_ion_code(i + 1, 0, 0, 0)
            self.list_code_split_pkl_file.append(my_code)

    def import_DLL(self, bird):

        # 这里是做函数的输入输出类型的定义
        bird.malloc_ion_index.argtype = (c_int)
        bird.malloc_ion_index.restype = (POINTER(ION_INDEX_MZ_LOCK))

        bird.free_ion_index.argtypes = (c_int, POINTER(ION_INDEX_MZ_LOCK))
        bird.free_ion_index.restype = (c_int)

        bird.malloc_record_score_list.argtype = (c_long, c_int)
        bird.malloc_record_score_list.restype = (POINTER(c_double))

        bird.free_record_score_list.argtype = (POINTER(c_double))
        bird.free_record_score_list.restype = (c_int)

        bird.get_ion_index_match_num.argtype = (POINTER(c_int), c_int, POINTER(ION_INDEX_MZ_LOCK), c_int, c_long, c_long, c_int, c_double, POINTER(c_double), c_int)
        bird.get_ion_index_match_num.restype = (c_int)

        bird.get_ion_index_match_num_percent.argtype = (POINTER(c_long), POINTER(c_double), c_int, c_int, c_int, c_long, c_long)
        bird.get_ion_index_match_num_percent.restype = (c_int)

        bird.get_ion_index_match_score.argtype = (POINTER(c_int), c_int, POINTER(c_double), POINTER(ION_INDEX_MZ_LOCK), c_int, c_long, c_long, c_int, c_double, POINTER(c_double), c_int, c_int)
        bird.get_ion_index_match_score.restype = (c_int)

        bird.get_sum_score.argtype = (POINTER(c_double), c_long, c_int, c_int, c_int)
        bird.get_sum_score.restype = (c_int)

        bird.read_ion_index_bin.argtype = (POINTER(c_char_p), POINTER(ION_INDEX_MZ_LOCK), c_int)
        bird.read_ion_index_bin.restype = (c_int)

        # 这里是做函数的输入输出类型的定义

    def search(self, list_coarse_data, coarse_peptide, list_useful_result):

        if self.link_type == 1:
            # 单肽
            self.search_singlepeptide(list_coarse_data, coarse_peptide, list_useful_result)

        elif self.link_type == 2:
            # loop-link
            self.search_loop_link(list_coarse_data, coarse_peptide, list_useful_result)

        else:

            logGetError("[Error] Some unknown error in one peptide searching !")

    def search_singlepeptide(self, list_coarse_data, coarse_peptide, list_useful_result):
        # 记粗打分结果，是编码的结果
        self.link_type = 1
        self.__captain_coarse_single_peptide_search(list_coarse_data, coarse_peptide, list_useful_result)
        # 要回详细结果
        self.__captain_get_single_peptide_data(list_useful_result)

    def __captain_coarse_single_peptide_search(self, list_coarse_data, coarse_peptide, list_useful_peptide):

        # 提前处理的数据
        if self.dp.myCFG.is_ppm_pre:
            self.is_ppm_pre = 1
            self.fragment_ppm = c_double(self.dp.myCFG.C2_PPM_TOL_FRAGMENT)
            self.precursor_ppm = c_double(self.dp.myCFG.C1_PPM_TOL_PRECURSOR)
        else:
            self.is_ppm_pre = 0
            self.fragment_ppm = c_double(self.dp.myCFG.C2_PPM_TOL_FRAGMENT)
            self.precursor_ppm = c_double(self.dp.myCFG.C1_PPM_TOL_PRECURSOR)

        # 需要用到的指针
        c_double_list_ion_num = (c_long * len(self.list_ion_num))()
        for n in range(len(self.list_ion_num)):
            c_double_list_ion_num[n] = self.list_ion_num[n]
        p_C_list_ion_num = byref(c_double_list_ion_num)

        # 导入动态链接库
        bird = CDLL(os.getcwd() + '\\' + DLL_SINGLE)
        self.import_DLL(bird)

        # 把python离子索引转成C的
        p_C_ion_index_start = bird.malloc_ion_index(self.dp.myCFG.max_mz_index_size)
        read_path = c_char_p(bytes(self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[1][self.link_type] + '\\' + FILE_NAME_ION_INDEX[1][self.link_type] + '.dc', 'utf-8'))
        bird.read_ion_index_bin(read_path, p_C_ion_index_start, self.dp.myCFG.max_mz_index_size)
        # 把python离子索引转成C的

        for res_index, one_spe in enumerate(list_coarse_data):
            # print('Finish num : ' + str(res_index + 1) + '/' + str(len(list_coarse_data)))
            if res_index % 10000 == 0:
                print("[Info] Finished single peptide coarse search : %.2f%%" %((res_index + 1) / len(list_coarse_data) * 100))
            if one_spe is None:
                continue
            # 一张张谱图来
            # 单肽没有跨交联的

            # 不跨交联的谱图也受母离子质量影响因为电荷数的不同
            no_linker_mz = one_spe[2]
            no_linker_inten = one_spe[3]

            # 不对强度加和所以原始的就是最高的
            raw_max_inten = 0.0
            for one_inten in no_linker_inten:
                if one_inten > raw_max_inten:
                    raw_max_inten = one_inten
            if raw_max_inten <= 0:
                continue
            # 检查母离子质量
            precursor_mass_start = one_spe[0]
            precursor_mass_end = one_spe[0]
            if self.dp.myCFG.is_ppm_pre:
                precursor_mass_start -= (precursor_mass_start * self.dp.myCFG.C1_PPM_TOL_PRECURSOR)
                precursor_mass_end += (precursor_mass_end * self.dp.myCFG.C1_PPM_TOL_PRECURSOR)
            else:
                precursor_mass_start -= self.dp.myCFG.C1_PPM_TOL_PRECURSOR
                precursor_mass_end += self.dp.myCFG.C1_PPM_TOL_PRECURSOR
            pep_mass_start = precursor_mass_start - self.dp.myINI.MASS_PROTON_MONO
            pep_mass_end = precursor_mass_end - self.dp.myINI.MASS_PROTON_MONO
            if pep_mass_end < self.dp.myCFG.D6_MASS_PEP_LOW:
                continue
            if pep_mass_start > self.dp.myCFG.D7_MASS_PEP_UP:
                continue
            if pep_mass_start < self.dp.myCFG.D6_MASS_PEP_LOW:
                pep_mass_start = self.dp.myCFG.D6_MASS_PEP_LOW
            if pep_mass_end > self.dp.myCFG.D7_MASS_PEP_UP:
                pep_mass_end = self.dp.myCFG.D7_MASS_PEP_UP
            pep_start_index = tool_binary_search_index(self.list_pep_mass, pep_mass_start)
            pep_end_index = tool_binary_search_index(self.list_pep_mass, pep_mass_end)
            if pep_end_index[0] == -1 or pep_start_index[0] == len(self.list_pep_mass):
                continue
            record_pep_start_index = pep_start_index[0]
            record_pep_end_index = pep_end_index[1]
            if pep_start_index[0] == -1:
                record_pep_start_index = 0
            candidate_pep_num = record_pep_end_index - record_pep_start_index
            if candidate_pep_num <= 0:
                continue

            p_C_record_list = bird.malloc_record_score_list(candidate_pep_num, 3)
            self.__captain_get_peptide_match_score(no_linker_mz, no_linker_inten, raw_max_inten, record_pep_start_index, record_pep_end_index, bird, p_C_ion_index_start, p_C_list_ion_num, p_C_record_list)

            list_match_num = p_C_record_list[0: candidate_pep_num * 1]
            list_BM25_score = p_C_record_list[candidate_pep_num * 2: candidate_pep_num * 3]

            bird.free_record_score_list(p_C_record_list)
            sort_result_index = np.argsort(list_BM25_score)  # 把结果分数排序
            self.__soldier_get_one_peptide_result(list_BM25_score, list_match_num, record_pep_start_index, sort_result_index, res_index, coarse_peptide, list_useful_peptide)

        print("[Info] Finished coarse search : 100%")
        bird.free_ion_index(self.dp.myCFG.max_mz_index_size, p_C_ion_index_start)

    def __captain_get_peptide_match_score(self, no_linker_mz, no_linker_inten, raw_max_inten, record_pep_start_index, record_pep_end_index, bird, p_C_ion_index_start, p_C_list_ion_num, p_C_record_list):

        spe_no_linker_mz = (c_long * len(no_linker_mz))()
        spe_no_linker_inten = (c_double * len(no_linker_inten))()

        for mz_index, mz in enumerate(no_linker_mz):
            spe_no_linker_mz[mz_index] = mz
            spe_no_linker_inten[mz_index] = no_linker_inten[mz_index] / raw_max_inten

        p_C_spe_no_linker_mz = byref(spe_no_linker_mz)
        p_C_spe_no_linker_inten = byref(spe_no_linker_inten)

        # 每一张谱图都要是C能读的，要转换一下
        bird.get_ion_index_match_num(p_C_spe_no_linker_mz, len(no_linker_mz), p_C_ion_index_start,
                                     self.dp.myCFG.max_mz_index_size, record_pep_start_index, record_pep_end_index,
                                     self.is_ppm_pre, self.fragment_ppm, p_C_record_list, 0)
        bird.get_ion_index_match_num_percent(p_C_list_ion_num, p_C_record_list, 0, 1, record_pep_start_index,
                                             record_pep_end_index)
        bird.get_ion_index_match_score(p_C_spe_no_linker_mz, len(no_linker_mz),
                                       p_C_spe_no_linker_inten, p_C_ion_index_start, self.dp.myCFG.max_mz_index_size,
                                       record_pep_start_index, record_pep_end_index, self.is_ppm_pre, self.fragment_ppm,
                                       p_C_record_list, 1, 2)

    def search_loop_link(self, list_coarse_data, coarse_peptide, list_useful_result):
        self.link_type = 2
        # 记录编码的结果
        self.__captain_coarse_looplink_search(list_coarse_data, coarse_peptide, list_useful_result)
        # 要回详细结果
        self.__captain_get_loop_link_peptide_data(list_useful_result)

    def __captain_coarse_looplink_search(self, list_coarse_data, coarse_peptide, list_useful_peptide):

        linker_mass = self.dp.myLINK.linker_1_data.loop_mass
        # 提前处理的数据
        if self.dp.myCFG.is_ppm_pre:
            self.is_ppm_pre = 1
            self.fragment_ppm = c_double(self.dp.myCFG.C2_PPM_TOL_FRAGMENT)
            self.precursor_ppm = c_double(self.dp.myCFG.C1_PPM_TOL_PRECURSOR)
        else:
            self.is_ppm_pre = 0
            self.fragment_ppm = c_double(self.dp.myCFG.C2_PPM_TOL_FRAGMENT)
            self.precursor_ppm = c_double(self.dp.myCFG.C1_PPM_TOL_PRECURSOR)

        # 需要用到的指针
        c_double_list_ion_num = (c_long * len(self.list_ion_num))()
        for n in range(len(self.list_ion_num)):
            c_double_list_ion_num[n] = self.list_ion_num[n]
        p_C_list_ion_num = byref(c_double_list_ion_num)

        # 导入动态链接库
        bird = CDLL(os.getcwd() + '\\' + DLL_SINGLE)
        self.import_DLL(bird)
        # 提前处理的数据

        # 把python离子索引转成C的
        p_C_ion_index_start = bird.malloc_ion_index(self.dp.myCFG.max_mz_index_size)
        read_path = c_char_p(bytes(self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[1][self.link_type] + '\\' + FILE_NAME_ION_INDEX[1][self.link_type] + '.dc', 'utf-8'))
        bird.read_ion_index_bin(read_path, p_C_ion_index_start, self.dp.myCFG.max_mz_index_size)
        # 把python离子索引转成C的
        for res_index, one_spe in enumerate(list_coarse_data):
            if res_index % 10000 == 0:
                print("[Info] Finished loop-link coarse search : %.2f%%" %((res_index + 1) / len(list_coarse_data) * 100))
            if one_spe is None:
                continue
            # 一张张谱图来
            # 单肽没有跨交联的

            # 不跨交联的谱图也受母离子质量影响因为电荷数的不同
            no_linker_mz = one_spe[2]
            no_linker_inten = one_spe[3]

            # 不对强度加和所以原始的就是最高的
            raw_max_inten = 0.0
            for one_inten in no_linker_inten:
                if one_inten > raw_max_inten:
                    raw_max_inten = one_inten
            if raw_max_inten <= 0:
                continue

            precursor_mass_start = one_spe[0]
            precursor_mass_end = one_spe[0]
            if self.dp.myCFG.is_ppm_pre:
                precursor_mass_start -= (precursor_mass_start * self.dp.myCFG.C1_PPM_TOL_PRECURSOR)
                precursor_mass_end += (precursor_mass_end * self.dp.myCFG.C1_PPM_TOL_PRECURSOR)
            else:
                precursor_mass_start -= self.dp.myCFG.C1_PPM_TOL_PRECURSOR
                precursor_mass_end += self.dp.myCFG.C1_PPM_TOL_PRECURSOR

            pep_mass_start = precursor_mass_start - self.dp.myINI.MASS_PROTON_MONO - linker_mass
            pep_mass_end = precursor_mass_end - self.dp.myINI.MASS_PROTON_MONO - linker_mass

            if pep_mass_end < self.dp.myCFG.D6_MASS_PEP_LOW:
                continue
            if pep_mass_start > self.dp.myCFG.D7_MASS_PEP_UP:
                continue

            if pep_mass_start < self.dp.myCFG.D6_MASS_PEP_LOW:
                pep_mass_start = self.dp.myCFG.D6_MASS_PEP_LOW
            if pep_mass_end > self.dp.myCFG.D7_MASS_PEP_UP:
                pep_mass_end = self.dp.myCFG.D7_MASS_PEP_UP

            pep_start_index = tool_binary_search_index(self.list_pep_mass, pep_mass_start)
            pep_end_index = tool_binary_search_index(self.list_pep_mass, pep_mass_end)

            record_pep_start_index = pep_start_index[0]
            record_pep_end_index = pep_end_index[1]
            if pep_start_index[0] == -1:
                record_pep_start_index = 0

            if pep_end_index[0] == -1 or pep_start_index[0] == len(self.list_pep_mass):
                continue

            candidate_pep_num = record_pep_end_index - record_pep_start_index
            if candidate_pep_num <= 0:
                continue

            p_C_record_list = bird.malloc_record_score_list(candidate_pep_num, 3)
            self.__captain_get_peptide_match_score(no_linker_mz, no_linker_inten, raw_max_inten, record_pep_start_index, record_pep_end_index, bird, p_C_ion_index_start, p_C_list_ion_num, p_C_record_list)

            list_match_num = p_C_record_list[0: candidate_pep_num * 1]
            list_BM25_score = p_C_record_list[candidate_pep_num * 2: candidate_pep_num * 3]
            bird.free_record_score_list(p_C_record_list)
            sort_result_index = np.argsort(list_BM25_score)  # 把结果分数排序
            self.__soldier_get_one_peptide_result(list_BM25_score, list_match_num, record_pep_start_index, sort_result_index, res_index, coarse_peptide, list_useful_peptide)

        print("[Info] Finished coarse search : 100%")
        bird.free_ion_index(self.dp.myCFG.max_mz_index_size, p_C_ion_index_start)

    def __soldier_get_one_peptide_result(self, list_BM25_score, list_match_num, pep_start_index, sort_result_index, res_index, coarse_peptide, list_candidate_pep_data):

        # 因为都是直接根据母离子质量限制，loop和单肽一样的
        candidate_peptdie_num = 0
        for j in range(len(sort_result_index) - 1, -1, -1):
            peptide_score = list_BM25_score[sort_result_index[j]]
            one_match_ion_num = list_match_num[sort_result_index[j]]
            if one_match_ion_num < 2:  # 这就没有匹配比较的意义了
                break
            if peptide_score < 1:
                break
            if candidate_peptdie_num > self.dp.myCFG.top_N:
                break
            # 到这一步说明这张谱图找到了结果
            coarse_peptide[res_index].append(sort_result_index[j] + pep_start_index)
            list_candidate_pep_data[sort_result_index[j] + pep_start_index] = 1
            candidate_peptdie_num += 1

    def __captain_get_single_peptide_data(self, list_useful_result):

        functionPickle = CFunctionPickle()
        path_protein_pkl_file = os.path.join(self.dp.myCFG.A4_PATH_FASTA_EXPORT, PRO_PKL_FILE)
        protein_list = functionPickle.load_pkl_to_data(path_protein_pkl_file)

        path_peptide2protein_pkl_file = os.path.join(self.dp.myCFG.A4_PATH_FASTA_EXPORT, PEP_PRO_INDEX_FILE)
        peptide2protein = functionPickle.load_pkl_to_data(path_peptide2protein_pkl_file)

        # 由编码获取肽段的详细信息
        pkl_index = 0
        stop_flag_code = self.list_code_split_pkl_file[pkl_index]

        peptide_list = functionPickle.load_pkl_to_data(self.dp.LIST_MASS_INDEX_FILE[pkl_index])
        one_pkl_ind_list = functionPickle.load_pkl_to_data(self.dp.LIST_MASS_IND_FILE[pkl_index])

        for i, one_pepcode in enumerate(self.list_pep_data_pkl):
            if i % 20000 == 0:
                print("[Info] Getting single peptides data %.4f%%" % ((i + 1) / len(self.list_pep_data_pkl) * 100))
            if list_useful_result[i] == -1:
               continue
            if one_pepcode >= stop_flag_code:
                while one_pepcode >= stop_flag_code:
                    pkl_index = pkl_index + 1
                    stop_flag_code = self.list_code_split_pkl_file[pkl_index]
                del peptide_list
                del one_pkl_ind_list
                peptide_list = functionPickle.load_pkl_to_data(self.dp.LIST_MASS_INDEX_FILE[pkl_index])
                one_pkl_ind_list = functionPickle.load_pkl_to_data(self.dp.LIST_MASS_IND_FILE[pkl_index])
            pep_data = tool_new_get_ion_data(one_pepcode)
            if pep_data[1] == 0:
                pep_info = peptide_list[pep_data[2]]
            else:
                pep_info = peptide_list[one_pkl_ind_list[pep_data[1] - 1] + pep_data[2]]
            coarse_data = CSinglePeptideCoarseData()
            op_fill_CSingleCoarseData(coarse_data, pep_info, protein_list, peptide2protein, self.dp.myMOD)
            list_useful_result[i] = coarse_data

    def __captain_get_loop_link_peptide_data(self, list_useful_result):

        functionPickle = CFunctionPickle()
        path_protein_pkl_file = os.path.join(self.dp.myCFG.A4_PATH_FASTA_EXPORT, PRO_PKL_FILE)
        protein_list = functionPickle.load_pkl_to_data(path_protein_pkl_file)

        path_peptide2protein_pkl_file = os.path.join(self.dp.myCFG.A4_PATH_FASTA_EXPORT, PEP_PRO_INDEX_FILE)
        peptide2protein = functionPickle.load_pkl_to_data(path_peptide2protein_pkl_file)

        # 由编码获取肽段的详细信息
        pkl_index = 0
        stop_flag_code = self.list_code_split_pkl_file[pkl_index]

        functionPickle = CFunctionPickle()
        peptide_list = functionPickle.load_pkl_to_data(self.dp.LIST_MASS_INDEX_FILE[pkl_index])
        one_pkl_ind_list = functionPickle.load_pkl_to_data(self.dp.LIST_MASS_IND_FILE[pkl_index])

        for i, one_pepcode in enumerate(self.list_pep_data_pkl):
            if i % 20000 == 0:
                print("[Info] Getting loop-link peptides data %.4f%%" % ((i + 1) / len(self.list_pep_data_pkl) * 100))
            if list_useful_result[i] == -1:
               continue
            if one_pepcode >= stop_flag_code:
                while one_pepcode >= stop_flag_code:
                    pkl_index = pkl_index + 1
                    stop_flag_code = self.list_code_split_pkl_file[pkl_index]
                del peptide_list
                del one_pkl_ind_list
                peptide_list = functionPickle.load_pkl_to_data(self.dp.LIST_MASS_INDEX_FILE[pkl_index])
                one_pkl_ind_list = functionPickle.load_pkl_to_data(self.dp.LIST_MASS_IND_FILE[pkl_index])
            pep_data = tool_new_get_ion_data(one_pepcode)
            if pep_data[1] == 0:
                pep_info = peptide_list[pep_data[2]]
            else:
                pep_info = peptide_list[one_pkl_ind_list[pep_data[1] - 1] + pep_data[2]]
            link_site = self.list_link_site[i]
            link_site_1 = link_site[0][0]
            link_site_2 = link_site[1][0]

            coarse_data = CSinglePeptideCoarseData()
            op_fill_CSingleCoarseData(coarse_data, pep_info, protein_list, peptide2protein, self.dp.myMOD, True, link_site_1, link_site_2)
            list_useful_result[i] = coarse_data

class CFunctionOnePeptideSearchFine:

    link_type = 1

    def __init__(self, inputDP):

        self.dp = inputDP

    def update_class_variable(self, input_type):

        self.link_type = input_type

    def search(self, protein_list, list_fine_data, coarse_peptide, list_useful_result, s2pep_res):

        if self.link_type == 1:

            self.search_singlepeptide(protein_list, list_fine_data, coarse_peptide, list_useful_result, s2pep_res)

        elif self.link_type == 2:

            self.search_loop_link(protein_list, list_fine_data, coarse_peptide, list_useful_result, s2pep_res)

        elif self.link_type == 3:

            pass

    def search_singlepeptide(self, protein_list, list_fine_data, coarse_peptide, list_useful_result, s2pep_res):

        for res_index, one_spectrum in enumerate(list_fine_data):
            one_spectrum_alpha_res = coarse_peptide[res_index]

            if res_index % 10000 == 0:
                print("[Info] Finished single peptide fine search : %.2f%%" % ((res_index + 1) / len(list_fine_data) * 100))
            if one_spectrum is None:
                continue
            candidate_spectrum_result = []

            self.__captain_singlepeptide_fine_search(protein_list, one_spectrum_alpha_res, list_useful_result, one_spectrum, candidate_spectrum_result)

            if len(candidate_spectrum_result) != 0:
                cmpfun = operator.attrgetter('pep_score')
                candidate_spectrum_result.sort(key=cmpfun, reverse=True)
                max_crosspeptide_score = candidate_spectrum_result[0].pep_score
                out_result_index = 0
                target_result_index = 0
                score_index = 0
                while(1):
                    if score_index >= len(candidate_spectrum_result) or max_crosspeptide_score - candidate_spectrum_result[score_index].pep_score > 0.01:
                        out_result_index = score_index
                        break
                    if candidate_spectrum_result[score_index].pep_data.is_target == 1:
                        target_result_index = score_index
                    score_index += 1

                if len(candidate_spectrum_result) > out_result_index:
                    delta_score = candidate_spectrum_result[target_result_index].pep_score - candidate_spectrum_result[out_result_index].pep_score
                else:
                    delta_score = candidate_spectrum_result[target_result_index].pep_score
                write_result = CWriteResultSinglePeptide(one_spectrum)
                op_fill_CWriteSinglePeptdidResult(write_result, candidate_spectrum_result[target_result_index], self.dp.myMOD, self.dp.myLINK, protein_list, delta_score)
                s2pep_res.append(write_result)
        print("[Info] Finished single peptide fine search : 100%%")

    def __captain_singlepeptide_fine_search(self, protein_list, one_spectrum_peptide_res, list_useful_result, one_spectrum, candidate_spectrum_result):

        for i in range(len(one_spectrum_peptide_res)):  # 与结果肽段对应索引

            pep_data_index = one_spectrum_peptide_res[i]
            pep_data = list_useful_result[pep_data_index]
            pep_sq = protein_list.sq[pep_data.pro_index][pep_data.start_pos: pep_data.end_pos]

            pep_aa_mass_list = self.__soldier_create_aa_list(pep_sq, pep_data.mod_site_list)
            self.__soldier_one_single_peptide_fine_score(pep_aa_mass_list, pep_data, one_spectrum, candidate_spectrum_result)

    def __soldier_create_aa_list(self, sq, mod_site_list):
        # create the mass list of peptide sequence
        # e.g., peptide is "ACEDFK with C+57 modification"
        # return [M(A), M(C)+57, M(E), ..., M(K)] in which M(X) is the mass of X
        aa = [0.0 for i in range(len(sq))]
        for aa_index, aa_mod_num in enumerate(mod_site_list):
            if aa_mod_num >= self.dp.myMOD.fix_mod_num:
                add_mod_mass = self.dp.myINI.DIC_MOD[self.dp.myMOD.var_mod_list[aa_mod_num - self.dp.myMOD.fix_mod_num]].mass
            elif aa_mod_num > -1:
                add_mod_mass = self.dp.myINI.DIC_MOD[self.dp.myMOD.fix_mod_list[aa_mod_num]].mass
            else:
                add_mod_mass = 0.0
            if aa_index == 0:
                aa[0] += add_mod_mass
            elif aa_index == len(sq) + 1:
                aa[aa_index - 2] += add_mod_mass
            else:
                aa[aa_index - 1] += self.dp.myINI.DIC_AA[sq[aa_index - 1]] + add_mod_mass
        return aa

    def __soldier_one_single_peptide_fine_score(self, pep_aa_mass_list, pep_data, one_spectrum, candidate_spectrum_result):

        ion_charge = self.dp.myCFG.max_fragment_charge
        pep_list_b_ion, pep_list_y_ion = tool_generate_singlepepdite_by_ion(pep_aa_mass_list, ion_charge)

        pep_ion, pep_same_num = [], []  # 将b y离子合并
        tool_set_two_list(sorted(pep_list_b_ion), sorted(pep_list_y_ion), 0.0005, pep_ion, pep_same_num)

        pep_continue_score, pep_cover = [], []
        pep_continue_data = op_get_peptide_by_continue_data(pep_list_b_ion, pep_list_y_ion, pep_data.sq_length - 1, one_spectrum, self.dp.myCFG.is_ppm_fra, self.dp.myCFG.C2_PPM_TOL_FRAGMENT, 1)
        op_get_peptide_continue_score(pep_continue_data, pep_continue_score, pep_cover)

        pep_score = op_get_match_ion_score(one_spectrum, pep_ion, pep_same_num, self.dp.myCFG.is_ppm_fra, self.dp.myCFG.C2_PPM_TOL_FRAGMENT, 1)

        # 有效离子肽谱匹配得分
        if len(pep_ion) == 0:
            pep_part_1_score = 0.0
        else:
            pep_part_1_score = pep_score.match_ion_score

        # 连续性分数和覆盖率得分
        pep_part_2_score = 0.0
        for i in range(len(pep_continue_score)):
            pep_part_2_score += pep_continue_score[i]

        if len(pep_ion) == 0:
            pep_part_3_score = 0.0
        else:
            pep_part_3_score = 0.0#pep_score.match_ion_score


        pep_final_score = pep_part_1_score + pep_part_2_score + pep_part_3_score

        pep_mass = pep_data.mass + self.dp.myINI.MASS_PROTON_MONO

        precursor_mass_erro_Da = one_spectrum.mass - pep_mass
        precursor_mass_erro_ppm = precursor_mass_erro_Da / pep_mass * 1e6


        if len(pep_score.list_ppm) == 0:
            pass
        else:
            ppm_sum = np.sum(pep_score.list_Da)
            ppm_average = np.mean(pep_score.list_Da)
            # 求方差
            ppm_var = np.var(pep_score.list_Da)

            try:
                pParseNum = int(one_spectrum.title.split('.')[4])
            except:
                pParseNum = 0

            peptide_feature = CRerankOnePeptideFeature(0,
                                                       pep_final_score,
                                                       ppm_sum, ppm_average, ppm_var,
                                                       pep_score,
                                                       pep_part_2_score,
                                                       abs(precursor_mass_erro_Da),
                                                       0,
                                                       pParseNum)

            if pep_final_score > 0:

                if self.dp.myCFG.is_ppm_pre:
                    if abs(precursor_mass_erro_ppm) < self.dp.myCFG.C1_PPM_TOL_PRECURSOR * 1e6:
                        one_result = CSinglePeptideResult(pep_data, peptide_feature)
                        op_fill_CSinglePeptideResult(one_result, pep_mass, pep_final_score, precursor_mass_erro_Da, precursor_mass_erro_ppm, [pep_score, pep_part_1_score, pep_part_2_score, pep_part_3_score])
                        candidate_spectrum_result.append(one_result)
                else:
                    if abs(precursor_mass_erro_Da) > self.dp.myCFG.C1_PPM_TOL_PRECURSOR * 1e6:
                        pass
                    else:
                        one_result = CSinglePeptideResult(pep_data, peptide_feature)
                        op_fill_CSinglePeptideResult(one_result, pep_mass, pep_final_score, precursor_mass_erro_Da, precursor_mass_erro_ppm, [pep_score, pep_part_1_score, pep_part_2_score, pep_part_3_score])
                        candidate_spectrum_result.append(one_result)

    def search_loop_link(self, protein_list, list_fine_data, coarse_peptide, list_useful_result, s2pep_res):

        for res_index, one_spectrum in enumerate(list_fine_data):
            one_spectrum_res = coarse_peptide[res_index]
            if res_index % 10000 == 0:
                print("[Info] Finished loop-link fine search : %.2f%%" % ((res_index + 1) / len(list_fine_data) * 100))
            if one_spectrum is None:
                continue
            candidate_spectrum_result = []

            self.__captain_loop_link_fine_search(protein_list, one_spectrum_res, list_useful_result, one_spectrum, candidate_spectrum_result)

            if len(candidate_spectrum_result) != 0:
                cmpfun = operator.attrgetter('pep_score')
                candidate_spectrum_result.sort(key=cmpfun, reverse=True)
                max_crosspeptide_score = candidate_spectrum_result[0].pep_score
                out_result_index = 0
                target_result_index = 0
                score_index = 0
                while (1):
                    if score_index >= len(candidate_spectrum_result) or max_crosspeptide_score - candidate_spectrum_result[score_index].pep_score > 0.01:
                        out_result_index = score_index
                        break
                    if candidate_spectrum_result[score_index].pep_data.is_target == 1:
                        target_result_index = score_index
                    score_index += 1

                if len(candidate_spectrum_result) > out_result_index:
                    delta_score = candidate_spectrum_result[target_result_index].pep_score - candidate_spectrum_result[out_result_index].pep_score
                else:
                    delta_score = candidate_spectrum_result[target_result_index].pep_score
                write_result = CWriteResultSinglePeptide(one_spectrum, True)
                op_fill_CWriteSinglePeptdidResult(write_result, candidate_spectrum_result[target_result_index], self.dp.myMOD, self.dp.myLINK, protein_list, delta_score)
                if self.dp.myLINK.linker_1_data.linker_name == 'TG' and write_result.peptide_sq.startswith('Q'):
                    pass
                else:
                    s2pep_res.append(write_result)
        print("[Info] Finished loop-link fine search : 100%%")

    def __captain_loop_link_fine_search(self, protein_list, one_spectrum_peptide_res, list_useful_result, one_spectrum, candidate_spectrum_result):

        for i in range(len(one_spectrum_peptide_res)):  # 与结果肽段对应索引

            pep_data_index = one_spectrum_peptide_res[i]
            pep_data = list_useful_result[pep_data_index]
            pep_sq = protein_list.sq[pep_data.pro_index][pep_data.start_pos: pep_data.end_pos]

            pep_aa_mass_list = self.__soldier_create_aa_list(pep_sq, pep_data.mod_site_list)

            self.__soldier_one_loop_link_fine_score(pep_aa_mass_list, pep_data, pep_data.loop_site_1, pep_data.loop_site_2, one_spectrum, candidate_spectrum_result)

    def __soldier_one_loop_link_fine_score(self, pep_aa_mass_list, pep_data, link_site_1, link_site_2, one_spectrum, candidate_spectrum_result):

        ion_charge = self.dp.myCFG.max_fragment_charge
        pep_list_b_ion, pep_list_y_ion = tool_generate_loop_type_by_ion(pep_aa_mass_list, link_site_1, link_site_2, self.dp.myLINK.linker_1_data.loop_mass, ion_charge)

        pep_ion, pep_same_num = [], []  # 将b y离子合并
        tool_set_two_list(sorted(pep_list_b_ion), sorted(pep_list_y_ion), 0.0005, pep_ion, pep_same_num)

        pep_continue_score, pep_cover = [], []
        ion_len = pep_data.sq_length - 1 - (link_site_2 - link_site_1)
        pep_continue_data = op_get_peptide_by_continue_data(pep_list_b_ion, pep_list_y_ion, ion_len, one_spectrum, self.dp.myCFG.is_ppm_fra, self.dp.myCFG.C2_PPM_TOL_FRAGMENT, 1)
        op_get_peptide_continue_score(pep_continue_data, pep_continue_score, pep_cover)

        pep_score = op_get_match_ion_score(one_spectrum, pep_ion, pep_same_num, self.dp.myCFG.is_ppm_fra, self.dp.myCFG.C2_PPM_TOL_FRAGMENT, 1)

        # 有效离子肽谱匹配得分
        if len(pep_ion) == 0:
            pep_part_1_score = 0.0
        else:
            pep_part_1_score = pep_score.match_ion_score * len(pep_score.list_ppm) / len(pep_ion)

        # 连续性分数和覆盖率得分
        pep_part_2_score = 0.0
        for i in range(len(pep_continue_score)):
            pep_part_2_score += pep_continue_score[i]

        if len(pep_ion) == 0:
            pep_part_3_score = 0.0
        else:
            pep_part_3_score = len(pep_score.list_ppm) / len(pep_ion) * pep_score.match_ion_score

        pep_final_score = pep_part_1_score + pep_part_2_score + pep_part_3_score

        pep_mass = pep_data.mass + self.dp.myINI.MASS_PROTON_MONO

        precursor_mass_erro_Da = one_spectrum.mass - pep_mass - self.dp.myLINK.linker_1_data.loop_mass
        precursor_mass_erro_ppm = precursor_mass_erro_Da / pep_mass * 1e6

        if len(pep_score.list_ppm) == 0:
            pass
        else:
            ppm_sum = np.sum(pep_score.list_Da)
            ppm_average = np.mean(pep_score.list_Da)
            # 求方差
            ppm_var = np.var(pep_score.list_Da)

            try:
                pParseNum = int(one_spectrum.title.split('.')[4])
            except:
                pParseNum = 0

            peptide_feature = CRerankOnePeptideFeature(0,
                                                       pep_final_score,
                                                       ppm_sum,ppm_average, ppm_var,
                                                       pep_score,
                                                       pep_part_2_score,
                                                       abs(precursor_mass_erro_ppm),
                                                       0,
                                                       pParseNum)

            if pep_final_score > 0:

                if self.dp.myCFG.is_ppm_pre:
                    if abs(precursor_mass_erro_ppm) < self.dp.myCFG.C1_PPM_TOL_PRECURSOR * 1e6:
                        one_result = CSinglePeptideResult(pep_data, peptide_feature, True, link_site_1, link_site_2)
                        op_fill_CSinglePeptideResult(one_result, pep_mass, pep_final_score, precursor_mass_erro_Da,
                                                     precursor_mass_erro_ppm,
                                                     [pep_score, pep_part_1_score, pep_part_2_score, pep_part_3_score])
                        candidate_spectrum_result.append(one_result)
                else:
                    if abs(precursor_mass_erro_Da) > self.dp.myCFG.C1_PPM_TOL_PRECURSOR * 1e6:
                        pass
                    else:
                        one_result = CSinglePeptideResult(pep_data, peptide_feature, True, link_site_1, link_site_2)
                        op_fill_CSinglePeptideResult(one_result, pep_mass, pep_final_score, precursor_mass_erro_Da,
                                                     precursor_mass_erro_ppm,
                                                     [pep_score, pep_part_1_score, pep_part_2_score, pep_part_3_score])
                        candidate_spectrum_result.append(one_result)
