# -*- mode: python ; coding: utf-8 -*-..
import multiprocessing
import pickle

from MSFunction import CFunctionPickle
from MSFunctionRerank import CFunctionOutTwoPeptidePSM
from MSFunctionMutiProcess import CFunctionThreadCoarseSearch
from MSData import CSpectrum, CPeak
from MSDataResult import COnlyCrossCoarseData, CRerankTwoPeptideFeature
from MSDataResult import COnlyCrossResult, CWriteResultOnlyCross
from MSDataResult import CCrossCrossResult, CWriteResultCrossCross
from MSDataResult import CCrossLoopResult, CWriteResultCrossLoop
from MSSysterm import PRO_PKL_FILE, PEP_PRO_INDEX_FILE
from MSSysterm import FILE_NAME_PSM_RESULT
from MSSysterm import FOLDER_INDEX, FILE_NAME_ION_INDEX, FILE_NAME_ION_NUM, FILE_NAME_PEP_MASS, FILE_NAME_PEP_PKL, \
    FILE_NAME_LINK_SITE
from MSOperator import op_create_peaks_index, op_fill_COnlyCrossCoarseData
from MSOperator import op_fill_COnlyCrossResult, op_fill_CWriteOnlyCrossResult
from MSOperator import op_fill_CCrossCrossResult, op_fill_CWriteCrossCrossResult
from MSOperator import op_fill_CCrossLoopResult, op_fill_CWriteCrossLoopResult

from MSOperatorScore import op_get_peptide_by_continue_data, op_get_peptide_continue_score, op_get_match_ion_score
from MSTimer import Timer
from MSTool import tool_generate_only_cross_by_ion, tool_set_two_list, tool_generate_cross_cross_by_ion, \
    tool_generate_cross_loop_by_ion
from MSTool_code import tool_new_get_ion_code, tool_new_get_ion_data
from MSLogging import logGetError, Logger

from ctypes import *
import time
import os
import numpy as np


class CFunctionTwoPeptideSearch:
    link_type = 1

    def __init__(self, inputDP):

        self.dp = inputDP

    def update_class_variable(self, link_type):

        self.link_type = link_type

    def search(self, id_psm_name=None, id_psm_scan=None):

        # 粗打分的参数初始化
        functionSearchCoarse = CFunctionTwoPeptideSearchCoarse_C(self.dp)
        pep_num = functionSearchCoarse.update_class_variable(self.link_type)
        # 细打分的参数初始化
        functionSearchFine = CFunctionTwoPeptideSearchFine(self.dp)
        functionSearchFine.update_class_variable(self.link_type)

        functionPickle = CFunctionPickle()
        all_match_result = []
        # 一次粗打分和细打分调用时给出一个mgf的检索结果
        for mgf_index in range(self.dp.myMS2.ms2_num):
            dataMS2 = functionPickle.load_pkl_to_data(self.dp.myMS2.trans_ms2_file[mgf_index])
            # 把谱图弄成列表形式
            list_coarse_spectrum_data = [None for i in range(dataMS2.MS2_NUM)]
            list_fine_spectrum_data = [None for i in range(dataMS2.MS2_NUM)]
            self.trans_dataMS2_to_CSpectrum_list(dataMS2, list_coarse_spectrum_data, list_fine_spectrum_data,
                                                 id_psm_name, id_psm_scan)
            # for i in range(len(list_coarse_spectrum_data)):
            #     if list_coarse_spectrum_data[i] is None:
            #         continue
            print("[Info]#Spectra in {0} is {1}".format(self.dp.myMS2.trans_ms2_file[mgf_index], dataMS2.MS2_NUM))
            del dataMS2
            time_1 = time.time()
            # 用来记录结果
            coarse_res = []
            list_useful_result = [-1 for i in range(pep_num)]
            # 粗打分
            functionSearchCoarse.search(list_coarse_spectrum_data, coarse_res, list_useful_result)
            time_2 = time.time()
            # 细打分
            s2pep_res = []
            functionSearchFine.search(list_fine_spectrum_data, coarse_res, list_useful_result, s2pep_res)
            all_match_result += s2pep_res
            time_3 = time.time()
            print("[Info] Coarse search using time : %.4f s ." % (time_2 - time_1))
            print("[Info] Fine search using time : %.4f s ." % (time_3 - time_2))

        functionOutTwoPeptidePSM = CFunctionOutTwoPeptidePSM(self.dp)

        result_folder = self.dp.myCFG.A5_PATH_RESULT_EXPORT
        result_folder_pkl = result_folder + '\\tmp\\'
        if not os.path.exists(result_folder_pkl):
            os.makedirs(result_folder_pkl)

        write_result_path = result_folder + FILE_NAME_PSM_RESULT[2][self.link_type] + '.txt'
        write_result_pkl_path = result_folder_pkl + FILE_NAME_PSM_RESULT[2][self.link_type] + '.pkl'
        functionOutTwoPeptidePSM.update_class_variable(1, 0, 0)  # FDR阈值，排序方式 清除反库
        functionOutTwoPeptidePSM.outputOnlyCrossPSM(all_match_result, write_result_pkl_path, write_result_path,
                                                    self.dp.myCFG.separate_FDR)

    def trans_dataMS2_to_CSpectrum_list(self, dataMS2, list_coarse_data, list_fine_data, id_psm_name=None,
                                        id_psm_scan=None):
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
        # record_scan = []
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
                        dataMS2.MATRIX_PEAK_MOZ_LINKER[ms2_scan][precursor_moz_index],
                        dataMS2.MATRIX_PEAK_INT_LINKER[ms2_scan][precursor_moz_index],
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
                            dataMS2.MATRIX_PEAK_MOZ_LINKER[ms2_scan][precursor_moz_index],
                            dataMS2.MATRIX_PEAK_INT_LINKER[ms2_scan][precursor_moz_index],
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
                        # record_scan.append(ms2_scan)
                    else:
                        if dataMS2.MATRIX_FILE_NAME[ms2_scan][precursor_moz_index] in id_scan_list[ms2_scan]:
                            pass
                        else:
                            list_coarse_data[spectrum_num] = [
                                dataMS2.MATRIX_PRECURSOR_MOZ[ms2_scan][precursor_moz_index],
                                dataMS2.MATRIX_PRECURSOR_CHARGE[ms2_scan][precursor_moz_index],
                                dataMS2.MATRIX_PEAK_MOZ_NO_LINKER[ms2_scan][precursor_moz_index],
                                dataMS2.MATRIX_PEAK_INT_NO_LINKER[ms2_scan][precursor_moz_index],
                                dataMS2.MATRIX_PEAK_MOZ_LINKER[ms2_scan][precursor_moz_index],
                                dataMS2.MATRIX_PEAK_INT_LINKER[ms2_scan][precursor_moz_index],
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
                            # record_scan.append(ms2_scan)
                    spectrum_num += 1
        # with open(self.dp.myCFG.A4_PATH_FASTA_EXPORT + 'cross_spe_scan.txt', 'w') as f:
        #     for scan in record_scan:
        #         f.write(str(scan))
        #         f.write('\n')
        # exit(0)

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


class CFunctionTwoPeptideSearchCoarse_C:
    # 类变量，即使多个它们只存一份就行，谁在使用之前就做一次更新
    # mz_index_data = []  # 这个是离子索引，mz-> [肽段编码]
    list_pep_data_pkl = []  # 这个是记录满足交联条件的肽段索引，index->肽段编码
    list_ion_num = []  # 这个是记录肽段的碎片离子数目，index->离子数
    list_pep_mass = []  # 这个是记录肽段质量的，不包括交联剂质量
    list_link_site = []

    list_code_split_pkl_file = []  # 这个是用来记录每个区间的肽段到哪个范围

    link_type = 1

    linker_mass = 0.0

    debug_pep_sq_data = []

    def __init__(self, inputDP):

        self.dp = inputDP

    def update_class_variable(self, link_type):
        # 初始化粗打分参数
        self.link_type = link_type
        self.__captain_load_index()
        self.__captain_get_pkl_file_split_code()
        pep_num = len(self.list_pep_data_pkl)
        return pep_num

    def __captain_load_index(self):

        functionPickle = CFunctionPickle()

        path = self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[2][self.link_type] + '\\' + \
               FILE_NAME_ION_NUM[2][self.link_type] + '.pkl'
        self.list_ion_num = functionPickle.load_pkl_to_data(path)

        path = self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[2][self.link_type] + '\\' + \
               FILE_NAME_PEP_MASS[2][self.link_type] + '.pkl'
        self.list_pep_mass = functionPickle.load_pkl_to_data(path)

        path = self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[2][self.link_type] + '\\' + \
               FILE_NAME_PEP_PKL[2][self.link_type] + '.pkl'
        self.list_pep_data_pkl = functionPickle.load_pkl_to_data(path)

        path = self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[2][self.link_type] + '\\' + \
               FILE_NAME_LINK_SITE[2][self.link_type] + '.pkl'
        self.list_link_site = functionPickle.load_pkl_to_data(path)

    def __captain_get_pkl_file_split_code(self):

        for i in range(int((
                                   self.dp.myCFG.D7_MASS_PEP_UP - self.dp.myCFG.D6_MASS_PEP_LOW) / self.dp.myCFG.D10_INDEX_SPLIT_MASS)):
            my_code = tool_new_get_ion_code(i + 1, 0, 0, 0)
            self.list_code_split_pkl_file.append(my_code)

    def search(self, list_coarse_data, coarse_res, list_useful_result):
        # 检索
        if self.link_type == 1:

            self.linker_mass = self.dp.myLINK.linker_1_data.cross_mass
            self.search_only_cross(list_coarse_data, coarse_res, list_useful_result)

        elif self.link_type == 2:

            pass

        elif self.link_type == 3:

            pass

    def search_only_cross(self, list_coarse_data, coarse_res, list_useful_result):

        # 记粗打分结果，coarse_alpha_res和coarse_beta_res是记录结果来自于所有结果里的哪个位置
        self.link_type = 1

        read_path = c_char_p(bytes(
            self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[2][self.link_type] + '\\' + FILE_NAME_ION_INDEX[2][
                self.link_type] + '.dc', 'utf-8'))

        list_pepcode2resultindex = [-1 for i in range(len(self.list_pep_data_pkl))]  # 记录记结果所在位置

        # functionMutiPreocessCoarseSearch = CFunctionThreadCoarseSearch(self.dp, self.list_ion_num, self.list_pep_mass, self.list_pep_data_pkl, self.list_link_site, self.linker_mass, list_pepcode2resultindex)
        # functionMutiPreocessCoarseSearch.MutiThreadOnlyCrossCoarseSearch(list_coarse_data, coarse_res, read_path, False)

        coarse_path = r"F:\wjq0619\pkl\spectrum_0_20000.npz_elink_result1111.pkl"
        # 要回详细结果
        self.__load_coarse_data(coarse_path, coarse_res, list_pepcode2resultindex)

        self.__captain_get_alpha_beta_data(list_pepcode2resultindex, list_useful_result)

    def __load_coarse_data(self, path, coarse_res, list_pepcode2resultindex):
        f = open(path, 'rb')
        coarse_data = pickle.load(f)
        f.close()
        Logger.debug('load coarse data')
        for one_coarse_data in coarse_data:
            for one_result in one_coarse_data:
                list_pepcode2resultindex[one_result[0]] = 1
                list_pepcode2resultindex[one_result[1]] = 1

        coarse_res += coarse_data

    def __captain_get_alpha_beta_data(self, list_pepcode2resultindex, list_useful_result):
        # 获取肽段的序列结果
        # 只要有用到的结果,这些结果只用一份，通过索引号获取
        time_start = time.time()

        functionPickle = CFunctionPickle()
        path_protein_pkl_file = os.path.join(self.dp.myCFG.A4_PATH_FASTA_EXPORT, PRO_PKL_FILE)
        protein_list = functionPickle.load_pkl_to_data(path_protein_pkl_file)

        path_peptide2protein_pkl_file = os.path.join(self.dp.myCFG.A4_PATH_FASTA_EXPORT, PEP_PRO_INDEX_FILE)
        peptide2protein = functionPickle.load_pkl_to_data(path_peptide2protein_pkl_file)
        # 这里是从小到大遍历编码，对一个pkl文件内的肽段信息要完分到对应位置后再到下一个文件
        # 由编码获取肽段的详细信息
        pkl_index = 0
        stop_flag_code = self.list_code_split_pkl_file[pkl_index]
        peptide_list = functionPickle.load_pkl_to_data(self.dp.LIST_MASS_INDEX_FILE[pkl_index])
        one_pkl_ind_list = functionPickle.load_pkl_to_data(self.dp.LIST_MASS_IND_FILE[pkl_index])
        cur_need_num = 0
        for i, one_pepcode in enumerate(self.list_pep_data_pkl):
            if i % 20000 == 0:
                print("[Info] Getting cross-link peptides data %.4f%%" % ((i + 1) / len(self.list_pep_data_pkl) * 100))
            if list_pepcode2resultindex[i] == -1:
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
            coarse_data = COnlyCrossCoarseData()
            op_fill_COnlyCrossCoarseData(coarse_data, pep_info, protein_list, peptide2protein, self.dp.myMOD)
            list_useful_result[i] = coarse_data
            cur_need_num += 1
            j = i
            if j > 0:
                while (self.list_pep_data_pkl[j] == self.list_pep_data_pkl[j - 1]):
                    list_useful_result[j - 1] = coarse_data
                    j -= 1
                    if j == 0:
                        break
            j = i
            if j < len(self.list_pep_data_pkl) - 1:
                while (self.list_pep_data_pkl[j] == self.list_pep_data_pkl[j + 1]):
                    list_useful_result[j + 1] = coarse_data
                    j += 1
                    if j == len(self.list_pep_data_pkl) - 1:
                        break

        time_end = time.time()
        print("[Info] Getting cross-link peptides data 100%")
        print("[Info] Get peptides data using time : %.4f s ." % (time_end - time_start))


class CFunctionTwoPeptideSearchFine:
    link_type = 2
    cross_mass = 0.0

    def __init__(self, inputDP):
        self.dp = inputDP
        # self.pool = multiprocessing.Pool(processes=16)

    def update_class_variable(self, link_type):

        self.link_type = link_type

        functionPickle = CFunctionPickle()
        path = self.dp.myCFG.A4_PATH_FASTA_EXPORT + '\\' + FOLDER_INDEX[2][self.link_type] + '\\' + \
               FILE_NAME_LINK_SITE[2][self.link_type] + '.pkl'
        self.list_link_site = functionPickle.load_pkl_to_data(path)

    def search(self, list_fine_data, coarse_res, list_useful_result, s2pep_res):

        if self.link_type == 1:

            self.search_only_cross(list_fine_data, coarse_res, list_useful_result, s2pep_res)

        elif self.link_type == 2:

            # self.search_cross_cross(list_fine_data, coarse_res, list_useful_result, s2pep_res)
            pass
        elif self.link_type == 3:
            pass
            # self.search_cross_loop(list_fine_data, coarse_res, list_useful_result, s2pep_res)

    def multi_process_search(self, list_fine_data, coarse_res, list_useful_result, s2pep_res):
        split = self.__split_fine_search_data(list_fine_data, coarse_res, list_useful_result)
        with multiprocessing.pool.ThreadPool(processes=16) as pool:
            pool.map(self.multi_search_only_cross, split)
            pool.join()

        for one_piece in split:
            s2pep_res += one_piece[3]

    def __split_fine_search_data(self, list_fine_data, coarse_res, list_useful_result):
        split = []

        chunk_size = len(list_fine_data) // self.dp.myCFG.thread_num
        # 把数据分成 thread_num 份
        for i in range(self.dp.myCFG.thread_num):
            split.append([list_fine_data[i * chunk_size: (i + 1) * chunk_size],
                          coarse_res[i * chunk_size: (i + 1) * chunk_size],
                          list_useful_result,
                          []])
        return split

    def multi_search_only_cross(self, input_data):
        list_fine_data = input_data[0]
        coarse_res = input_data[1]
        list_useful_result = input_data[2]
        s2pep_res = input_data[3]
        self.search_only_cross(list_fine_data, coarse_res, list_useful_result, s2pep_res)

    def search_only_cross(self, list_fine_data, coarse_res, list_useful_result, s2pep_res):

        self.cross_mass = self.dp.myLINK.linker_1_data.cross_mass

        functionPickle = CFunctionPickle()
        path_protein_pkl_file = os.path.join(self.dp.myCFG.A4_PATH_FASTA_EXPORT, PRO_PKL_FILE)
        protein_list = functionPickle.load_pkl_to_data(path_protein_pkl_file)

        timer = Timer("Fine Search")
        for res_index, one_spectrum in enumerate(list_fine_data):
            if res_index % 500 == 0:
                print(
                    "[Info] Finished cross-link fine search : %.2f%%" % ((res_index + 1) / (len(list_fine_data)) * 100))
                timer.print()
                timer.init()

            candidate_spectrum_result = []
            one_spectrum_res = coarse_res[res_index]
            if one_spectrum is None:
                continue

            self.__captain_only_cross_fine_search(protein_list, list_useful_result, one_spectrum_res, one_spectrum,
                                                  candidate_spectrum_result, timer)

            # record_sq = []
            # for one_candidate_spectrum_result in candidate_spectrum_result:
            #     alpha_sq = protein_list.sq[one_candidate_spectrum_result.alpha_pep_data.pro_index][one_candidate_spectrum_result.alpha_pep_data.start_pos: one_candidate_spectrum_result.alpha_pep_data.end_pos]
            #     beta_sq = protein_list.sq[one_candidate_spectrum_result.beta_pep_data.pro_index][one_candidate_spectrum_result.beta_pep_data.start_pos: one_candidate_spectrum_result.beta_pep_data.end_pos]
            #     record_sq.append([alpha_sq, beta_sq, one_candidate_spectrum_result.pep_score])

            if len(candidate_spectrum_result) > 0:
                candidate_spectrum_result.sort(key=lambda candidate_scan_result: candidate_scan_result.pep_score,
                                               reverse=True)
                out_result_index, second_index = self.__captain_get_final_result_index(candidate_spectrum_result)
                if second_index == - 1:
                    delta_score = candidate_spectrum_result[0].pep_score
                else:
                    delta_score = candidate_spectrum_result[out_result_index].pep_score - candidate_spectrum_result[
                        second_index].pep_score
                write_result = CWriteResultOnlyCross(one_spectrum)
                op_fill_CWriteOnlyCrossResult(write_result, candidate_spectrum_result[out_result_index], self.dp.myMOD,
                                              protein_list, self.dp.myLINK.linker_1_data.linker_name)
                write_result.feature.delta_score = delta_score
                s2pep_res.append(write_result)
        print("[Info] Finished cross-link fine search : 100%")

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

    def __captain_only_cross_fine_search(self, protein_list, list_useful_result, one_spectrum_res, one_spectrum,
                                         all_fine_score, timer):

        # record_sq = []
        # for i in range(len(one_spectrum_res)):
        #     alpha_pep_data_index = one_spectrum_res[i][0]
        #     alpha_pep_data = list_useful_result[alpha_pep_data_index]
        #     alpha_pep_sq = protein_list.sq[alpha_pep_data.pro_index][alpha_pep_data.start_pos: alpha_pep_data.end_pos]
        #
        #     beta_pep_data_index = one_spectrum_res[i][1]
        #     beta_pep_data = list_useful_result[beta_pep_data_index]
        #     beta_pep_sq = protein_list.sq[beta_pep_data.pro_index][beta_pep_data.start_pos: beta_pep_data.end_pos]
        #     record_sq.append([alpha_pep_sq, beta_pep_sq, one_spectrum_res[i][2]])
        args = []
        for i in range(len(one_spectrum_res)):  # 与结果肽段对应索引,因为序列一样所以暂时不用管位点

            alpha_pep_data_index = one_spectrum_res[i][0]
            alpha_pep_data = list_useful_result[alpha_pep_data_index]
            alpha_pep_sq = protein_list.sq[alpha_pep_data.pro_index][alpha_pep_data.start_pos: alpha_pep_data.end_pos]
            alpha_pep_aa_mass_list = self.__soldier_create_aa_list(alpha_pep_sq, alpha_pep_data.mod_site_list)
            alpha_pep_site = self.list_link_site[alpha_pep_data_index][0]

            beta_pep_data_index = one_spectrum_res[i][1]
            beta_pep_data = list_useful_result[beta_pep_data_index]
            beta_pep_sq = protein_list.sq[beta_pep_data.pro_index][beta_pep_data.start_pos: beta_pep_data.end_pos]
            beta_pep_aa_mass_list = self.__soldier_create_aa_list(beta_pep_sq, beta_pep_data.mod_site_list)
            beta_pep_site = self.list_link_site[beta_pep_data_index][0]
            args.append((alpha_pep_aa_mass_list, beta_pep_aa_mass_list, alpha_pep_data, beta_pep_data, alpha_pep_site,
                         beta_pep_site, one_spectrum, [], timer, alpha_pep_sq, beta_pep_sq))
            # self.__soldier_one_only_cross_fine_score(alpha_pep_aa_mass_list, beta_pep_aa_mass_list, alpha_pep_data,
            #                                          beta_pep_data, alpha_pep_site, beta_pep_site, one_spectrum,
            #                                          all_fine_score, timer, alpha_pep_sq, beta_pep_sq)
        with multiprocessing.Pool(processes=16) as pool:
            res = pool.starmap(self.soldier_one_only_cross_fine_score, args)
            # self.pool.join()
            # pool.join()
            all_fine_score += res

    def __soldier_create_aa_list(self, sq, mod_site_list):
        # create the mass list of peptide sequence
        # e.g., peptide is "ACEDFK with C+57 modification"
        # return [M(A), M(C)+57, M(E), ..., M(K)] in which M(X) is the mass of X
        aa = [0.0 for i in range(len(sq))]
        for aa_index, aa_mod_num in enumerate(mod_site_list):
            if aa_mod_num >= self.dp.myMOD.fix_mod_num:
                add_mod_mass = self.dp.myINI.DIC_MOD[
                    self.dp.myMOD.var_mod_list[aa_mod_num - self.dp.myMOD.fix_mod_num]].mass
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

    def soldier_one_only_cross_fine_score(self, alpha_pep_aa_mass_list, beta_pep_aa_mass_list, alpha_pep_data,
                                          beta_pep_data, alpha_pep_cross_site, beta_pep_cross_site, one_spectrum,
                                          all_fine_score, timer: Timer, alpha_pep_sq='', beta_pep_sq=''):
        # 这是针对交联剂两端位点一致的情况, beta的位点是空的就是一样的，不空就跳转到beta的位点凑
        ion_charge = self.dp.myCFG.max_fragment_charge
        timer.reset()

        alpha_pep_list_b_ion, alpha_pep_list_y_ion = tool_generate_only_cross_by_ion(alpha_pep_aa_mass_list, ion_charge,
                                                                                     alpha_pep_cross_site,
                                                                                     beta_pep_data.mass + self.cross_mass)
        beta_pep_list_b_ion, beta_pep_list_y_ion = tool_generate_only_cross_by_ion(beta_pep_aa_mass_list, ion_charge,
                                                                                   beta_pep_cross_site,
                                                                                   alpha_pep_data.mass + self.cross_mass)
        timer.elapsed_and_reset("gen pep ion")
        alpha_ion, alpha_ion_same_num = [], []  # 将b y离子合并
        tool_set_two_list(sorted(alpha_pep_list_b_ion), sorted(alpha_pep_list_y_ion), 0.0005, alpha_ion,
                          alpha_ion_same_num)

        # Logger.debug("beta_pep_list_b_ion: " + str(beta_pep_list_b_ion))
        beta_ion, beta_ion_same_num = [], []  # 将b y离子合并
        tool_set_two_list(sorted(beta_pep_list_b_ion), sorted(beta_pep_list_y_ion), 0.0005, beta_ion, beta_ion_same_num)

        all_ion, all_ion_same_num = [], []
        delta_alpha_ion, delta_beta_ion = [], []
        tool_set_two_list(alpha_ion, beta_ion, 0.0005, all_ion, all_ion_same_num, alpha_ion_same_num, beta_ion_same_num,
                          delta_alpha_ion, delta_beta_ion)
        timer.elapsed_and_reset("merge pep ion")
        alpha_continue_score, alpha_cover = [], []
        alpha_continue_data = op_get_peptide_by_continue_data(alpha_pep_list_b_ion, alpha_pep_list_y_ion,
                                                              alpha_pep_data.sq_length - 1, one_spectrum,
                                                              self.dp.myCFG.is_ppm_fra,
                                                              self.dp.myCFG.C2_PPM_TOL_FRAGMENT, 1)
        op_get_peptide_continue_score(alpha_continue_data, alpha_continue_score, alpha_cover)

        beta_continue_score, beta_cover = [], []
        beta_continue_data = op_get_peptide_by_continue_data(beta_pep_list_b_ion, beta_pep_list_y_ion,
                                                             beta_pep_data.sq_length - 1, one_spectrum,
                                                             self.dp.myCFG.is_ppm_fra,
                                                             self.dp.myCFG.C2_PPM_TOL_FRAGMENT, 1)
        op_get_peptide_continue_score(beta_continue_data, beta_continue_score, beta_cover)
        timer.elapsed_and_reset("get continue score")

        alpha_score = op_get_match_ion_score(one_spectrum, alpha_ion, alpha_ion_same_num, self.dp.myCFG.is_ppm_fra,
                                             self.dp.myCFG.C2_PPM_TOL_FRAGMENT, 1)
        beta_score = op_get_match_ion_score(one_spectrum, beta_ion, beta_ion_same_num, self.dp.myCFG.is_ppm_fra,
                                            self.dp.myCFG.C2_PPM_TOL_FRAGMENT, 1)
        delta_beta_ion_same_num = [1 for i in range(len(delta_beta_ion))]
        delta_beta_ion_score = op_get_match_ion_score(one_spectrum, delta_beta_ion, delta_beta_ion_same_num,
                                                      self.dp.myCFG.is_ppm_fra, self.dp.myCFG.C2_PPM_TOL_FRAGMENT, 1)
        all_score = op_get_match_ion_score(one_spectrum, all_ion, all_ion_same_num, self.dp.myCFG.is_ppm_fra,
                                           self.dp.myCFG.C2_PPM_TOL_FRAGMENT, 1)
        timer.elapsed_and_reset("get match score")
        # 有效离子肽谱匹配得分
        if len(all_ion) == 0:
            pep_part_1_score = 0.0
        else:
            pep_part_1_score = all_score.match_ion_score * len(all_score.list_ppm) / len(all_ion)

        # 连续性分数
        pep_part_2_score = 0.0
        pep_part_3_score = 0.0
        for i in range(len(alpha_continue_score)):
            pep_part_2_score += alpha_continue_score[i]
        for i in range(len(beta_continue_score)):
            pep_part_3_score += beta_continue_score[i]
        # 连续性分数

        # 肽段覆盖率
        if len(alpha_ion) == 0:
            pep_part_4_score = 0.0
        else:
            pep_part_4_score = len(alpha_score.list_ppm) / len(alpha_ion) * alpha_score.match_ion_score
        if len(beta_ion) == 0:
            pep_part_5_score = 0.0
        else:
            pep_part_5_score = len(beta_score.list_ppm) / len(beta_ion) * beta_score.match_ion_score

        if len(alpha_score.list_ppm) == 0:
            alpha_pep_score = 0.0
            alpha_percent = 0.0
        else:
            alpha_percent = len(alpha_score.list_ppm) / len(alpha_ion)
            alpha_pep_score = alpha_score.match_ion_score * alpha_percent

        if len(beta_score.list_ppm) == 0:
            beta_pep_score = 0.0
            beta_pep_percent = 0.0
        else:
            beta_pep_percent = len(alpha_score.list_ppm) / len(alpha_ion)
            beta_pep_score = beta_score.match_ion_score * beta_pep_percent

        pep_score = pep_part_1_score + pep_part_2_score + pep_part_3_score + pep_part_4_score + pep_part_5_score
        timer.elapsed_and_reset("get all score")
        cross_link_pep_mass = alpha_pep_data.mass + self.dp.myLINK.linker_1_data.cross_mass + beta_pep_data.mass + self.dp.myINI.MASS_PROTON_MONO
        precursor_mass_erro_Da = one_spectrum.mass - cross_link_pep_mass
        precursor_mass_erro_ppm = precursor_mass_erro_Da / cross_link_pep_mass * 1e6

        if len(all_score.list_ppm) < 3 or len(alpha_score.list_ppm) < 1 or len(beta_score.list_ppm) < 1:
            pass
        else:

            ppm_sum = np.sum(np.abs(np.array(all_score.list_ppm)))
            ppm_average = np.mean(np.abs(np.array(all_score.list_ppm)))
            # 求方差
            ppm_var = np.var(np.abs(np.array(all_score.list_ppm)))

            try:
                pParseNum = int(one_spectrum.title.split('.')[4])
            except:
                pParseNum = 0

            if one_spectrum.mobility is None:
                mobility = 0
            else:
                mobility = one_spectrum.mobility

            if alpha_score.match_ion_score > 0 and beta_score.match_ion_score > 0:
                if alpha_score.match_ion_score > beta_score.match_ion_score:

                    crosslink_delta_score = delta_beta_ion_score.match_ion_score
                    peptide_feature = CRerankTwoPeptideFeature(0,
                                                               pep_score, alpha_pep_score, beta_pep_score,
                                                               ppm_sum, ppm_average, ppm_var,
                                                               all_score,
                                                               pep_part_2_score, pep_part_3_score,
                                                               one_spectrum,
                                                               alpha_pep_data.sq_length, beta_pep_data.sq_length,
                                                               precursor_mass_erro_Da, crosslink_delta_score,
                                                               mobility,
                                                               0, pParseNum, 0)
                    if self.dp.myCFG.is_ppm_pre:
                        if abs(precursor_mass_erro_ppm) < self.dp.myCFG.C1_PPM_TOL_PRECURSOR * 1e6:
                            one_COnlyCrossResult = COnlyCrossResult(one_spectrum, alpha_pep_data, beta_pep_data,
                                                                    peptide_feature, alpha_pep_cross_site,
                                                                    beta_pep_cross_site)
                            op_fill_COnlyCrossResult(one_COnlyCrossResult, cross_link_pep_mass,
                                                     pep_score, alpha_pep_score, beta_pep_score,
                                                     precursor_mass_erro_Da,
                                                     precursor_mass_erro_ppm,
                                                     alpha_score, beta_score,
                                                     [all_score,
                                                      pep_part_1_score, pep_part_2_score, pep_part_3_score,
                                                      pep_part_4_score, pep_part_5_score]
                                                     )
                            return one_COnlyCrossResult
                            # all_fine_score.append(one_COnlyCrossResult)
                    else:
                        if abs(precursor_mass_erro_Da) < self.dp.myCFG.C1_PPM_TOL_PRECURSOR:
                            one_COnlyCrossResult = COnlyCrossResult(one_spectrum, alpha_pep_data, beta_pep_data,
                                                                    peptide_feature, alpha_pep_cross_site,
                                                                    beta_pep_cross_site)
                            op_fill_COnlyCrossResult(one_COnlyCrossResult, cross_link_pep_mass,
                                                     pep_score, alpha_pep_score, beta_pep_score,
                                                     precursor_mass_erro_Da,
                                                     precursor_mass_erro_ppm,
                                                     alpha_score, beta_score,
                                                     [all_score,
                                                      pep_part_1_score, pep_part_2_score, pep_part_3_score,
                                                      pep_part_4_score, pep_part_5_score]
                                                     )
                            return one_COnlyCrossResult
                            # all_fine_score.append(one_COnlyCrossResult)
                else:
                    crosslink_delta_score = delta_beta_ion_score.match_ion_score

                    peptide_feature = CRerankTwoPeptideFeature(0,
                                                               pep_score, beta_pep_score, alpha_pep_score,
                                                               ppm_sum, ppm_average, ppm_var,
                                                               all_score,
                                                               pep_part_3_score, pep_part_2_score,
                                                               one_spectrum,
                                                               beta_pep_data.sq_length, alpha_pep_data.sq_length,
                                                               precursor_mass_erro_Da, crosslink_delta_score,
                                                               mobility,
                                                               0, pParseNum, 0)
                    if self.dp.myCFG.is_ppm_pre:
                        if abs(precursor_mass_erro_ppm) < self.dp.myCFG.C1_PPM_TOL_PRECURSOR * 1e6:
                            one_COnlyCrossResult = COnlyCrossResult(one_spectrum, beta_pep_data, alpha_pep_data,
                                                                    peptide_feature, beta_pep_cross_site,
                                                                    alpha_pep_cross_site)
                            op_fill_COnlyCrossResult(one_COnlyCrossResult, cross_link_pep_mass,
                                                     pep_score, beta_pep_score, alpha_pep_score,
                                                     precursor_mass_erro_Da,
                                                     precursor_mass_erro_ppm,
                                                     beta_score, alpha_score,
                                                     [all_score,
                                                      pep_part_1_score, pep_part_3_score, pep_part_2_score,
                                                      pep_part_4_score, pep_part_5_score]
                                                     )
                            return one_COnlyCrossResult
                            # all_fine_score.append(one_COnlyCrossResult)
                    else:
                        if abs(precursor_mass_erro_Da) < self.dp.myCFG.C1_PPM_TOL_PRECURSOR:
                            one_COnlyCrossResult = COnlyCrossResult(one_spectrum, beta_pep_data, alpha_pep_data,
                                                                    peptide_feature, beta_pep_cross_site,
                                                                    alpha_pep_cross_site)
                            op_fill_COnlyCrossResult(one_COnlyCrossResult, cross_link_pep_mass,
                                                     pep_score, beta_pep_score, alpha_pep_score,
                                                     precursor_mass_erro_Da,
                                                     precursor_mass_erro_ppm,
                                                     beta_score, alpha_score,
                                                     [all_score,
                                                      pep_part_1_score, pep_part_3_score, pep_part_2_score,
                                                      pep_part_4_score, pep_part_5_score]
                                                     )
                            return one_COnlyCrossResult
                            # all_fine_score.append(one_COnlyCrossResult)
        timer.elapsed_and_reset("get res")

    def __captain_get_final_result_index(self, candidate_spectrum_result):
        # 获取最高分结果和第二名结果（这个是序列不一样的）
        # 先确定最高分
        max_crosspeptide_score = candidate_spectrum_result[0].pep_score

        score_index = 1
        second_index = -1
        out_result_index = 0
        out_result_flag = False
        # 这事为了处理target和decoy同分的情况
        while True:
            if (score_index == len(candidate_spectrum_result)):
                out_result_index = 0
                break
            if (max_crosspeptide_score - candidate_spectrum_result[score_index].pep_score > 1e-4):
                second_index = score_index
                break
            if candidate_spectrum_result[score_index].alpha_pep_data.is_target + candidate_spectrum_result[
                score_index].beta_pep_data.is_target > 1:
                if not out_result_flag:
                    out_result_index = score_index
                    out_result_flag = True
            score_index += 1
        return out_result_index, second_index
