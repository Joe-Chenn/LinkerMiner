# -*- mode: python ; coding: utf-8 -*-..
from MSData import CDataPack
from MSDataResult import CWriteResultOnlyCross, SpectrumData
from MSFunctionFilterTG import CFunctionReadMS1, CFunctionFilterTG
from MSLogging import logGetError, Logger
from MSSysterm import PRO_PKL_FILE, PEPTIDE_INDEX_DATA_PKL_FILE, FILE_NAME_PSM_RESULT
from MSFunction import CFunctionLoadINIFile, CFunctionFillConfigData, CFunctionPickle, CFunctionSpectrum, \
    CFunctionCheck, CFunctionExtractMS2
from MSFunctionCreateIndex import CFunctionCreatePepIndex, CFunctionCreateOnePeptideIonIndex, \
    CFunctionCreateTwoPeptideIonIndex
from MSFunctionSearchTwoPeptide import CFunctionTwoPeptideSearch
from MSFunctionSearchOnePeptide import CFunctionOnePeptideSearch
from MSFunctionRerank import CFunctionRerankOnePeptide, CFunctionRerankTwoPeptide
from MSFunctionResult import CFunctionReadReasult
from MSOperator import tool_get_peptide_info

import os
import time


class CTaskCheck:

    # 有很多东西可以check：内、外存、文件等等；

    def __init__(self, inputDP):
        self.dp = inputDP

    def work(self):
        # 检查路径
        # 检查关键文件是否被打开，确保能写进去
        functionCheck = CFunctionCheck(self.dp)
        functionCheck.check_file()


class CTaskFillConfig:

    def __init__(self, inputDP):

        self.dp = inputDP

    def work(self):

        # 下面都是从这两个function调
        functionLoadINIFile = CFunctionLoadINIFile(self.dp)
        functionFillConfigData = CFunctionFillConfigData(self.dp)
        # 因为不同流程需要的参数不同，所以在这里进行分类，只加载这个流程需要用到的参数
        if self.dp.myCFG.D13_TYPE_FLOW in [1]:

            # 这两个一个是读INI文件，一个是读fasta文件并存成pkl
            functionLoadINIFile.get_INI_data()
            functionLoadINIFile.get_protein_data()
            # 只跑其中一个，靠D14_TYPE_PEPTIDE和D15_TYPE_LINK控制
            # 加载需要用到的参数
            functionFillConfigData.generate_mod_config()
            functionFillConfigData.generate_link_config()
            if self.dp.myCFG.D14_TYPE_PEPTIDE in [1] and self.dp.myCFG.D15_TYPE_LINK in [1]:
                pass
            else:
                # 需要把交联的mono-link加上
                functionFillConfigData.add_mono_var_mod(self.dp.myCFG.B7_LINKER_1_NAME)

            functionFillConfigData.generate_enzyme_config()
            functionFillConfigData.generate_G_matrix()
            functionFillConfigData.generate_massindex_config()
            functionFillConfigData.generate_ms2_config()

        elif self.dp.myCFG.D13_TYPE_FLOW in [2, 3]:
            # 这是跑常规的，单肽、mono-link、loop-link、cross-link
            # 这两个一个是读INI文件，一个是读fasta文件并存成pkl

            functionLoadINIFile.get_INI_data()
            functionLoadINIFile.get_protein_data()
            # 加载需要用到的参数
            functionFillConfigData.generate_mod_config()
            functionFillConfigData.generate_link_config()

            # 需要把交联的mono-link加上
            functionFillConfigData.add_mono_var_mod(self.dp.myCFG.B7_LINKER_1_NAME)

            functionFillConfigData.generate_enzyme_config()
            functionFillConfigData.generate_G_matrix()
            functionFillConfigData.generate_massindex_config()
            functionFillConfigData.generate_ms2_config()
        else:

            logGetError("[Info]The Flow may not exiting !")


class CTaskCreatePeptideIndex:

    def __init__(self, inputDP):
        self.dp = inputDP

    def work(self):
        functionPickle = CFunctionPickle()
        functionCreatePepIndex = CFunctionCreatePepIndex(self.dp)

        path_protein_pkl_file = os.path.join(self.dp.myCFG.A4_PATH_FASTA_EXPORT, PRO_PKL_FILE)
        protein_list = functionPickle.load_pkl_to_data(path_protein_pkl_file)
        # 产生肽段存储的文件
        pkl_file_list = []
        functionCreatePepIndex.generate_peptide_index_file(pkl_file_list)

        # 产生候选肽段
        peptide_matrix = []
        peptide_ind_matrix = []
        functionCreatePepIndex.create_peptide_index(protein_list, peptide_matrix, peptide_ind_matrix)

        functionCreatePepIndex.output_peptide_data(peptide_matrix, peptide_ind_matrix)
        del peptide_ind_matrix
        del peptide_matrix


class CTaskCreateIonIndex:

    def __init__(self, inputDP):

        self.dp = inputDP

    def work(self):

        functionPickle = CFunctionPickle()

        path_protein_pkl_file = os.path.join(self.dp.myCFG.A4_PATH_FASTA_EXPORT, PRO_PKL_FILE)
        protein_list = functionPickle.load_pkl_to_data(path_protein_pkl_file)

        # 生成质量区间pkl文件

        if self.dp.myCFG.D14_TYPE_PEPTIDE in [2]:
            functionCreateIonIndex = CFunctionCreateTwoPeptideIonIndex(self.dp)
            functionCreateIonIndex.create_ion_index(protein_list)

        if self.dp.myCFG.D14_TYPE_PEPTIDE in [1]:
            functionCreateIonIndex = CFunctionCreateOnePeptideIonIndex(self.dp)
            functionCreateIonIndex.create_ion_index(protein_list)


class CTaskPreProcessSpectrum:

    def __init__(self, inputDP):

        self.dp = inputDP

    def work(self):

        funcionExtractMS2 = CFunctionExtractMS2(self.dp)

        if self.dp.myCFG.A1_TYPE_MS2 in ['d']:
            # .d的需要先跑timsTOFReader
            for file_index in range(len(self.dp.myMS2.original_file_name)):
                funcionExtractMS2.extract_tims(file_index)

        elif self.dp.myCFG.A1_TYPE_MS2 in ['raw']:
            # raw的需要先跑pParse
            for file_index in range(len(self.dp.myMS2.original_file_name)):
                funcionExtractMS2.extract_raw(file_index)

        elif self.dp.myCFG.A1_TYPE_MS2 in ['mgf']:
            pass

        functionPreProcessCrossSpectrum = CFunctionSpectrum(self.dp)
        functionPreProcessCrossSpectrum.trans_ms2()


class CTaskSearch:

    def __init__(self, inputDP):

        self.dp = inputDP

    def work(self, id_psm_name=None, id_psm_scan=None):

        # 读入蛋白质序列
        functionPickle = CFunctionPickle()
        path_protein_pkl_file = os.path.join(self.dp.myCFG.A4_PATH_FASTA_EXPORT, PRO_PKL_FILE)
        protein_list = functionPickle.load_pkl_to_data(path_protein_pkl_file)

        # 不同流程需要加载的参数不一样，单肽不一样要单独来
        if self.dp.myCFG.D14_TYPE_PEPTIDE in [2]:

            functionTwoPeptideSearch = CFunctionTwoPeptideSearch(self.dp)

            if self.dp.myCFG.D15_TYPE_LINK in [1, 2, 3]:
                functionTwoPeptideSearch.update_class_variable(self.dp.myCFG.D15_TYPE_LINK)
                functionTwoPeptideSearch.search(id_psm_name, id_psm_scan)

        elif self.dp.myCFG.D14_TYPE_PEPTIDE in [1]:

            functionOnePeptideSearch = CFunctionOnePeptideSearch(self.dp)

            if self.dp.myCFG.D15_TYPE_LINK in [1, 2]:
                functionOnePeptideSearch.update_class_variable(self.dp.myCFG.D15_TYPE_LINK)
                functionOnePeptideSearch.search(protein_list, id_psm_name, id_psm_scan)

        else:

            pass


class CTaskRerank:

    def __init__(self, inputDP):

        self.dp = inputDP

    def work(self):

        if self.dp.myCFG.D14_TYPE_PEPTIDE in [2]:
            functionRerankTwopeptide = CFunctionRerankTwoPeptide(self.dp)
            functionRerankTwopeptide.updata_class_variable(0.05, 1)
            if self.dp.myCFG.D15_TYPE_LINK in [1]:
                functionRerankTwopeptide.rerank_onlycross()

        elif self.dp.myCFG.D14_TYPE_PEPTIDE in [1]:

            functionRerankOnepeptide = CFunctionRerankOnePeptide(self.dp)
            functionRerankOnepeptide.updata_class_variable(0.05, 1)
            if self.dp.myCFG.D15_TYPE_LINK in [1]:
                functionRerankOnepeptide.rerank_singlepeptide()

            if self.dp.myCFG.D15_TYPE_LINK in [2]:
                functionRerankOnepeptide.rerank_looplink()


class CTaskGetIdPSM:

    def __init__(self, inputDP):
        self.dp = inputDP

    def work(self, id_psm_name, id_psm_scan):
        functionReadReasult = CFunctionReadReasult(self.dp)
        functionReadReasult.get_id_spectrum(self.dp.myCFG.E1_FDR_PSM, id_psm_name, id_psm_scan)


class CTaskFilterTG:

    def __init__(self, inputDP):
        self.dp = inputDP

    def work(self):
        if self.dp.myCFG.D14_TYPE_PEPTIDE in [2]:
            taskReadMs1 = CTaskReadMS1(self.dp)
            taskReadMs1.work()

            if self.dp.myCFG.D13_TYPE_FLOW == 3:
                search_result = self.readResultFromPlink(
                    r"F:\syb_olsen\pLink_task_2024.04.09.19.58.05\reports\little_human_2024.04.09.filtered_cross-linked_spectra.csv")

                pass
            else:
                psm_result_file = self.dp.myCFG.A5_PATH_RESULT_EXPORT + "\\tmp\\" + FILE_NAME_PSM_RESULT[2][1] + '.pkl'
                functionPickle = CFunctionPickle()
                search_result = functionPickle.load_pkl_to_data(psm_result_file)
            functionFilterTG = CFunctionFilterTG(self.dp)
            functionFilterTG.filter(search_result)
        else:
            pass

    def readResultFromPlink(self, path):
        def convertMod(pepSeq: str, modifications: str):
            pep1Len = len(pepSeq.split("-")[0])
            pep2Len = len(pepSeq.split("-")[1])

            if modifications == "null":
                return "-"
            mods = modifications.split(";")

            converted = []
            for mod in mods:
                m = mod.split("(")[0]
                site = int(mod.split("(")[1][:-1])
                if site > pep1Len:
                    site -= pep1Len
                    site *= -1
                converted.append("{}:{}".format(site, m))
            return ";".join(converted)

        search_result = []
        with open(path, 'r') as f:
            header_line = next(f)
            headers = header_line.strip().split(',')
            title_index = headers.index("Title")
            charge_index = headers.index("Charge")
            peptide_index = headers.index("Peptide")
            linker_index = headers.index("Linker")
            modifications_index = headers.index("Modifications")
            proteins_index = headers.index("Proteins")
            protein_type_index = headers.index("Protein_Type")
            peptide_mass_index = headers.index("Peptide_Mass")
            for line in f:
                line = line.strip().split(',')
                spectra_data = SpectrumData(line[title_index],
                                            int(line[charge_index]),
                                            float(line[peptide_mass_index]),
                                            0)
                one_result = CWriteResultOnlyCross(spectra_data)
                one_result.peptide_sq = line[peptide_index]
                one_result.peptide_modification = convertMod(one_result.peptide_sq, line[modifications_index])
                one_result.peptide_protein = line[proteins_index].replace("/", "\\")
                one_result.peptide_type = "cross"
                one_result.spectrum_charge = line[charge_index]
                one_result.scan = line[title_index].split(".")[1]
                one_result.linker = line[linker_index]
                search_result.append(one_result)
        return search_result


class CTaskReadMS1:

    def __init__(self, inputDP: CDataPack):

        self.dp = inputDP

    def work(self):
        if self.dp.myCFG.A1_TYPE_MS2 in ['raw']:
            # 转化cfg的raw文件名到ms1文件名，继续读取ms1数据
            LIST_PATH_RAW = []
            for file_name in self.dp.myMS2.original_file_name:
                LIST_PATH_RAW.append(file_name)
            self.dp.LIST_PATH_MS1 = [i_raw.replace('.raw', '.ms1') for i_raw in LIST_PATH_RAW]
            for path in self.dp.LIST_PATH_MS1:
                print("Reading ms1: " + path)
                functionMS1 = CFunctionReadMS1(self.dp)
                functionMS1.read(path)
        # 提取路径中的文件名，去掉先导路径和.后缀
        self.dp.LIST_MS1_NAME = [os.path.basename(i_ms1).replace('.ms1', '') for i_ms1 in self.dp.LIST_PATH_MS1]
