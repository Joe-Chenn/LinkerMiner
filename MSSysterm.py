# -*- mode: python ; coding: utf-8 -*-..
# from MSLogging import LogLevel

AA_NUM = 26
VALUE_ILLEGAL = -7.16
VALUE_MAX_SCAN = 2000000

SITE_DIC = {'A': 0, 'B': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0,
            'H': 0, 'I': 0, 'J': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0,
            'O': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'U': 0,
            'V': 0, 'W': 0, 'X': 0, 'Y': 0, 'Z': 0,
            -2: 0, -1: 0, 1: 0, 2: 0}

DLL_FUNCTION = "C_eLink_function.dll"
DLL_INDEX = "C_eLink_index.dll"
DLL_CROSS = "C_eLink_cross.dll"
DLL_SINGLE = "C_eLink_single.dll"

PEPTIDE_INDEX_DATA_PKL_FILE = "peptide_index.pkl"

PRO_PKL_FILE = "protein.pkl"
PEP_PRO_INDEX_FILE = "peptide2protein.pkl"
MOD_PKL_FILE = "modification.pkl"
VRA_PKL_FILE = 'var_peptide.pkl'

EXPIRATION_TIME = {'Year': 2025, 'Month': 12, 'Day': 31}
# 参数文件
FILE_NAME_CONFIG = 'eLink_config.txt'
# 文件输出的后缀
FILE_TYPE = {1: '.txt', 2: '.pkl'}
FOLDER_MASS = 'mass_index'
# 新建文件夹存索引
FOLDER_INDEX = [{},
                {1: 'single_peptide', 2: 'only_loop'},
                {1: 'only_cross', 2: 'cross_cross', 3: 'cross_loop'}]
# 这里是索引文件名
FILE_NAME_RECORED_ION_INDEX = [{},
                               {1: 'record_ion_index_single_peptide', 2: 'record_ion_index_only_loop'},
                               {1: 'record_ion_index_only_cross', 2: 'record_ion_index_result_cross_cross',
                                3: 'record_ion_index_result_cross_loop'}]
# 离子索引
FILE_NAME_ION_INDEX = [{},
                       {1: 'ion_index_single_peptide', 2: 'ion_index_only_loop'},
                       {1: 'ion_index_only_cross', 2: 'ion_index_result_cross_cross', 3: 'ion_index_result_cross_loop'}]
# 存排好序的肽段编码
FILE_NAME_PEP_DATA = [{},
                      {1: 'pep_data_single_peptide', 2: 'pep_data_only_loop'},
                      {1: 'pep_data_only_cross', 2: 'pep_data_cross_cross', 3: 'pep_data_cross_loop'}]
# 存该肽段有多少个碎片离子
FILE_NAME_ION_NUM = [{},
                     {1: 'ion_num_single_peptide', 2: 'ion_num_only_loop'},
                     {1: 'ion_num_only_cross', 2: 'ion_num_cross_cross', 3: 'ion_num_cross_loop'}]
# 存肽段所属的肽段质量
FILE_NAME_PEP_MASS = [{},
                      {1: 'pep_mass_single_peptide', 2: 'pep_mass_only_loop'},
                      {1: 'pep_mass_only_cross', 2: 'pep_mass_cross_cross', 3: 'pep_mass_cross_loop'}]
# 存肽段所在的肽段编码位置
FILE_NAME_PEP_PKL = [{},
                     {1: 'pep_pkl_single_peptide', 2: 'pep_pkl_only_loop'},
                     {1: 'pep_pkl_only_cross', 2: 'pep_pkl_cross_cross', 3: 'pep_pkl_cross_loop'}]
# 交联位点,单肽不涉及所以没有
FILE_NAME_LINK_SITE = [{},
                       {2: 'link_site_only_loop'},
                       {1: 'link_site_only_cross', 2: 'link_site_cross_cross', 3: 'link_site_cross_loop'}]
# 这里输出的结果文件名
FILE_NAME_PSM_RESULT = [{},
                        {1: 'result_single_peptide', 2: 'result_only_loop'},
                        {1: 'result_only_cross', 2: 'complex_result_cross_cross', 3: 'complex_result_cross_loop'}]

# 重打分的
FILE_NAME_RERANK = [{},
                    {1: 'rerank_result_single_peptide', 2: 'rerank_result_only_loop'},
                    {1: 'rerank_result_only_cross', 2: 'rerank_complex_result_cross_cross',
                     3: 'rerank_complex_result_cross_loop'}]

TMP_FOLDER = '\\tmp\\'

# 和libsvm相关的
PATH_LIBSVM = '\\libsvm-3.24\\windows'
SCALE_EXE = 'svm-scale.exe'
TRAIN_EXE = 'svm-train.exe'
PREDICT_EXE = 'svm-predict.exe'

PPARSE_EXE = "pParse.exe"
TIMSTOFREADER_EXE = "timsTOFReader.exe"

TimsTOFReadConfigName = "timsTOFReader_config.txt"

MARK_LABEL_INFO = ('NONE', 'AA:', 'MOD:', 'CLAA1', 'CLAA2', 'LINK:', 'GLYCAN')

LOG_LEVEL = 3
