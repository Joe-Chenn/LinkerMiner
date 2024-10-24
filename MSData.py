# -*- mode: python ; coding: utf-8 -*-
from MSSysterm import VALUE_ILLEGAL, VALUE_MAX_SCAN


class CConfig:
    # [data]
    A1_TYPE_MS2 = "mgf"
    A2_PATH_MS2 = " #VIP"
    A3_PATH_FASTA = " #VIP"
    A4_PATH_FASTA_EXPORT = "./protein_index/ #VIP"
    A5_PATH_RESULT_EXPORT = "./result/ #VIP"

    # [biology]
    B1_NAME_ENZYME = "trypsin KR C #use ';' to set multiply enzymes"
    B2_TYPE_DIGEST = "0 #0 for specific; 1 for semi-specific; 2 for non-specific"
    B3_NUMBER_MAX_MISS_CLV = "3"
    B4_NAME_MOD_FIX = "Carbamidomethyl[C] #VIP, use ';' to set multiply fixed modifications"
    B5_NAME_MOD_VAR = "Oxidation[M] #VIP, use ';' to set multiply variable modifications"
    B6_NUMBER_MAX_MOD = "3 #Maximum of variable modification in one peptide sequence (not consider the fixed modifications)"
    B7_LINKER_1_NAME = 'BS3'
    B8_LINKER_2_NAME = 'B3'
    B9_MAX_PEPNUM = '9000000'

    # [mass spectrometry]
    C1_PPM_TOL_PRECURSOR = "20ppm"
    C2_PPM_TOL_FRAGMENT = "20ppm"
    C3_activation_type = "HCD"

    # [performance]
    D1_NUMBER_THREAD = "8"
    D2_TYPE_THREAD = "0 #0 is for multi-process (high speed but use more memory); 1 is for multi-thread (low speed but use less memory)"
    D3_NUMBER_SELECT_PEAK = "200"
    D4_NUMBER_SPECTRUM = "200000"
    D5_LEN_MAX_PROTEIN = "100000"
    D6_MASS_PEP_LOW = "500"
    D7_MASS_PEP_UP = "6000"
    D8_LEN_PEP_LOW = "5"
    D9_LEN_PEP_UP = "60"
    D10_INDEX_SPLIT_MASS = "100 #create one pkl file for each 100Da ([0, 100], [100, 200], ..., [9900, 10000])"
    D11_NUMBER_TOP_RESULT = "10 #output top-10 peptides for each spectrum"
    D12_MULTI_MASS = "100 #use mass hash (mass*MUTLI_MASS) to retrive peptide, spectrum or peak, this value of creating peptide index and searching mgf must be same."
    D13_TYPE_FLOW = "1"
    D14_TYPE_PEPTIDE = '1'
    D15_TYPE_LINK = '1'
    D16_FILTER_LOOP = 0

    # [filter]
    E1_FDR_PSM = "0.05"
    E2_FDR_SEPARATE = '1'
    E3_MOBILITY_RERANK = '1'
    # [ini]
    F1_PATH_INI_ELEMENT = "./ini/element.ini"
    F2_PATH_INI_AA = "./ini/aa.ini"
    F3_PATH_INI_MOD = "./ini/modification.ini"
    F4_PATH_INI_XLINK = './ini/xlink.ini'

    is_ppm_fra = True

    using_decoy = True
    max_fragment_charge = 2
    max_mz = 6000
    max_mz_index_size = 600000
    top_N = 10

    typeC = 1
    single_rerank_time = 1
    rerank_time = 5
    cross_cross_rerank_time = 2
    separate_FDR = 1

    start_mass = 0
    end_mass = 0

    S1_SKIP_CREATE_PEPTIDE_INDEX = 0
    S2_SKIP_CREATE_ION_INDEX = 0
    S3_SKIP_SEARCH = 1
    S4_SKIP_RERANK = 1


class CINI:
    MASS_ELECTRON = 0.0005485799
    MASS_H2O = 18.0105647
    MASS_PROTON_MONO = 1.00727645224  # 1.00782503214-0.0005485799
    MASS_PROTON_ARVG = 1.0025

    MASS_CH2 = 14.0156492
    MASS_CO = 27.994914100000003
    MASS_NH = 15.010897799999999

    DIC_ELEMENT_MASS = {}
    DIC_ELEMENT_ALL_MASS = {}
    DIC_ELEMENT_MASS_PERCENTAGE = {}
    DIC_MOD = {}
    DIC_AA = {}
    DIC_LINKER = {}

    DIC_AA_COMPOSITION = {}



class CLinkerData:
    def __init__(self, linker_name, linker_site, dic_site, alpha_site, beta_site, same_site, cross_mass, loop_mass,
                 mono_mass, composition):
        self.linker_name = linker_name
        self.linker_site = linker_site
        self.dic_site = dic_site
        self.alpha_site = alpha_site
        self.beta_site = beta_site
        self.same_site = same_site
        self.cross_mass = cross_mass
        self.loop_mass = loop_mass
        self.mono_mass = mono_mass
        self.composition = composition


class CProtein:
    ac = []  # 名称
    de = []  # 描述
    sq = []  # 序列
    sq_original = []


class CMod:
    fix_mod_dic = {}
    fix_mod_dic_pep_N = {}
    fix_mod_dic_pep_C = {}
    fix_mod_dic_pro_N = {}
    fix_mod_dic_pro_C = {}

    var_mod_dic = {}
    var_mod_dic_pep_N = {}
    var_mod_dic_pep_C = {}
    var_mod_dic_pro_N = {}
    var_mod_dic_pro_C = {}

    fix_mod_num = 0
    var_mod_num = 0
    mono_mod_num = 0

    fixmodname2index = {}  # 修饰对应在mod_list中的位置
    fix_mod_list = []

    varmodname2index = {}
    var_mod_list = []  # 设置的修饰全在这个list中


class CModSite:

    def __init__(self, site, fix_or_var, mod_index):
        self.site = site  # the modification occured in the index of peptide sequence
        self.fix_or_var = fix_or_var  # 0是固定修饰，1是可变修饰
        self.mod_index = mod_index  # the mass of modification


class PeptideData:
    __slots__ = ["pep_code", "gdm", "mass", "fix_mod_code"]

    def __int__(self):
        self.pep_code = []
        self.gdm = []
        self.mass = []
        self.fix_mod_code = []


class CPeptidePKL:

    def __init__(self):
        self.pep_code = None
        self.gdm = None
        self.mass = None
        self.fix_mod_code = None  # 这个用64位二进制记录固定修饰，也就限制了肽段长度最长为62长
        self.var_mod_num = None  # 记录这个肽段上有多少个可变修饰，用来解下面这个list的值
        self.var_mod_code = None  # 这个是一个list存二进制数，其中用12位二进制记一个可变修饰，6位为肽段位置，6位为修饰索引号，
    # self.pep_code_list = None


class CPeak:
    __slots__ = ['mz', 'inten']

    def __init__(self, mz, inten):
        self.mz = mz
        self.inten = inten


class CSpectrum:
    __slots__ = ['title', 'mass', 'charge', 'max_int', 'all_int', 'peaks', 'peaks_index', 'mobility']

    def __init__(self, title, mass, charge, spe_max_int, spe_all_int, peaks_list, peaks_index, mobility=None):
        # the mass of precursor ion is added by (H2O+P)
        self.title = title
        self.mass = mass
        self.charge = charge
        self.max_int = spe_max_int
        self.all_int = spe_all_int
        self.peaks = peaks_list
        self.peaks_index = peaks_index
        self.mobility = mobility


class CLink:
    linker_1_data = None
    linker_2_data = None
    same_linker = 0


class CEnzyme:
    enzyme_name = []
    enzyme_aa = []
    enzyme_flag = []

    enzyme_C = False
    enzyme_N = False

    aa_enzyme_C = None
    aa_enzyme_N = None


class CMS2:
    ms2_num = 0
    original_file_name = []
    ms2_file = []
    trans_ms2_file = []
    trans_ms2_file_no_linker_ion = []
    trans_ms2_file_linker_ion = []


class CDataPack:
    # 全局维护的
    myCFG = CConfig()
    myINI = CINI()
    myMOD = CMod()
    myENZ = CEnzyme()
    myLINK = CLink()
    myMS2 = CMS2()

    G_matrix = []
    aa_num = 0
    max_sq_len = 0

    PROTEIN_DATA = None

    LIST_FW_MASS_INDEX = []
    LIST_MASS_INDEX_FILE = []
    LIST_MASS_IND_FILE = []
    LIST_PATH_MS1 = []
    LIST_MS1_NAME = []


class CModData:

    def __init__(self, name, mass, type, aa, composition=""):
        self.name = name  # modification name
        self.mass = mass  # modification mass
        self.type = type  # modification type (NORMAL, PEP_N, PEP_C, PRO_N, PRO_C)
        self.aa = aa  # modification can be occured in which amino acids
        self.composition = composition


class CMatchOnePeptideScore:

    def __init__(self, ion_continue_score, match_ion_score, match_inten_score, match_ion_spe_score, peptide_coverage,
                 list_ppm=[], match_ion_intensity=0.0, other_data=[]):
        self.ion_continue_score = ion_continue_score
        self.match_ion_score = match_ion_score
        self.match_inten_score = match_inten_score
        self.match_ion_spe_score = match_ion_spe_score
        self.peptide_coverage = peptide_coverage
        self.list_ppm = list_ppm
        self.other_data = other_data
        self.match_ion_num = len(list_ppm)
        self.match_ion_intensity = match_ion_intensity


class CMatchIonScore:

    def __init__(self, match_ion_score, list_Da, list_ppm, match_ion_num_percent, match_ion_inten_percent,
                 match_spe_ion_percent, match_spe_inten_percent, match_ion_num, match_inten_sum, spe_inten_sum,
                 other_data=None):
        self.match_ion_score = match_ion_score
        self.list_Da = list_Da
        self.list_ppm = list_ppm
        self.match_ion_num_percent = match_ion_num_percent
        self.match_ion_inten_percent = match_ion_inten_percent
        self.match_spe_ion_percent = match_spe_ion_percent
        self.match_spe_inten_percent = match_spe_inten_percent
        self.match_ion_num = match_ion_num
        self.match_inten_sum = match_inten_sum
        self.spe_inten_sum = spe_inten_sum
        self.other_data = other_data


class CFileMS2:  # 注意：这是按列存储，每个属性都是list。如果按列搞，这个行是一个对象，不太好管理。

    MS2_NUM = 0

    INDEX_SCAN = []  # 这个大小和大家不一样，方便快速索引
    INDEX_RT = []

    MATRIX_FILE_NAME = []
    LIST_PEAK_MOZ = []  # 每一行是个list
    LIST_PEAK_INT = []
    # 这个是不跨交联的
    MATRIX_PEAK_MOZ_NO_LINKER = []  # 每一行是个list，多个母离子的也是只存一个就好
    MATRIX_PEAK_INT_NO_LINKER = []
    # 这个是跨交联的， 因为母离子不同，所以一行是一张谱图，与MATRIX_PRECURSOR_CHARGE数目对应
    MATRIX_PEAK_MOZ_LINKER = []
    MATRIX_PEAK_INT_LINKER = []

    MATRIX_PRECURSOR_CHARGE = []  # 相同的scan，可能有多个母离子状态（质量+电荷）
    MATRIX_PRECURSOR_MOZ = []  # 相同的scan，可能有多个母离子状态（质量+电荷）
    MATRIX_MOBILITY = []

    LIST_RET_TIME = []  # 暂时与我无关
    LIST_ION_INJECTION_TIME = []  # 暂时与我无关
    LIST_ACTIVATION_CENTER = []  # 暂时与我无关
    LIST_PRECURSOR_SCAN = []


# RES_INDEX = # 这个是记录第0个母离子对应的谱图号


class CTimsTOFReaderConfig:
    # [data]
    A1_FOLDER_TIMS_DATA = ""
    A2_FOLDER_EXPORT = ""
    A3_MS_NAME = ""

    # [performance]
    B1_EXTRACT_MS1 = '1'
    B2_EXTRACT_MS2 = '1'

    B3_EXTRACT_MGF = '1'

    B4_MERGE_MS2 = '1'

    B5_MS1_INTENSITY_THRESHOLD = '100'
    B6_MS2_INTENSITY_THRESHOLD = '0'

    B7_MS1_FILTER = '0'
    B8_MS2_FILTER = '0'

    # [task]
    C1_TYPE_FLOW = '1'

    # 这些参数后续可以考虑加到参数文件中
    S1_PEAK_WINDOW_SIZE = 5  # 这个必须是整数因为用索引号处理


class CFileMS1:
    INDEX_SCAN = []  # 这两个大小和大家不一样，方便快速索引
    INDEX_RT = []

    # LIST_RET_TIME = [] 理论上应该有，但是一直没用过
    LIST_ION_INJECTION_TIME = []

    MATRIX_PEAK_MOZ = []  # 每一行是个list
    MATRIX_PEAK_INT = []

    @staticmethod
    def newInstance():
        instance = CFileMS1()
        instance.INDEX_SCAN = []
        instance.INDEX_RT = []

        instance.LIST_ION_INJECTION_TIME = [VALUE_ILLEGAL] * VALUE_MAX_SCAN

        instance.MATRIX_PEAK_MOZ = [[] * 1] * VALUE_MAX_SCAN  # 每一行是个list
        instance.MATRIX_PEAK_INT = [[] * 1] * VALUE_MAX_SCAN
        return instance


class CSeed:  # 为得到Evidence准备的，Labeling只需要这就可以了。LabelFree还需要个Reference。

    MID_SCAN = 0  # 根据需求不同，使用scan或rt
    MID_RT = 0

    DICT_COMPOSITION = {}  # 这个有时也需要

    DIS_ISO_MOZ_CLC = []  # 理论质量
    DIS_ISO_INT_CLC = []  # 理论强度

    VALUE_MOBILITY = 0.0

    INDEX_MONO = 0  # 上面几根峰，哪根是mono

    WET_INT = 1.0  # 强度的权重

    SCORE2 = 0.0


class CEvidence:
    # 这是第二批要填的信息
    MATRIX_PROFILE = []  # 每一行是个曲线，一样长

    # 这是第三批要填的信息
    LIST_RET_TIME = []  # 曲线上每个点的保留时间
    LIST_SCAN = []  # 曲线上每个点的scan

    LIST_I_START = []  # 每条曲线单独确定起止
    LIST_I_END = []

    I_START = 0  # 最后确定总的起止。在老pQuant里面，这一步很关键。所有profile要取最短的。
    I_END = -1

    # 这是第四批要填的信息
    DIS_ISO_MOZ_EXP = []  # 实验质量
    DIS_ISO_INT_EXP = []  # 实验强度
    PROFILE_ALL = []  # 一行；是所有同位素峰加和的曲线

    SCORE_IS_PEPTIDE = 1.0
