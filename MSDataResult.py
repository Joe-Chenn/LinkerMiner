# -*- mode: python ; coding: utf-8 -*-
class COnlyCrossCoarseData():

    def __init__(self):
        self.start_pos = None
        self.end_pos = None
        self.sq_length = None
        self.pro_index = None
        self.protein_list = None
        self.protein_start_pos_list = None

        self.mass = None

        self.mod_site_list = None
        self.pep2pro_is_C_term = None
        self.pep2pro_is_N_term = None
        self.is_target = None


class COnlyCrossResult():
    def __init__(self, one_spectrum_data, alpha_pep_data, beta_pep_data, peptide_feature, alpha_pep_cross_site,
                 beta_pep_cross_site):
        self.spectrum_title = one_spectrum_data.title
        self.spectrum_charge = one_spectrum_data.charge
        self.spectrum_mass = one_spectrum_data.mass
        self.mobility = one_spectrum_data.mobility

        self.alpha_pep_data = alpha_pep_data
        self.alpha_pep_cross_site = alpha_pep_cross_site
        self.beta_pep_data = beta_pep_data
        self.beta_pep_cross_site = beta_pep_cross_site
        self.pep_mass = 0.0
        self.pep_precursor_bias_Da = 0.0
        self.pep_precursor_bias_ppm = 0.0
        self.pep_score = 0.0
        self.alpha_pep_score = 0.0
        self.beta_pep_score = 0.0
        self.alpha_pep_match_data = []
        self.beta_pep_match_data = []
        self.feature = peptide_feature


class CWriteResultOnlyCross():
    def __init__(self, one_spectrum_data, FDR=0, rerank_score=0.0):
        self.spectrum_title = one_spectrum_data.title
        self.spectrum_charge = str(one_spectrum_data.charge)
        self.spectrum_mass = str(one_spectrum_data.mass)
        self.mobility = one_spectrum_data.mobility
        self.peptide_sq = ''
        self.peptide_type = ''
        self.linker = ''
        self.peptide_mass = ''
        self.peptide_modification = ''
        self.Evalue = str(0.0)
        self.peptide_score = 0.0
        self.precursor_mass_error_Da = ''
        self.precursor_mass_error_ppm = ''
        self.peptide_protein = ''
        self.protein_type = ''
        self.File_ID = ''
        self.Label_ID = ''
        self.alpha_matched = ''
        self.beta_matched = ''
        self.alpha_evalue = ''
        self.beta_evalue = ''
        self.feature = ''
        self.FDR = FDR
        self.rerank_score = rerank_score

        self.inter_protein = 0
        self.inter_FDR = 0
        self.intra_FDR = 0
        self.scan = int(one_spectrum_data.title.split('.')[2])


class SpectrumData:
    def __init__(self, title, charge, mass, mobility):
        self.title = title
        self.charge = charge
        self.mass = mass
        self.mobility = mobility

    # ================================

class CSinglePeptideCoarseData():

    def __init__(self):
        self.start_pos = None
        self.end_pos = None
        self.sq_length = None
        self.pro_index = None
        self.protein_list = None
        self.protein_start_pos_list = None
        self.mass = None
        self.mod_site_list = None
        self.pep2pro_is_C_term = None
        self.pep2pro_is_N_term = None
        self.is_target = None

        self.loop_link = False
        self.loop_site_1 = None
        self.loop_site_2 = None


class CSinglePeptideResult():

    def __init__(self, pep_data, peptide_feature, loop_link=False, loop_site_1=False, loop_site_2=False):
        self.pep_data = pep_data
        self.pep_mass = 0.0
        self.pep_precursor_bias_Da = 0.0
        self.pep_precursor_bias_ppm = 0.0
        self.pep_score = 0.0
        self.match_data = []
        self.feature = peptide_feature

        self.loop_link = loop_link
        self.loop_site_1 = loop_site_1
        self.loop_site_2 = loop_site_2


class CWriteResultSinglePeptide():

    def __init__(self, one_spectrum_data, loop_link=False, FDR=1, rerank_score=0.0):
        self.spectrum_title = one_spectrum_data.title
        self.spectrum_charge = str(one_spectrum_data.charge)
        self.spectrum_mass = str(one_spectrum_data.mass)
        self.mobility = one_spectrum_data.mobility
        self.peptide_sq = ''
        self.peptide_type = ''
        self.peptide_mass = ''
        self.peptide_modification = ''
        self.Evalue = str(0.0)
        self.peptide_score = 0
        self.precursor_mass_error_Da = ''
        self.precursor_mass_error_ppm = ''
        self.peptide_protein = ''
        self.protein_type = ''
        self.File_ID = ''
        self.Label_ID = ''
        self.pep_matched = ''
        self.pep_evalue = ''
        self.feature = ''
        self.FDR = FDR
        self.rerank_score = rerank_score

        self.scan = int(one_spectrum_data.title.split('.')[2])

        self.momo_link = False
        self.peptide_mono = ''

        self.loop_link = loop_link
        self.peptide_loop = ''
        self.linker_name = ''


# =================================

class CCrossCrossResult():

    def __init__(self, one_spectrum_data, alpha_pep_data, beta_pep_data, peptide_feature, alpha_pep_cross_site_data,
                 beta_pep_cross_site_data):
        self.spectrum_title = one_spectrum_data.title
        self.spectrum_charge = one_spectrum_data.charge
        self.spectrum_mass = one_spectrum_data.mass

        self.alpha_pep_data = alpha_pep_data
        self.alpha_pep_site_data = alpha_pep_cross_site_data
        self.beta_pep_data = beta_pep_data
        self.beta_pep_site_data = beta_pep_cross_site_data
        self.pep_mass = 0.0
        self.pep_precursor_bias_Da = 0.0
        self.pep_precursor_bias_ppm = 0.0
        self.pep_score = 0.0
        self.alpha_pep_score = 0.0
        self.beta_pep_score = 0.0
        self.alpha_pep_match_data = []
        self.beta_pep_match_data = []
        self.feature = peptide_feature


class CWriteResultCrossCross():

    def __init__(self, one_spectrum_data, FDR=0, rerank_score=0.0):
        self.spectrum_title = one_spectrum_data.title
        self.spectrum_charge = str(one_spectrum_data.charge)
        self.spectrum_mass = str(one_spectrum_data.mass)

        self.peptide_sq = ''
        self.peptide_type = ''
        self.linker = ''
        self.peptide_mass = ''
        self.peptide_modification = ''
        self.Evalue = str(0.0)
        self.peptide_score = 0
        self.precursor_mass_error_Da = ''
        self.precursor_mass_error_ppm = ''
        self.peptide_protein = ''
        self.protein_type = ''
        self.protein_type_math = None
        self.File_ID = ''
        self.Label_ID = ''
        self.alpha_matched = ''
        self.beta_matched = ''
        self.alpha_evalue = ''
        self.beta_evalue = ''
        self.feature = ''
        self.FDR = FDR
        self.rerank_score = rerank_score

        self.inter_protein = 0
        self.inter_FDR = 0
        self.intra_FDR = 0
        self.scan = int(one_spectrum_data.title.split('.')[2])

        '''
            str_alpha_mod = ''
            str_beta_mod = ''
            if len(final_spectrum_res.alpha_pep_mod) != 0:
                for i in range(len(final_spectrum_res.alpha_pep_mod)):
                    str_alpha_mod += str(final_spectrum_res.alpha_pep_mod[i].site) + ':' + str(
                        final_spectrum_res.alpha_pep_mod[i].mod_name) + ';'
            if len(final_spectrum_res.beta_pep_mod) != 0:
                for i in range(len(final_spectrum_res.beta_pep_mod)):
                    str_beta_mod += str(final_spectrum_res.beta_pep_mod[i].site) + ':' + str(
                        final_spectrum_res.beta_pep_mod[i].mod_name) + ';'
            str_protein = ''
            for alpha_protein_name in final_spectrum_res.alpha_pep_from_protein:
                for beta_protein_name in final_spectrum_res.beta_pep_from_protein:
                    str_protein += alpha_protein_name + '-' + beta_protein_name + '\\'
            peptide_type = str(final_spectrum_res.alpha_is_target) + ',' + str(final_spectrum_res.beta_is_target)
            if final_spectrum_res.alpha_is_target + final_spectrum_res.beta_is_target == 2:
                peptide_type_math = 1
            elif final_spectrum_res.alpha_is_target + final_spectrum_res.beta_is_target == 0:
                peptide_type_math = -1
            else:
                peptide_type_math = 0
            final_spectrum_res.feature.delta_score = delta_score
            self.spectrum_title = one_spectrum_data.title
            self.spectrum_charge = str(one_spectrum_data.charge)
            self.spectrum_mass = '%.4f' % one_spectrum_data.mass
            self.peptide_sq = final_spectrum_res.alpha_pep_sq + '(linker_1:' + str(final_spectrum_res.alpha_pep_cross_site_linker_1) + ',linker_2:' + str(final_spectrum_res.alpha_pep_cross_site_linker_2) + ')-' + final_spectrum_res.beta_pep_sq + '(linker_1:' + str(final_spectrum_res.beta_pep_cross_site_linker_1) + ',linker_2:' + str(final_spectrum_res.beta_pep_cross_site_linker_2) + ')'
            self.peptide_type = 'cross_cross'
            self.linker = linker_1_name + ';' + linker_2_name
            self.peptide_mass = '%.4f' % final_spectrum_res.pep_mass
            self.peptide_modification = str_alpha_mod + '-' + str_beta_mod
            self.Evalue = str(0.0)
            self.peptide_score = final_spectrum_res.pep_score
            self.precursor_mass_error_Da = '%.4f' % final_spectrum_res.pep_precursor_bias_Da
            self.precursor_mass_error_ppm = '%.4f' % final_spectrum_res.pep_precursor_bias_ppm
            self.peptide_protein = str_protein
            # 这个是肽段类型
            self.peptide_type_math = peptide_type_math
            self.peptide_type = peptide_type

            self.File_ID = 'File_ID'
            self.Label_ID = 'Label_ID'
            self.alpha_matched = ''
            self.beta_matched = ''
            self.alpha_evalue = ''
            self.beta_evalue = ''
            self.feature = final_spectrum_res.feature
            self.FDR = FDR
            self.rerank_score = rerank_score
        '''


# =================================

class CCrossLoopResult():

    def __init__(self, one_spectrum_data, alpha_pep_data, beta_pep_data, peptide_feature, alpha_pep_site_data,
                 beta_pep_site_data):
        self.spectrum_title = one_spectrum_data.title
        self.spectrum_charge = one_spectrum_data.charge
        self.spectrum_mass = one_spectrum_data.mass

        self.alpha_pep_data = alpha_pep_data
        self.alpha_pep_site_data = alpha_pep_site_data
        self.beta_pep_data = beta_pep_data
        self.beta_pep_site_data = beta_pep_site_data
        self.pep_mass = 0.0
        self.pep_precursor_bias_Da = 0.0
        self.pep_precursor_bias_ppm = 0.0
        self.pep_score = 0.0
        self.alpha_pep_score = 0.0
        self.beta_pep_score = 0.0
        self.alpha_pep_match_data = []
        self.beta_pep_match_data = []
        self.feature = peptide_feature


class CWriteResultCrossLoop():

    def __init__(self, one_spectrum_data, FDR=0, rerank_score=0.0):
        self.spectrum_title = one_spectrum_data.title
        self.spectrum_charge = str(one_spectrum_data.charge)
        self.spectrum_mass = str(one_spectrum_data.mass)
        self.mobility = one_spectrum_data.mobility

        self.peptide_sq = ''
        self.peptide_type = ''
        self.linker = ''
        self.peptide_mass = ''
        self.peptide_modification = ''
        self.Evalue = str(0.0)
        self.peptide_score = 0
        self.precursor_mass_error_Da = ''
        self.precursor_mass_error_ppm = ''
        self.peptide_protein = ''
        self.protein_type = ''
        self.protein_type_math = None
        self.File_ID = ''
        self.Label_ID = ''
        self.alpha_matched = ''
        self.beta_matched = ''
        self.alpha_evalue = ''
        self.beta_evalue = ''
        self.feature = ''
        self.FDR = FDR
        self.rerank_score = rerank_score

        self.inter_protein = 0
        self.inter_FDR = 0
        self.intra_FDR = 0
        self.scan = int(one_spectrum_data.title.split('.')[2])

        '''
            str_alpha_mod = ''
            str_beta_mod = ''
            if len(final_spectrum_res.alpha_pep_mod) != 0:
                for i in range(len(final_spectrum_res.alpha_pep_mod)):
                    str_alpha_mod += str(final_spectrum_res.alpha_pep_mod[i].site) + ':' + str(
                        final_spectrum_res.alpha_pep_mod[i].mod_name) + ';'
            if len(final_spectrum_res.beta_pep_mod) != 0:
                for i in range(len(final_spectrum_res.beta_pep_mod)):
                    str_beta_mod += str(final_spectrum_res.beta_pep_mod[i].site) + ':' + str(
                        final_spectrum_res.beta_pep_mod[i].mod_name) + ';'
            str_protein = ''
            for alpha_protein_name in final_spectrum_res.alpha_pep_from_protein:
                for beta_protein_name in final_spectrum_res.beta_pep_from_protein:
                    str_protein += alpha_protein_name + '-' + beta_protein_name + '\\'
            peptide_type = str(final_spectrum_res.alpha_is_target) + ',' + str(final_spectrum_res.beta_is_target)
            if final_spectrum_res.alpha_is_target + final_spectrum_res.beta_is_target == 2:
                peptide_type_math = 1
            elif final_spectrum_res.alpha_is_target + final_spectrum_res.beta_is_target == 0:
                peptide_type_math = -1
            else:
                peptide_type_math = 0
            final_spectrum_res.feature.delta_score = delta_score
            self.spectrum_title = one_spectrum_data.title
            self.spectrum_charge = str(one_spectrum_data.charge)
            self.spectrum_mass = '%.4f' % one_spectrum_data.mass
            self.peptide_sq = final_spectrum_res.alpha_pep_sq + '(linker_1:' + str(final_spectrum_res.alpha_pep_cross_site_linker_1) + ',linker_2:' + str(final_spectrum_res.alpha_pep_cross_site_linker_2) + ')-' + final_spectrum_res.beta_pep_sq + '(linker_1:' + str(final_spectrum_res.beta_pep_cross_site_linker_1) + ',linker_2:' + str(final_spectrum_res.beta_pep_cross_site_linker_2) + ')'
            self.peptide_type = 'cross_cross'
            self.linker = linker_1_name + ';' + linker_2_name
            self.peptide_mass = '%.4f' % final_spectrum_res.pep_mass
            self.peptide_modification = str_alpha_mod + '-' + str_beta_mod
            self.Evalue = str(0.0)
            self.peptide_score = final_spectrum_res.pep_score
            self.precursor_mass_error_Da = '%.4f' % final_spectrum_res.pep_precursor_bias_Da
            self.precursor_mass_error_ppm = '%.4f' % final_spectrum_res.pep_precursor_bias_ppm
            self.peptide_protein = str_protein
            # 这个是肽段类型
            self.peptide_type_math = peptide_type_math
            self.peptide_type = peptide_type

            self.File_ID = 'File_ID'
            self.Label_ID = 'Label_ID'
            self.alpha_matched = ''
            self.beta_matched = ''
            self.alpha_evalue = ''
            self.beta_evalue = ''
            self.feature = final_spectrum_res.feature
            self.FDR = FDR
            self.rerank_score = rerank_score
        '''


# =================================

class CRerankTwoPeptideFeature:
    def __init__(self, rerank_score,
                 match_score, match_alpha_score, match_beta_score,
                 match_error_sum, match_error_average, match_error_var,
                 CMATCH_ion_score,
                 continue_alpha_score, continue_beta_score,
                 spectrum_data,
                 alpha_pep_len, beta_pep_len,
                 precursor_bias, crosslink_delta_score,
                 mobility,
                 delta_score=0, pParseNum=0, FDR=0):
        # 分数类的特征
        self.rerank_score = rerank_score
        self.match_score = match_score
        self.match_alpha_score = match_alpha_score
        self.match_beta_score = match_beta_score
        # 匹配碎片离子偏差类特征
        self.match_error_sum = match_error_sum
        self.match_error_average = match_error_average
        self.match_error_var = match_error_var
        # 匹配上数目和总强度类特征
        self.match_ion_score = CMATCH_ion_score.match_ion_score
        self.match_ion_num = CMATCH_ion_score.match_ion_num
        self.match_ion_intensity = CMATCH_ion_score.match_inten_sum
        # 匹配上离子的百分比特征
        self.match_ion_num_percent = CMATCH_ion_score.match_ion_num_percent
        self.match_ion_intensity_percent = CMATCH_ion_score.match_ion_inten_percent
        # 匹配对于谱图的百分比特征
        self.match_spe_ion_percent = CMATCH_ion_score.match_spe_ion_percent
        self.match_spe_inten_percent = CMATCH_ion_score.match_spe_inten_percent
        # 连续性的特征
        self.continue_alpha_score = continue_alpha_score
        self.continue_beta_score = continue_beta_score
        # 其他和分数以及匹配情况不太相关的特征

        # 谱图的特征
        if spectrum_data.max_int == 0 or len(spectrum_data.peaks) == 0:
            self.spectrum_average_intensity = 0.0
        else:
            self.spectrum_average_intensity = spectrum_data.all_int / spectrum_data.max_int / len(spectrum_data.peaks)

        # 肽段的特征
        self.alpha_pep_len = alpha_pep_len
        self.beta_pep_len = beta_pep_len

        # 和肽段母离子相关的
        self.precursor_bias = precursor_bias
        self.delta_score = delta_score
        self.pParseNum = pParseNum
        self.mobility = mobility
        # 两条肽段差别的分数
        self.crosslink_delta_score = crosslink_delta_score

        self.FDR = FDR


class CRerankOnePeptideFeature:

    def __init__(self, rerank_score,
                 match_score,
                 match_error_sum, match_error_average, match_error_var,
                 CMATCH_ion_score,
                 continue_peptide_score,
                 precursor_bias, delta_score=0,
                 pParseNum=0, FDR=0):
        # 分数类的特征
        self.rerank_score = rerank_score
        self.match_score = match_score
        # 匹配碎片离子偏差类特征
        self.match_error_sum = match_error_sum
        self.match_error_average = match_error_average
        self.match_error_var = match_error_var
        # 匹配上数目和总强度类特征
        self.match_ion_num = CMATCH_ion_score.match_ion_num
        self.match_ion_intensity = CMATCH_ion_score.match_inten_sum
        # 匹配上离子的百分比特征
        self.match_ion_num_percent = CMATCH_ion_score.match_ion_num_percent
        self.match_ion_intensity_percent = CMATCH_ion_score.match_ion_inten_percent
        # 匹配对于谱图的百分比特征
        self.match_spe_ion_percent = CMATCH_ion_score.match_spe_ion_percent
        self.match_spe_inten_percent = CMATCH_ion_score.match_spe_inten_percent
        # 连续性的特征
        self.continue_peptide_score = continue_peptide_score
        # 其他和分数以及匹配情况不太相关的特征
        self.precursor_bias = precursor_bias
        self.delta_score = delta_score
        self.pParseNum = pParseNum
        self.FDR = FDR
