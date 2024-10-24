# -*- mode: python ; coding: utf-8 -*-..
from MSSysterm import VALUE_ILLEGAL, VALUE_MAX_SCAN
from MSTool_code import tool_get_peptide_info, tool_get_var_mod_data


def op_get_fix_mod_site_list(pep_sq, peptide_length, fix_mod_code, input_CMod, input_mod_site_list):
    for site in range(peptide_length + 2):
        site_data = fix_mod_code & 1
        if site_data == 0:
            input_mod_site_list.append(-1)
        else:
            if site == 0:
                if pep_sq[0] in input_CMod.fix_mod_dic_pro_N.keys():
                    input_mod_site_list.append(input_CMod.fix_mod_dic_pro_N[pep_sq[0]])
                else:
                    if pep_sq[0] in input_CMod.fix_mod_dic_pep_N.keys():
                        input_mod_site_list.append(input_CMod.fix_mod_dic_pep_N[pep_sq[0]])
            elif site == peptide_length + 1:
                if pep_sq[-1] in input_CMod.fix_mod_dic_pro_C.keys():
                    input_mod_site_list.append(input_CMod.fix_mod_dic_pro_C[pep_sq[-1]])
                else:
                    if pep_sq[-1] in input_CMod.fix_mod_dic_pep_C.keys():
                        input_mod_site_list.append(input_CMod.fix_mod_dic_pep_C[pep_sq[-1]])
                    else:
                        input_mod_site_list.append(-1)
            else:
                input_mod_site_list.append(input_CMod.fix_mod_dic[pep_sq[site - 1]])
        fix_mod_code = fix_mod_code >> 1


def op_get_fix_mod_code(input_mod_site_list):
    # 这个是为了和上面配套，就在operator里也加了这个
    fix_mod_code = 0
    for site, site_data in enumerate(input_mod_site_list):
        if site_data != -1:
            fix_mod_code = fix_mod_code | (1 << site)
    return fix_mod_code


def op_restore_fix_mod_site_list(input_mod_site_list, pep_sq, input_CMod):
    # 因为固定修饰是用二进制存，只能记录是否发生固定修饰，不能确定是什么修饰
    for i in range(len(input_mod_site_list)):
        if input_mod_site_list[i] == 0:
            if i == 0:
                if pep_sq[0] in input_CMod.fix_mod_dic_pro_N.keys():
                    input_mod_site_list[i] = input_CMod.fix_mod_dic_pro_N[pep_sq[0]]
                else:
                    if pep_sq[0] in input_CMod.fix_mod_dic_pep_N.keys():
                        input_mod_site_list[i] = input_CMod.fix_mod_dic_pep_N[pep_sq[0]]
            elif i == len(input_mod_site_list) - 1:
                if pep_sq[-1] in input_CMod.fix_mod_dic_pro_C.keys():
                    input_mod_site_list[i] = input_CMod.fix_mod_dic_pro_C[pep_sq[-1]]
                else:
                    if pep_sq[-1] in input_CMod.fix_mod_dic_pep_C.keys():
                        input_mod_site_list[i] = input_CMod.fix_mod_dic_pep_C[pep_sq[-1]]
            else:
                input_mod_site_list[i] = input_CMod.fix_mod_dic[pep_sq[i - 1]]


def op_modify_timsReaderconfig(inputConfig, path_d_file, path_out, ms_name):
    inputConfig.A1_FOLDER_TIMS_DATA = path_d_file
    inputConfig.A2_FOLDER_EXPORT = path_out
    inputConfig.A3_MS_NAME = ms_name


def op_INIT_CFILE_MS2(inputMS2):
    inputMS2.INDEX_SCAN = []
    inputMS2.INDEX_RT = []

    inputMS2.LIST_RET_TIME = [VALUE_ILLEGAL] * VALUE_MAX_SCAN
    inputMS2.LIST_ION_INJECTION_TIME = [VALUE_ILLEGAL] * VALUE_MAX_SCAN
    inputMS2.LIST_ACTIVATION_CENTER = [VALUE_ILLEGAL] * VALUE_MAX_SCAN

    inputMS2.LIST_PEAK_MOZ = [[] * 1] * VALUE_MAX_SCAN  # 每一行是个list
    inputMS2.LIST_PEAK_INT = [[] * 1] * VALUE_MAX_SCAN

    inputMS2.MATRIX_PEAK_MOZ_NO_LINKER = [[] * 1] * VALUE_MAX_SCAN  # 每一行是个list
    inputMS2.MATRIX_PEAK_INT_NO_LINKER = [[] * 1] * VALUE_MAX_SCAN

    inputMS2.MATRIX_PEAK_MOZ_LINKER = [[] * 1] * VALUE_MAX_SCAN  # 每一行是个list
    inputMS2.MATRIX_PEAK_INT_LINKER = [[] * 1] * VALUE_MAX_SCAN

    inputMS2.MATRIX_FILE_NAME = [[] * 1] * VALUE_MAX_SCAN  # 相同的scan，可能有多个母离子状态（质量+电荷）
    inputMS2.MATRIX_PRECURSOR_CHARGE = [[] * 1] * VALUE_MAX_SCAN  # 相同的scan，可能有多个母离子状态（质量+电荷）
    inputMS2.MATRIX_PRECURSOR_MOZ = [[] * 1] * VALUE_MAX_SCAN  # 相同的scan，可能有多个母离子状态（质量+电荷）
    inputMS2.MATRIX_MOBILITY = [[] * 1] * VALUE_MAX_SCAN


def op_copy_CPeptide(new_CPeptide, old_CPeptide):
    new_CPeptide.mass = old_CPeptide.mass
    new_CPeptide.gdm = old_CPeptide.gdm
    new_CPeptide.pep_code = old_CPeptide.pep_code
    new_CPeptide.fix_mod_code = old_CPeptide.fix_mod_code

    new_CPeptide.var_mod_code = []

    if old_CPeptide.var_mod_num is None:
        new_CPeptide.var_mod_num = 0
    else:
        new_CPeptide.var_mod_num = old_CPeptide.var_mod_num
        for item in old_CPeptide.var_mod_code:
            new_CPeptide.var_mod_code.append(item)

    # new_CPeptide.pep_code_list = []
    # if old_CPeptide.pep_code_list is None:
    #     new_CPeptide.pep_code_list = []
    # else:
    #     for item in old_CPeptide.pep_code_list:
    #         new_CPeptide.pep_code_list.append(item)


def op_fill_COnlyCrossCoarseData(inputCOnlyCoarseData, pep_info, protein_list, peptide2protein, CMod):
    is_target, protein_index, peptide_start, peptide_length = tool_get_peptide_info(pep_info[1])

    inputCOnlyCoarseData.is_target = is_target
    inputCOnlyCoarseData.pro_index = protein_index
    inputCOnlyCoarseData.start_pos = peptide_start
    inputCOnlyCoarseData.sq_length = peptide_length
    inputCOnlyCoarseData.end_pos = peptide_start + peptide_length

    inputCOnlyCoarseData.protein_list = []
    inputCOnlyCoarseData.protein_start_pos_list = []
    for other_pep_code in peptide2protein[pep_info[3]]:
        pep_data = tool_get_peptide_info(other_pep_code)
        inputCOnlyCoarseData.protein_list.append(pep_data[1])
        inputCOnlyCoarseData.protein_start_pos_list.append(pep_data[2])

    inputCOnlyCoarseData.mass = pep_info[0]

    mod_site_list = []
    pep_sq = protein_list.sq[protein_index][peptide_start: peptide_start + peptide_length]
    op_get_fix_mod_site_list(pep_sq, peptide_length, pep_info[2], CMod, mod_site_list)
    if len(pep_info) == 6:
        tool_get_var_mod_data(pep_info[5], mod_site_list)
    inputCOnlyCoarseData.mod_site_list = mod_site_list

    if peptide_start == 0:
        inputCOnlyCoarseData.pep2pro_is_N_term = True
    if peptide_start + peptide_length == len(protein_list.sq[protein_index]):
        inputCOnlyCoarseData.pep2pro_is_C_term = True


def op_fill_CSingleCoarseData(inputCSingleCoarseData, pep_info, protein_list, peptide2protein, CMod, loop_link=False,
                              loop_site_1=None, loop_site_2=None):
    is_target, protein_index, peptide_start, peptide_length = tool_get_peptide_info(pep_info[1])

    inputCSingleCoarseData.is_target = is_target
    inputCSingleCoarseData.pro_index = protein_index
    inputCSingleCoarseData.start_pos = peptide_start
    inputCSingleCoarseData.sq_length = peptide_length
    inputCSingleCoarseData.end_pos = peptide_start + peptide_length

    inputCSingleCoarseData.protein_list = []
    inputCSingleCoarseData.protein_start_pos_list = []
    for other_pep_code in peptide2protein[pep_info[3]]:
        pep_data = tool_get_peptide_info(other_pep_code)
        inputCSingleCoarseData.protein_list.append(pep_data[1])
        inputCSingleCoarseData.protein_start_pos_list.append(pep_data[2])

    inputCSingleCoarseData.mass = pep_info[0]

    mod_site_list = []
    pep_sq = protein_list.sq[protein_index][peptide_start: peptide_start + peptide_length]
    op_get_fix_mod_site_list(pep_sq, peptide_length, pep_info[2], CMod, mod_site_list)
    if len(pep_info) == 6:
        tool_get_var_mod_data(pep_info[5], mod_site_list)
    inputCSingleCoarseData.mod_site_list = mod_site_list

    if peptide_start == 0:
        inputCSingleCoarseData.pep2pro_is_N_term = True
    if peptide_start + peptide_length == len(protein_list.sq[protein_index]):
        inputCSingleCoarseData.pep2pro_is_C_term = True
    if loop_link:
        inputCSingleCoarseData.loop_link = True
        inputCSingleCoarseData.loop_site_1 = loop_site_1
        inputCSingleCoarseData.loop_site_2 = loop_site_2


def op_create_peaks_index(peak_list, multi):
    end_mass = peak_list[-1].mz
    num = int(end_mass * multi) + 1
    index_list = [-1 for i in range(num)]
    for i in range(len(peak_list)):
        p = peak_list[i]
        index = int(p.mz * multi)
        if index >= num: index = num - 1
        if index_list[index] == -1:
            index_list[index] = i
    end_val = len(peak_list)
    for i in range(num)[::-1]:
        if index_list[i] == -1:
            index_list[i] = end_val
        else:
            end_val = index_list[i]
    return index_list


# ===============================

def op_fill_COnlyCrossResult(input_COnlyCrossResult, mass, pep_score, alpha_pep_score, beta_pep_score,
                             pep_precursor_bias_Da, pep_precursor_bias_ppm, alpha_match_data=[], beta_match_data=[],
                             other_data=[]):
    input_COnlyCrossResult.pep_mass = mass
    input_COnlyCrossResult.pep_precursor_bias_Da = pep_precursor_bias_Da
    input_COnlyCrossResult.pep_precursor_bias_ppm = pep_precursor_bias_ppm
    input_COnlyCrossResult.pep_score = pep_score
    input_COnlyCrossResult.alpha_pep_score = alpha_pep_score
    input_COnlyCrossResult.beta_pep_score = beta_pep_score

    input_COnlyCrossResult.other_data = other_data
    input_COnlyCrossResult.alpha_pep_match_data = alpha_match_data
    input_COnlyCrossResult.beta_pep_match_data = beta_match_data


def op_fill_CWriteOnlyCrossResult(input_CWriteOnlyCrossResult, input_COnlyCrossResult, myMOD, protein_list, linker_name,
                                  FDR=0, rerank_score=0.0):
    str_alpha_mod = ''
    str_beta_mod = ''
    for site, mod_index in enumerate(input_COnlyCrossResult.alpha_pep_data.mod_site_list):
        if mod_index >= myMOD.fix_mod_num:
            mod_name = myMOD.var_mod_list[mod_index - myMOD.fix_mod_num]
            str_alpha_mod += str(site) + ':' + str(mod_name) + ';'
        elif mod_index > -1:
            mod_name = myMOD.fix_mod_list[mod_index]
            str_alpha_mod += str(site) + ':' + str(mod_name) + ';'
        else:
            pass
    for site, mod_index in enumerate(input_COnlyCrossResult.beta_pep_data.mod_site_list):
        if mod_index >= myMOD.fix_mod_num:
            mod_name = myMOD.var_mod_list[mod_index - myMOD.fix_mod_num]
            str_beta_mod += str(site) + ':' + str(mod_name) + ';'
        elif mod_index > -1:
            mod_name = myMOD.fix_mod_list[mod_index]
            str_beta_mod += str(site) + ':' + str(mod_name) + ';'
        else:
            pass

    protein_1_site = input_COnlyCrossResult.alpha_pep_data.start_pos + input_COnlyCrossResult.alpha_pep_cross_site
    protein_2_site = input_COnlyCrossResult.beta_pep_data.start_pos + input_COnlyCrossResult.beta_pep_cross_site

    str_protein = protein_list.ac[input_COnlyCrossResult.alpha_pep_data.pro_index] + '(' + str(protein_1_site) + ')-' + \
                  protein_list.ac[input_COnlyCrossResult.beta_pep_data.pro_index] + '(' + str(protein_2_site) + ')\\'
    for alpha_index, alpha_protein_index in enumerate(input_COnlyCrossResult.alpha_pep_data.protein_list):
        for beta_index, beta_protein_index in enumerate(input_COnlyCrossResult.beta_pep_data.protein_list):
            if (alpha_protein_index == input_COnlyCrossResult.alpha_pep_data.pro_index
                    and input_COnlyCrossResult.alpha_pep_data.protein_start_pos_list[
                        alpha_index] == input_COnlyCrossResult.alpha_pep_data.start_pos
                    and beta_protein_index == input_COnlyCrossResult.beta_pep_data.pro_index
                    and input_COnlyCrossResult.beta_pep_data.protein_start_pos_list[
                        beta_index] == input_COnlyCrossResult.beta_pep_data.start_pos):
                pass
            else:
                protein_1_site = input_COnlyCrossResult.alpha_pep_data.protein_start_pos_list[
                                     alpha_index] + input_COnlyCrossResult.alpha_pep_cross_site
                protein_2_site = input_COnlyCrossResult.beta_pep_data.protein_start_pos_list[
                                     beta_index] + input_COnlyCrossResult.beta_pep_cross_site
                str_protein += protein_list.ac[alpha_protein_index] + '(' + str(protein_1_site) + ')-' + \
                               protein_list.ac[beta_protein_index] + '(' + str(protein_2_site) + ')\\'
    protein_type = str(input_COnlyCrossResult.alpha_pep_data.is_target) + ',' + str(
        input_COnlyCrossResult.beta_pep_data.is_target)

    if input_COnlyCrossResult.alpha_pep_data.is_target + input_COnlyCrossResult.beta_pep_data.is_target == 2:
        protein_type_math = 1
    elif input_COnlyCrossResult.alpha_pep_data.is_target + input_COnlyCrossResult.beta_pep_data.is_target == 0:
        protein_type_math = -1
    else:
        protein_type_math = 0

    input_CWriteOnlyCrossResult.spectrum_title = input_COnlyCrossResult.spectrum_title
    input_CWriteOnlyCrossResult.spectrum_charge = str(input_COnlyCrossResult.spectrum_charge)
    input_CWriteOnlyCrossResult.spectrum_mass = '%.4f' % input_COnlyCrossResult.spectrum_mass

    alpha_sq = protein_list.sq_original[input_COnlyCrossResult.alpha_pep_data.pro_index][
               input_COnlyCrossResult.alpha_pep_data.start_pos: input_COnlyCrossResult.alpha_pep_data.end_pos]
    beta_sq = protein_list.sq_original[input_COnlyCrossResult.beta_pep_data.pro_index][
              input_COnlyCrossResult.beta_pep_data.start_pos: input_COnlyCrossResult.beta_pep_data.end_pos]
    input_CWriteOnlyCrossResult.peptide_sq = alpha_sq + '(' + str(
        input_COnlyCrossResult.alpha_pep_cross_site) + ')-' + beta_sq + '(' + str(
        input_COnlyCrossResult.beta_pep_cross_site) + ')'
    input_CWriteOnlyCrossResult.peptide_type = 'cross'
    input_CWriteOnlyCrossResult.linker = linker_name
    input_CWriteOnlyCrossResult.peptide_mass = '%.4f' % input_COnlyCrossResult.pep_mass
    input_CWriteOnlyCrossResult.peptide_modification = str_alpha_mod + '-' + str_beta_mod
    input_CWriteOnlyCrossResult.peptide_score = input_COnlyCrossResult.pep_score
    input_CWriteOnlyCrossResult.precursor_mass_error_Da = '%.4f' % input_COnlyCrossResult.pep_precursor_bias_Da
    input_CWriteOnlyCrossResult.precursor_mass_error_ppm = '%.4f' % input_COnlyCrossResult.pep_precursor_bias_ppm
    input_CWriteOnlyCrossResult.peptide_protein = str_protein
    input_CWriteOnlyCrossResult.protein_type_math = protein_type_math
    input_CWriteOnlyCrossResult.protein_type = protein_type

    # 判断是inter还是intra，默认intra的，因为第一个蛋白就能确定是不是intra的了
    # end_pos不包括这个位置
    if (input_COnlyCrossResult.alpha_pep_data.pro_index == input_COnlyCrossResult.beta_pep_data.pro_index):
        # 应该说来自同一个蛋白质的两个肽段有重叠部分就算inter的
        if input_COnlyCrossResult.alpha_pep_data.start_pos >= input_COnlyCrossResult.beta_pep_data.start_pos and input_COnlyCrossResult.alpha_pep_data.start_pos < input_COnlyCrossResult.beta_pep_data.end_pos:
            # 两条肽段来源位置有交集
            input_CWriteOnlyCrossResult.inter_protein = 1
        elif input_COnlyCrossResult.beta_pep_data.start_pos >= input_COnlyCrossResult.alpha_pep_data.start_pos and input_COnlyCrossResult.beta_pep_data.start_pos < input_COnlyCrossResult.alpha_pep_data.end_pos:
            # 两条肽段来源位置有交集
            input_CWriteOnlyCrossResult.inter_protein = 1
    elif input_COnlyCrossResult.alpha_pep_data.pro_index != input_COnlyCrossResult.beta_pep_data.pro_index:
        if input_COnlyCrossResult.alpha_pep_data.pro_index % 2 == 0:
            # 同一个蛋白的正库和反库就是intra的
            if input_COnlyCrossResult.alpha_pep_data.pro_index + 1 == input_COnlyCrossResult.beta_pep_data.pro_index:
                pass
            else:
                input_CWriteOnlyCrossResult.inter_protein = 1
        elif input_COnlyCrossResult.beta_pep_data.pro_index % 2 == 0:
            # 同一个蛋白的正库和反库就是intra的
            if input_COnlyCrossResult.alpha_pep_data.pro_index == input_COnlyCrossResult.beta_pep_data.pro_index + 1:
                pass
            else:
                input_CWriteOnlyCrossResult.inter_protein = 1
        else:
            input_CWriteOnlyCrossResult.inter_protein = 1

    input_CWriteOnlyCrossResult.feature = input_COnlyCrossResult.feature
    input_CWriteOnlyCrossResult.FDR = FDR
    input_CWriteOnlyCrossResult.rerank_score = rerank_score


# ===============================

def op_fill_CSinglePeptideResult(input_CSinglePeptideResult, mass, pep_score, pep_precursor_bias_Da,
                                 pep_precursor_bias_ppm, other_data=[]):
    input_CSinglePeptideResult.pep_mass = mass
    input_CSinglePeptideResult.pep_precursor_bias_Da = pep_precursor_bias_Da
    input_CSinglePeptideResult.pep_precursor_bias_ppm = pep_precursor_bias_ppm
    input_CSinglePeptideResult.pep_score = pep_score

    input_CSinglePeptideResult.other_data = other_data


def op_fill_CWriteSinglePeptdidResult(input_CWriteResultSinglePeptide, input_CSinglePeptideResult, myMOD, myLINK,
                                      protein_list, delta_score, FDR=0, rerank_score=0.0):
    str_mod = ''
    str_mono = ''

    if input_CSinglePeptideResult.loop_link:
        str_loop = myLINK.linker_1_data.linker_name
        if input_CSinglePeptideResult.loop_link:
            str_loop += "(" + str(input_CSinglePeptideResult.loop_site_1) + ");(" + str(
                input_CSinglePeptideResult.loop_site_2) + ")"
    else:
        str_loop = ''
    for site, mod_index in enumerate(input_CSinglePeptideResult.pep_data.mod_site_list):
        if mod_index >= myMOD.fix_mod_num + myMOD.var_mod_num:
            if input_CSinglePeptideResult.loop_link:
                pass
            else:
                input_CWriteResultSinglePeptide.momo_link = True
                mono_index = mod_index - myMOD.fix_mod_num - myMOD.var_mod_num
                if mono_index == 0:
                    str_mono = myLINK.linker_1_data.linker_name + "(" + str(site) + ");"
                elif mono_index == 0:
                    str_mono = myLINK.linker_2_data.linker_name + "(" + str(site) + ");"
        elif mod_index >= myMOD.fix_mod_num:
            mod_name = myMOD.var_mod_list[mod_index - myMOD.fix_mod_num]
            str_mod += str(site) + ':' + str(mod_name) + ';'
        elif mod_index > -1:
            mod_name = myMOD.fix_mod_list[mod_index]
            str_mod += str(site) + ':' + str(mod_name) + ';'
        else:
            pass

    str_protein = protein_list.ac[input_CSinglePeptideResult.pep_data.pro_index].split('\t')[0] + '\\'
    for protein_index in input_CSinglePeptideResult.pep_data.protein_list:
        str_protein += protein_list.ac[protein_index].split('\t')[0] + '\\'
    protein_type = str(input_CSinglePeptideResult.pep_data.is_target)

    if input_CSinglePeptideResult.pep_data.is_target == 1:
        protein_type_math = 1
    else:
        protein_type_math = -1
    input_CSinglePeptideResult.feature.delta_score = delta_score
    pep_sq = protein_list.sq_original[input_CSinglePeptideResult.pep_data.pro_index][
             input_CSinglePeptideResult.pep_data.start_pos: input_CSinglePeptideResult.pep_data.end_pos]
    input_CWriteResultSinglePeptide.peptide_sq = pep_sq
    if input_CWriteResultSinglePeptide.loop_link:
        input_CWriteResultSinglePeptide.peptide_type = 'loop-link'
    elif input_CWriteResultSinglePeptide.momo_link:
        input_CWriteResultSinglePeptide.peptide_type = 'mono-link'
    else:
        input_CWriteResultSinglePeptide.peptide_type = 'single'
    input_CWriteResultSinglePeptide.peptide_mass = '%.4f' % input_CSinglePeptideResult.pep_mass
    input_CWriteResultSinglePeptide.peptide_modification = str_mod
    input_CWriteResultSinglePeptide.peptide_score = input_CSinglePeptideResult.pep_score
    input_CWriteResultSinglePeptide.precursor_mass_error_Da = '%.4f' % input_CSinglePeptideResult.pep_precursor_bias_Da
    input_CWriteResultSinglePeptide.precursor_mass_error_ppm = '%.4f' % input_CSinglePeptideResult.pep_precursor_bias_ppm
    input_CWriteResultSinglePeptide.peptide_protein = str_protein
    input_CWriteResultSinglePeptide.protein_type_math = protein_type_math
    input_CWriteResultSinglePeptide.protein_type = protein_type

    input_CWriteResultSinglePeptide.peptide_mono = str_mono
    input_CWriteResultSinglePeptide.peptide_loop = str_loop

    input_CWriteResultSinglePeptide.feature = input_CSinglePeptideResult.feature
    input_CWriteResultSinglePeptide.FDR = FDR
    input_CWriteResultSinglePeptide.rerank_score = rerank_score


# ===============================

def op_fill_CCrossCrossResult(input_CCrossCrossResult, mass, pep_score, alpha_pep_score, beta_pep_score,
                              pep_precursor_bias_Da, pep_precursor_bias_ppm, alpha_match_data=[], beta_match_data=[],
                              other_data=[]):
    input_CCrossCrossResult.pep_mass = mass
    input_CCrossCrossResult.pep_precursor_bias_Da = pep_precursor_bias_Da
    input_CCrossCrossResult.pep_precursor_bias_ppm = pep_precursor_bias_ppm
    input_CCrossCrossResult.pep_score = pep_score
    input_CCrossCrossResult.alpha_pep_score = alpha_pep_score
    input_CCrossCrossResult.beta_pep_score = beta_pep_score

    input_CCrossCrossResult.other_data = other_data
    input_CCrossCrossResult.alpha_pep_match_data = alpha_match_data
    input_CCrossCrossResult.beta_pep_match_data = beta_match_data


def op_fill_CWriteCrossCrossResult(input_CWriteCrossCrossResult, input_CCrossCrossResult, myMOD, protein_list,
                                   linker_name, FDR=0, rerank_score=0.0):
    str_alpha_mod = ''
    str_beta_mod = ''
    for site, mod_index in enumerate(input_CCrossCrossResult.alpha_pep_data.mod_site_list):
        if mod_index >= myMOD.fix_mod_num:
            mod_name = myMOD.var_mod_list[mod_index - myMOD.fix_mod_num]
            str_alpha_mod += str(site) + ':' + str(mod_name) + ';'
        elif mod_index > -1:
            mod_name = myMOD.fix_mod_list[mod_index]
            str_alpha_mod += str(site) + ':' + str(mod_name) + ';'
        else:
            pass
    for site, mod_index in enumerate(input_CCrossCrossResult.beta_pep_data.mod_site_list):
        if mod_index >= myMOD.fix_mod_num:
            mod_name = myMOD.var_mod_list[mod_index - myMOD.fix_mod_num]
            str_beta_mod += str(site) + ':' + str(mod_name) + ';'
        elif mod_index > -1:
            mod_name = myMOD.fix_mod_list[mod_index]
            str_beta_mod += str(site) + ':' + str(mod_name) + ';'
        else:
            pass

    str_protein = ''
    for alpha_index, alpha_protein_index in enumerate(input_CCrossCrossResult.alpha_pep_data.protein_list):
        for beta_index, beta_protein_index in enumerate(input_CCrossCrossResult.beta_pep_data.protein_list):
            protein_alpha_1_site = input_CCrossCrossResult.alpha_pep_data.protein_start_pos_list[alpha_index] + \
                                   input_CCrossCrossResult.alpha_pep_site_data[0][0]
            protein_alpha_2_site = input_CCrossCrossResult.alpha_pep_data.protein_start_pos_list[alpha_index] + \
                                   input_CCrossCrossResult.alpha_pep_site_data[1][0]

            protein_beta_1_site = input_CCrossCrossResult.beta_pep_data.protein_start_pos_list[beta_index] + \
                                  input_CCrossCrossResult.beta_pep_site_data[0][0]
            protein_beta_2_site = input_CCrossCrossResult.beta_pep_data.protein_start_pos_list[beta_index] + \
                                  input_CCrossCrossResult.beta_pep_site_data[1][0]

            str_protein += protein_list.ac[alpha_protein_index] + '(' + str(protein_alpha_1_site) + ',' + str(
                protein_alpha_2_site) + ')-' + protein_list.ac[beta_protein_index] + '(' + str(
                protein_beta_1_site) + ',' + str(protein_beta_2_site) + ')\\'

    protein_type = str(input_CCrossCrossResult.alpha_pep_data.is_target) + ',' + str(
        input_CCrossCrossResult.beta_pep_data.is_target)

    if input_CCrossCrossResult.alpha_pep_data.is_target + input_CCrossCrossResult.beta_pep_data.is_target == 2:
        protein_type_math = 1
    elif input_CCrossCrossResult.alpha_pep_data.is_target + input_CCrossCrossResult.beta_pep_data.is_target == 0:
        protein_type_math = -1
    else:
        protein_type_math = 0

    input_CWriteCrossCrossResult.spectrum_title = input_CCrossCrossResult.spectrum_title
    input_CWriteCrossCrossResult.spectrum_charge = str(input_CCrossCrossResult.spectrum_charge)
    input_CWriteCrossCrossResult.spectrum_mass = '%.4f' % input_CCrossCrossResult.spectrum_mass

    alpha_sq = protein_list.sq_original[input_CCrossCrossResult.alpha_pep_data.protein_list[0]][
               input_CCrossCrossResult.alpha_pep_data.start_pos: input_CCrossCrossResult.alpha_pep_data.end_pos]
    beta_sq = protein_list.sq_original[input_CCrossCrossResult.beta_pep_data.protein_list[0]][
              input_CCrossCrossResult.beta_pep_data.start_pos: input_CCrossCrossResult.beta_pep_data.end_pos]

    if input_CCrossCrossResult.alpha_pep_site_data[0][2] == 1:
        alpha_pep_site = '1-' + linker_name[0] + ':' + str(input_CCrossCrossResult.alpha_pep_site_data[0][0]) + ';2-' + \
                         linker_name[1] + ':' + str(input_CCrossCrossResult.alpha_pep_site_data[1][0]) + ';'
    else:
        alpha_pep_site = '2-' + linker_name[1] + ':' + str(input_CCrossCrossResult.alpha_pep_site_data[0][0]) + ';1-' + \
                         linker_name[0] + ':' + str(input_CCrossCrossResult.alpha_pep_site_data[1][0]) + ';'

    if input_CCrossCrossResult.beta_pep_site_data[0][2] == 1:
        beta_pep_site = '1-' + linker_name[0] + ':' + str(input_CCrossCrossResult.beta_pep_site_data[0][0]) + ';2-' + \
                        linker_name[1] + ':' + str(input_CCrossCrossResult.beta_pep_site_data[1][0]) + ';'
    else:
        beta_pep_site = '2-' + linker_name[1] + ':' + str(input_CCrossCrossResult.beta_pep_site_data[0][0]) + ';1-' + \
                        linker_name[0] + ':' + str(input_CCrossCrossResult.beta_pep_site_data[1][0]) + ';'

    input_CWriteCrossCrossResult.peptide_sq = alpha_sq + '(' + alpha_pep_site + ')-' + beta_sq + '(' + beta_pep_site + ')'

    input_CWriteCrossCrossResult.peptide_type = 'cross-cross'
    if input_CCrossCrossResult.alpha_pep_site_data[0][2] == 1:
        input_CWriteCrossCrossResult.linker = '1-' + linker_name[0] + ';2-' + linker_name[1]
    else:
        input_CWriteCrossCrossResult.linker = '2-' + linker_name[1] + ';1-' + linker_name[0]
    input_CWriteCrossCrossResult.peptide_mass = '%.4f' % input_CCrossCrossResult.pep_mass
    input_CWriteCrossCrossResult.peptide_modification = str_alpha_mod + '-' + str_beta_mod
    input_CWriteCrossCrossResult.peptide_score = input_CCrossCrossResult.pep_score
    input_CWriteCrossCrossResult.precursor_mass_error_Da = '%.4f' % input_CCrossCrossResult.pep_precursor_bias_Da
    input_CWriteCrossCrossResult.precursor_mass_error_ppm = '%.4f' % input_CCrossCrossResult.pep_precursor_bias_ppm
    input_CWriteCrossCrossResult.peptide_protein = str_protein
    input_CWriteCrossCrossResult.protein_type_math = protein_type_math
    input_CWriteCrossCrossResult.protein_type = protein_type

    # 判断是inter还是intra，默认intra的，因为第一个蛋白就能确定是不是intra的了
    if (input_CCrossCrossResult.alpha_pep_data.protein_list[0] == input_CCrossCrossResult.beta_pep_data.protein_list[
        0]):
        # 应该说来自同一个蛋白质的两个肽段有重叠部分就算inter的
        if input_CCrossCrossResult.alpha_pep_data.start_pos >= input_CCrossCrossResult.beta_pep_data.start_pos and input_CCrossCrossResult.alpha_pep_data.start_pos < input_CCrossCrossResult.beta_pep_data.end_pos:
            # 两条肽段来源位置有交集
            input_CWriteCrossCrossResult.inter_protein = 1
        elif input_CCrossCrossResult.beta_pep_data.start_pos >= input_CCrossCrossResult.alpha_pep_data.start_pos and input_CCrossCrossResult.beta_pep_data.start_pos < input_CCrossCrossResult.alpha_pep_data.end_pos:
            # 两条肽段来源位置有交集
            input_CWriteCrossCrossResult.inter_protein = 1
    elif input_CCrossCrossResult.alpha_pep_data.protein_list[0] != input_CCrossCrossResult.beta_pep_data.protein_list[
        0]:
        if input_CCrossCrossResult.alpha_pep_data.protein_list[0] % 2 == 0:
            # 同一个蛋白的正库和反库就是intra的
            if input_CCrossCrossResult.alpha_pep_data.protein_list[0] + 1 == \
                    input_CCrossCrossResult.beta_pep_data.protein_list[0]:
                pass
            else:
                input_CWriteCrossCrossResult.inter_protein = 1
        elif input_CCrossCrossResult.beta_pep_data.protein_list[0] % 2 == 0:
            # 同一个蛋白的正库和反库就是intra的
            if input_CCrossCrossResult.alpha_pep_data.protein_list[0] == \
                    input_CCrossCrossResult.beta_pep_data.protein_list[0] + 1:
                pass
            else:
                input_CWriteCrossCrossResult.inter_protein = 1
        else:
            input_CWriteCrossCrossResult.inter_protein = 1

    input_CWriteCrossCrossResult.feature = input_CCrossCrossResult.feature
    input_CWriteCrossCrossResult.FDR = FDR
    input_CWriteCrossCrossResult.rerank_score = rerank_score


# ===============================

def op_fill_CCrossLoopResult(input_CCrossLoopResult, mass, pep_score, alpha_pep_score, beta_pep_score,
                             pep_precursor_bias_Da, pep_precursor_bias_ppm, alpha_match_data=[], beta_match_data=[],
                             other_data=[]):
    input_CCrossLoopResult.pep_mass = mass
    input_CCrossLoopResult.pep_precursor_bias_Da = pep_precursor_bias_Da
    input_CCrossLoopResult.pep_precursor_bias_ppm = pep_precursor_bias_ppm
    input_CCrossLoopResult.pep_score = pep_score
    input_CCrossLoopResult.alpha_pep_score = alpha_pep_score
    input_CCrossLoopResult.beta_pep_score = beta_pep_score

    input_CCrossLoopResult.other_data = other_data
    input_CCrossLoopResult.alpha_pep_match_data = alpha_match_data
    input_CCrossLoopResult.beta_pep_match_data = beta_match_data


def op_fill_CWriteCrossLoopResult(input_CWriteCrossLoopResult, input_CCrossLoopResult, myMOD, protein_list, linker_name,
                                  FDR=0, rerank_score=0.0):
    str_alpha_mod = ''
    str_beta_mod = ''
    for site, mod_index in enumerate(input_CCrossLoopResult.alpha_pep_data.mod_site_list):
        if mod_index >= myMOD.fix_mod_num:
            mod_name = myMOD.var_mod_list[mod_index - myMOD.fix_mod_num]
            str_alpha_mod += str(site) + ':' + str(mod_name) + ';'
        elif mod_index > -1:
            mod_name = myMOD.fix_mod_list[mod_index]
            str_alpha_mod += str(site) + ':' + str(mod_name) + ';'
        else:
            pass
    for site, mod_index in enumerate(input_CCrossLoopResult.beta_pep_data.mod_site_list):
        if mod_index >= myMOD.fix_mod_num:
            mod_name = myMOD.var_mod_list[mod_index - myMOD.fix_mod_num]
            str_beta_mod += str(site) + ':' + str(mod_name) + ';'
        elif mod_index > -1:
            mod_name = myMOD.fix_mod_list[mod_index]
            str_beta_mod += str(site) + ':' + str(mod_name) + ';'
        else:
            pass

    str_protein = ''
    for alpha_index, alpha_protein_index in enumerate(input_CCrossLoopResult.alpha_pep_data.protein_list):
        for beta_index, beta_protein_index in enumerate(input_CCrossLoopResult.beta_pep_data.protein_list):
            # alpha_protein_site =
            if input_CCrossLoopResult.alpha_pep_site_data[0][2] == 1:
                alpha_protein_site = input_CCrossLoopResult.alpha_pep_site_data[0][0] + \
                                     input_CWriteCrossLoopResult.input_CCrossLoopResult.alpha_pep_data.protein_start_pos_list[
                                         alpha_index]

                beta_protein_site_1 = input_CCrossLoopResult.beta_pep_site_data[0][0] + \
                                      input_CWriteCrossLoopResult.input_CCrossLoopResult.beta_pep_data.protein_start_pos_list[
                                          beta_index]
                beta_protein_site_2 = input_CCrossLoopResult.beta_pep_site_data[1][0] + \
                                      input_CWriteCrossLoopResult.input_CCrossLoopResult.beta_pep_data.protein_start_pos_list[
                                          beta_index]
                beta_protein_site_3 = input_CCrossLoopResult.beta_pep_site_data[2][0] + \
                                      input_CWriteCrossLoopResult.input_CCrossLoopResult.beta_pep_data.protein_start_pos_list[
                                          beta_index]
                str_protein += protein_list.ac[alpha_protein_index] + '(' + str(alpha_protein_site) + ')-' + \
                               protein_list.ac[beta_protein_index] + '(' + str(beta_protein_site_1) + ',' + str(
                    beta_protein_site_2) + ',' + str(beta_protein_site_3) + ')\\'
            else:

                beta_protein_site = input_CCrossLoopResult.beta_pep_site_data[0][0] + \
                                    input_CWriteCrossLoopResult.input_CCrossLoopResult.beta_pep_data.protein_start_pos_list[
                                        beta_index]

                alpha_protein_site_1 = input_CCrossLoopResult.alpha_pep_site_data[0][0] + \
                                       input_CWriteCrossLoopResult.input_CCrossLoopResult.beta_pep_data.protein_start_pos_list[
                                           alpha_index]
                alpha_protein_site_2 = input_CCrossLoopResult.alpha_pep_site_data[1][0] + \
                                       input_CWriteCrossLoopResult.input_CCrossLoopResult.beta_pep_data.protein_start_pos_list[
                                           alpha_index]
                alpha_protein_site_3 = input_CCrossLoopResult.alpha_pep_site_data[2][0] + \
                                       input_CWriteCrossLoopResult.input_CCrossLoopResult.beta_pep_data.protein_start_pos_list[
                                           alpha_index]
                str_protein += protein_list.ac[alpha_protein_index] + '(' + str(alpha_protein_site_1) + ',' + str(
                    alpha_protein_site_2) + ',' + str(alpha_protein_site_3) + ')-' + protein_list.ac[
                                   beta_protein_index] + '(' + str(beta_protein_site) + ')\\'

    protein_type = str(input_CCrossLoopResult.alpha_pep_data.is_target) + ',' + str(
        input_CCrossLoopResult.beta_pep_data.is_target)

    if input_CCrossLoopResult.alpha_pep_data.is_target + input_CCrossLoopResult.beta_pep_data.is_target == 2:
        protein_type_math = 1
    elif input_CCrossLoopResult.alpha_pep_data.is_target + input_CCrossLoopResult.beta_pep_data.is_target == 0:
        protein_type_math = -1
    else:
        protein_type_math = 0

    input_CWriteCrossLoopResult.spectrum_title = input_CCrossLoopResult.spectrum_title
    input_CWriteCrossLoopResult.spectrum_charge = str(input_CCrossLoopResult.spectrum_charge)
    input_CWriteCrossLoopResult.spectrum_mass = '%.4f' % input_CCrossLoopResult.spectrum_mass

    alpha_sq = protein_list.sq_original[input_CCrossLoopResult.alpha_pep_data.protein_list[0]][
               input_CCrossLoopResult.alpha_pep_data.start_pos: input_CCrossLoopResult.alpha_pep_data.end_pos]
    beta_sq = protein_list.sq_original[input_CCrossLoopResult.beta_pep_data.protein_list[0]][
              input_CCrossLoopResult.beta_pep_data.start_pos: input_CCrossLoopResult.beta_pep_data.end_pos]

    if input_CCrossLoopResult.alpha_pep_site_data[0][2] == 1:
        alpha_pep_site = linker_name[0] + ':' + str(input_CCrossLoopResult.alpha_pep_site_data[0][0]) + ';'
        if len(input_CCrossLoopResult.alpha_pep_site_data) == 3:
            alpha_pep_site += linker_name[1] + ':' + str(input_CCrossLoopResult.alpha_pep_site_data[1][0]) + ';' + \
                              linker_name[1] + ':' + str(input_CCrossLoopResult.alpha_pep_site_data[2][0]) + ';'
    else:
        alpha_pep_site = linker_name[1] + ':' + str(input_CCrossLoopResult.alpha_pep_site_data[1][0]) + ';'
        if len(input_CCrossLoopResult.alpha_pep_site_data) == 3:
            alpha_pep_site += linker_name[0] + ':' + str(input_CCrossLoopResult.alpha_pep_site_data[1][0]) + ';' + \
                              linker_name[1] + ':' + str(input_CCrossLoopResult.alpha_pep_site_data[2][0]) + ';'

    if input_CCrossLoopResult.beta_pep_site_data[0][2] == 1:
        beta_pep_site = linker_name[0] + ':' + str(input_CCrossLoopResult.beta_pep_site_data[0][0]) + ';'
        if len(input_CCrossLoopResult.beta_pep_site_data) == 3:
            beta_pep_site += linker_name[1] + ':' + str(input_CCrossLoopResult.beta_pep_site_data[1][0]) + ';' + \
                             linker_name[1] + ':' + str(input_CCrossLoopResult.beta_pep_site_data[2][0]) + ';'
    else:
        beta_pep_site = linker_name[1] + ':' + str(input_CCrossLoopResult.beta_pep_site_data[1][0]) + ';'
        if len(input_CCrossLoopResult.beta_pep_site_data) == 3:
            beta_pep_site += linker_name[0] + ':' + str(input_CCrossLoopResult.beta_pep_site_data[1][0]) + ';' + \
                             linker_name[0] + ':' + str(input_CCrossLoopResult.beta_pep_site_data[2][0]) + ';'

    input_CWriteCrossLoopResult.peptide_sq = alpha_sq + '(' + alpha_pep_site + ')-' + beta_sq + '(' + beta_pep_site + ')'

    input_CWriteCrossLoopResult.peptide_type = 'cross-loop'
    input_CWriteCrossLoopResult.linker = linker_name[0] + '-' + linker_name[1]
    input_CWriteCrossLoopResult.peptide_mass = '%.4f' % input_CCrossLoopResult.pep_mass
    input_CWriteCrossLoopResult.peptide_modification = str_alpha_mod + '-' + str_beta_mod
    input_CWriteCrossLoopResult.peptide_score = input_CCrossLoopResult.pep_score
    input_CWriteCrossLoopResult.precursor_mass_error_Da = '%.4f' % input_CCrossLoopResult.pep_precursor_bias_Da
    input_CWriteCrossLoopResult.precursor_mass_error_ppm = '%.4f' % input_CCrossLoopResult.pep_precursor_bias_ppm
    input_CWriteCrossLoopResult.peptide_protein = str_protein
    input_CWriteCrossLoopResult.protein_type_math = protein_type_math
    input_CWriteCrossLoopResult.protein_type = protein_type

    # 判断是inter还是intra，默认intra的，因为第一个蛋白就能确定是不是intra的了
    if (input_CCrossLoopResult.alpha_pep_data.protein_list[0] == input_CCrossLoopResult.beta_pep_data.protein_list[0]):
        # 应该说来自同一个蛋白质的两个肽段有重叠部分就算inter的
        if input_CCrossLoopResult.alpha_pep_data.start_pos >= input_CCrossLoopResult.beta_pep_data.start_pos and input_CCrossLoopResult.alpha_pep_data.start_pos < input_CCrossLoopResult.beta_pep_data.end_pos:
            # 两条肽段来源位置有交集
            input_CWriteCrossLoopResult.inter_protein = 1
        elif input_CCrossLoopResult.beta_pep_data.start_pos >= input_CCrossLoopResult.alpha_pep_data.start_pos and input_CCrossLoopResult.beta_pep_data.start_pos < input_CCrossLoopResult.alpha_pep_data.end_pos:
            # 两条肽段来源位置有交集
            input_CWriteCrossLoopResult.inter_protein = 1
    elif input_CCrossLoopResult.alpha_pep_data.protein_list[0] != input_CCrossLoopResult.beta_pep_data.protein_list[0]:
        if input_CCrossLoopResult.alpha_pep_data.protein_list[0] % 2 == 0:
            # 同一个蛋白的正库和反库就是intra的
            if input_CCrossLoopResult.alpha_pep_data.protein_list[0] + 1 == \
                    input_CCrossLoopResult.beta_pep_data.protein_list[0]:
                pass
            else:
                input_CWriteCrossLoopResult.inter_protein = 1
        elif input_CCrossLoopResult.beta_pep_data.protein_list[0] % 2 == 0:
            # 同一个蛋白的正库和反库就是intra的
            if input_CCrossLoopResult.alpha_pep_data.protein_list[0] == \
                    input_CCrossLoopResult.beta_pep_data.protein_list[0] + 1:
                pass
            else:
                input_CWriteCrossLoopResult.inter_protein = 1
        else:
            input_CWriteCrossLoopResult.inter_protein = 1

    input_CWriteCrossLoopResult.feature = input_CCrossLoopResult.feature
    input_CWriteCrossLoopResult.FDR = FDR
    input_CWriteCrossLoopResult.rerank_score = rerank_score


def opGetStartAndEndForProfile(input_profile, input_seed, input_cutoff, input_n_hole):  # 就这个名字格式特殊

    nHoleLeft = 0
    nHoleRight = 0

    i_middle = input_seed

    if i_middle > 0:
        i_left = i_middle - 1
    else:
        i_left = i_middle

    i_right = i_middle

    result = [i_left, i_right]  # start and end

    int_left = input_profile[i_left]
    int_right = input_profile[i_right]

    int_max = int_left
    int_thr = int_max * input_cutoff

    walkLeft = True
    walkRight = True

    while True:

        if walkLeft or walkRight:
            pass
        else:
            break

        if i_left > 0:

            if walkLeft:
                i_left = i_left - 1

        else:

            walkLeft = False

        if i_right < len(input_profile) - 1:

            if walkRight:
                i_right = i_right + 1

        else:

            walkRight = False

        int_left = input_profile[i_left]
        int_right = input_profile[i_right]

        # max
        if int_max < int_left:
            int_max = int_left

        if int_max < int_right:
            int_max = int_right

        int_thr = int_max * input_cutoff

        # hole
        if int_left < int_thr and walkLeft:
            nHoleLeft = nHoleLeft + 1

        if int_right < int_thr and walkRight:
            nHoleRight = nHoleRight + 1

        if nHoleLeft == input_n_hole:
            walkLeft = False

        if nHoleRight == input_n_hole:
            walkRight = False

    result[0] = i_left
    result[1] = i_right

    return result
