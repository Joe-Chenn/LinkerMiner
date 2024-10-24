# -*- mode: python ; coding: utf-8 -*-
from MSData import CMatchOnePeptideScore, CMatchIonScore

import math


def op_get_peptide_continue_score(input_continue_data, out_continue_data, out_pep_cover):

    TAG_LEN = 3
    for one_continue_data in input_continue_data:
        max_tag_num = len(one_continue_data) // TAG_LEN - 1
        tmp_list = []
        tmp_continue_len = 0
        tmp_continue_score = 0.0
        all_match_num = 0.0
        cover_tag_num = 0.0
        for one_continue_ion_score in one_continue_data:
            if one_continue_ion_score > 0.001:
                tmp_continue_len += 1
                all_match_num += 1
                tmp_continue_score += one_continue_ion_score
            else:
                if tmp_continue_len > 0:
                    cover_tag_num += 1
                if tmp_continue_len < TAG_LEN:
                    pass
                else:
                    tmp_list.append(tmp_continue_score)
                tmp_continue_score = 0.0
                tmp_continue_len = 0
        if tmp_continue_len > 0:
            cover_tag_num += 1
        if tmp_continue_len < TAG_LEN:
            pass
        else:
            tmp_list.append(tmp_continue_score)

        if len(tmp_list) == 0:
            out_continue_data.append(0.0)
        else:
            tmp_continue_score = 0.0
            tag_num = 0
            for tmp_index, one_continue_score in enumerate(tmp_list):
                tmp_continue_score += one_continue_score / TAG_LEN
                tag_num += 1
            if tag_num == 0:
                continue_score = tmp_continue_score
            else:
                continue_score = tmp_continue_score / tag_num
            out_continue_data.append(continue_score)
        if cover_tag_num == 0:
            out_pep_cover.append(0.0)
        else:
            out_pep_cover.append(all_match_num / len(one_continue_data) / cover_tag_num)

def op_get_peptide_by_continue_data(input_b_ion_list, input_y_ion_list, ion_len, spectrum, is_ppm_fra, tol_fra, multi):

    # 返回b、y和by离子连续性
    b_continue = [0 for i in range(ion_len)]
    y_continue = [0 for i in range(ion_len)]
    by_continue = [0 for i in range(ion_len)]

    if len(input_b_ion_list) > len(input_y_ion_list):
        input_list_len = len(input_y_ion_list)
    else:
        input_list_len = len(input_b_ion_list)

    b_continue_data = [0 for i in range(input_list_len)]
    y_continue_data = [0 for i in range(input_list_len)]

    # 这里是算误差允许范围，ppm的max_me就是20，max_me是0.02
    if is_ppm_fra:
        max_bias = tol_fra * 1e6
    else:
        max_bias = tol_fra
    # 遍历候选肽段的碎片离子
    index_num = 0
    for input_list_index in range(input_list_len):
        if b_continue[index_num] == 1:
            pass
        else:
            b_ion = input_b_ion_list[input_list_index]
            is_match = op_check_one_ion_is_match(b_ion, spectrum, is_ppm_fra, tol_fra, max_bias, multi)
            if is_match == 0:
                pass
            else:
                b_continue_data[index_num] = b_ion
                b_continue[index_num] = 1
                by_continue[index_num] = 1

        if y_continue[index_num] == 1:
            pass
        else:
            y_ion = input_y_ion_list[input_list_index]
            is_match = op_check_one_ion_is_match(y_ion, spectrum, is_ppm_fra, tol_fra, max_bias, multi)
            if is_match == 0:
                pass
            else:
                y_continue_data[index_num] = y_ion
                y_continue[index_num] = 1
                by_continue[ion_len - index_num - 1] = 1
        index_num += 1
        if index_num == ion_len:
            index_num = 0
    return [b_continue, y_continue, by_continue]

def op_check_one_ion_is_match(input_ion, spectrum, is_ppm_fra, tol_fra, max_bias, multi):

    flag = 0
    if is_ppm_fra:
        tol_mass = input_ion * tol_fra
    else:
        tol_mass = tol_fra
    start_mass = input_ion - tol_mass
    end_mass = input_ion + tol_mass
    # 对这个偏差范围超出限制的就跳过下一个
    if end_mass <= 0 or start_mass > spectrum.peaks[-1].mz:
        pass
    else:
        if start_mass < 0:
            start_mass = 0.0
        if end_mass > spectrum.peaks[-1].mz:
            end_mass = spectrum.peaks[-1].mz
        s_start_ind = spectrum.peaks_index[int(start_mass * multi)]
        # 确定肽段的这个碎片离子所对应谱图离子所在区间
        while s_start_ind < len(spectrum.peaks) and spectrum.peaks[s_start_ind].mz < start_mass:
            s_start_ind += 1
        if int(end_mass * multi) >= len(spectrum.peaks_index):
            s_end_ind = spectrum.peaks_index[-1]
        else:
            s_end_ind = spectrum.peaks_index[int(end_mass * multi)]
        while s_end_ind < len(spectrum.peaks) and spectrum.peaks[s_end_ind].mz <= end_mass:
            s_end_ind += 1
        if s_start_ind >= s_end_ind:
            pass
        else:
            for i in range(s_start_ind, s_end_ind):
                if is_ppm_fra:
                    one_mz_bias = (spectrum.peaks[i].mz - input_ion) * 1e6 / input_ion
                else:
                    one_mz_bias = spectrum.peaks[i].mz - input_ion
                abs_one_mz_bias = abs(one_mz_bias) / max_bias
                if abs_one_mz_bias < max_bias:
                    # 这就算匹配上了
                    flag = 1
                    break
    return flag

def op_get_match_ion_score(spectrum, ions, ion_num, is_ppm_fra, tol_fra, multi):
    '''
    # 这个是两个列表指针的搜索方法
    :param spectrum:CSpectrum
    :param ions: [ion]
    :param ion_num: 对应每个离子有多少个相同的
    :param is_ppm_fra: 偏差类型 1：ppm 0: Da
    :param tol_fra: 偏差范围
    :param multi: 可能没用了
    :return:
    '''
    # 这里要求输入的离子是排好序得
    match_num = 0
    match_inten_sum = 0

    match_moz = []
    match_ion_score = 0
    list_Da = []
    list_ppm = []
    # 这里是算误差允许范围，ppm的max_me就是20，max_me是0.02
    record_ions = []

    max_bias_Da = tol_fra
    max_bias_ppm = tol_fra * 1e6

    if is_ppm_fra:
        max_bias = max_bias_ppm
    else:
        max_bias = max_bias_Da
    # 遍历候选肽段的碎片离子
    ion_index = 0
    peak_index = 0
    peaks = spectrum.peaks
    # 双指针查找
    while True:
        if ion_index == len(ions) or peak_index == len(peaks):
            break
        else:
            cur_ion_mz = ions[ion_index]
            cur_peak_mz = peaks[peak_index].mz
            if is_ppm_fra:
                tol_mass = cur_ion_mz * tol_fra
            else:
                tol_mass = tol_fra
            # 算现在现在这个碎片离子的误差范围
            cur_ion_start_mz = cur_ion_mz - tol_mass
            cur_ion_end_mz = cur_ion_mz + tol_mass
            if cur_peak_mz > cur_ion_start_mz:
                # 当前谱图碎片离子大于查找碎片离子最小要求
                if cur_peak_mz > cur_ion_end_mz:
                    # 当前谱图碎片离子大于碎片离子最大要求
                    ion_index += 1
                else:
                    # 当前谱图碎片离子正好在碎片离子的最大最小要求内
                    tmp_peak_index = peak_index
                    # 记录当前查找范围的最大值
                    record_cur_match_ion_max_score = 0.0
                    record_cur_match_ion_max_score_index = -1
                    record_cur_match_ion_max_score_Da = max_bias_Da
                    record_cur_match_ion_max_score_ppm = max_bias_ppm

                    one_ion_score = 0.0
                    # 开始往这个谱图碎片离子右边走找最大的分数的
                    while True:
                        if tmp_peak_index > len(peaks) - 1:
                            break
                        if peaks[tmp_peak_index].mz > cur_ion_end_mz:
                            break
                        # 获取离子偏差
                        bias_Da = peaks[tmp_peak_index].mz - cur_ion_mz
                        bias_ppm = bias_Da * 1e6 / cur_ion_mz
                        abs_bias_ppm = abs(bias_ppm)
                        if is_ppm_fra:
                            if abs_bias_ppm > max_bias:
                                tmp_peak_index += 1
                                continue
                        else:
                            if bias_Da > max_bias:
                                tmp_peak_index += 1
                                continue
                        # compute Score
                        one_ionbias_score = math.sin(0.785398 * (1 + abs(1 - abs_bias_ppm / max_bias)))
                        one_ionint_score = math.sin(0.785398 * (1 + peaks[tmp_peak_index].inten / spectrum.max_int))
                        one_ion_score = one_ionint_score + one_ionbias_score
                        # 记录当前最高分的谱图碎片离子
                        if one_ion_score > record_cur_match_ion_max_score:
                            record_cur_match_ion_max_score = one_ion_score
                            record_cur_match_ion_max_score_index = tmp_peak_index
                            record_cur_match_ion_max_score_ppm = bias_ppm
                            record_cur_match_ion_max_score_Da = bias_Da
                        tmp_peak_index += 1
                    if record_cur_match_ion_max_score_index == -1:
                        peak_index += 1
                    else:
                        record_ions.append(ions[ion_index])
                        match_ion_score += record_cur_match_ion_max_score
                        match_inten_sum += spectrum.peaks[record_cur_match_ion_max_score_index].inten
                        list_Da.append(record_cur_match_ion_max_score_Da)
                        list_ppm.append(record_cur_match_ion_max_score_ppm)
                        match_num += ion_num[ion_index]
                        # 退出循环时是当前离子在谱图中的最高分离子，下一个离子肯定得从这个谱图中的最高分离子往右找
                        peak_index = record_cur_match_ion_max_score_index + 1
                    ion_index += 1

            else:
                # 当前谱图碎片离子小于查找碎片离子最小要求
                peak_index += 1

    if len(ions) == 0:
        match_ion_num_percent = 0
    else:
        match_ion_num_percent = len(list_ppm) / len(ions)

    if spectrum.max_int == 0:
        match_ion_inten_percent = 0
    else:
        match_ion_inten_percent = match_inten_sum / spectrum.max_int

    if len(spectrum.peaks) == 0:
        match_spe_ion_percent = 0
    else:
        match_spe_ion_percent = len(list_ppm) / len(spectrum.peaks)

    if spectrum.all_int == 0:
        match_spe_inten_percent = 0
    else:
        match_spe_inten_percent = match_inten_sum / spectrum.all_int

    spe_inten_sum = spectrum.all_int

    return CMatchIonScore(match_ion_score, list_Da, list_ppm,
                          match_ion_num_percent, match_ion_inten_percent,
                          match_spe_ion_percent, match_spe_inten_percent,
                          match_num, match_inten_sum, spe_inten_sum,
                          match_moz)

'''
def op_get_match_ion_score(specturm, ions, ion_num, is_ppm_fra, tol_fra, multi):
    # peak_index是建一个0-谱图碎片离子最大质量的区间，1mz为一个区间，反映谱图碎片离子分布，之后匹配直接用这个索引看这里面的碎片离子进行匹配
    # multi是放缩谱图质荷比区间大小，默认是1mz
    # 把离子去重后都输入到这里面，只管匹配得分
    match_num = 0
    match_inten_sum = 0.0
    match_peak_index = []
    match_moz = []
    match_ion_score = 0.0
    list_ppm = []

    # 这里是算误差允许范围，ppm的max_me就是20，max_me是0.02
    if is_ppm_fra:
        max_bias = tol_fra * 1e6
    else:
        max_bias = tol_fra
    # 遍历候选肽段的碎片离子

    for mi, mz in enumerate(ions):
        # 确定碎片允许的偏差范围

        if is_ppm_fra:
            tol_mass = mz * tol_fra
        else:
            tol_mass = tol_fra
        start_mass = mz - tol_mass
        end_mass = mz + tol_mass
        # 对这个偏差范围超出限制的就跳过下一个
        if end_mass <= 0:
            continue
        if start_mass > specturm.peaks[-1].mz:
            continue
        if start_mass < 0:
            start_mass = 0.0
        if end_mass > specturm.peaks[-1].mz:
            end_mass = specturm.peaks[-1].mz
        s_start_ind = specturm.peaks_index[int(start_mass * multi)]
        # 确定肽段的这个碎片离子所对应谱图离子所在区间
        while s_start_ind < len(specturm.peaks) and specturm.peaks[s_start_ind].mz < start_mass:
            s_start_ind += 1
        if int(end_mass * multi) >= len(specturm.peaks_index):
            s_end_ind = specturm.peaks_index[-1]
        else:
            s_end_ind = specturm.peaks_index[int(end_mass * multi)]
        while s_end_ind < len(specturm.peaks) and specturm.peaks[s_end_ind].mz <= end_mass:
            s_end_ind += 1
        if s_start_ind >= s_end_ind:
            continue
        cur_ion_max_score = 0.0
        cur_ion_max_score_index = -1
        cur_ion_max_score_bias = max_bias
        cur_ion_max_score_inten = 0.0
        # 找这里面分最高的峰来计分
        for i in range(s_start_ind, s_end_ind):
            if is_ppm_fra:
                one_mz_bias = (specturm.peaks[i].mz - mz) * 1e6 / mz
            else:
                one_mz_bias = specturm.peaks[i].mz - mz
            # 获取离子偏差
            abs_one_mz_bias = abs(one_mz_bias) / max_bias
            if abs_one_mz_bias > max_bias:
                continue
            else:
                # compute Score
                one_ionbias_score = math.sin(0.785398 * (1 + 1 - abs_one_mz_bias))
                one_ionint_score = math.sin(0.785398 * (1 + specturm.peaks[i].inten / specturm.max_int))
                one_ion_score = math.sin(0.785398 * (1 + one_ionbias_score * one_ionint_score))
                if cur_ion_max_score <= one_ion_score:
                    cur_ion_max_score = one_ion_score
                    cur_ion_max_score_index = i
                    cur_ion_max_score_bias = abs_one_mz_bias
                    cur_ion_max_score_inten = specturm.peaks[i].inten
                else:
                    pass
        if cur_ion_max_score_index == -1:
            continue
        else:
            # 最大分的这个
            match_num += ion_num[mi]
            list_ppm.append(cur_ion_max_score_bias)
            match_moz.append(mz)
            # 进到这里算匹配上
            # 记录峰强度占谱图总强度的比例
            # 这根峰被用过，离子分数仍然加上去
            if cur_ion_max_score_index not in match_peak_index:
                match_peak_index.append(cur_ion_max_score_index)
                match_ion_score += cur_ion_max_score
                match_inten_sum += cur_ion_max_score_inten

    match_spe_inten_percent = match_inten_sum / specturm.all_int
    match_spe_ion_percent = len(match_peak_index) / len(specturm.peaks)
    match_spe_inten_score = match_inten_sum / specturm.max_int
    match_spe_score = match_spe_inten_score * match_spe_ion_percent
    spe_inten_sum = specturm.all_int
    spe_inten_score = specturm.all_int / specturm.max_int
    return CMatchIonScore(match_ion_score, match_spe_score, match_spe_inten_score, match_spe_ion_percent, match_spe_inten_percent, list_ppm, match_num, match_inten_sum, spe_inten_sum, spe_inten_score, match_moz)
'''

def op_get_continue_score(peak_index, spe_peaks, spe_max_int, spe_all_int, input_ResPeptidedata, is_ppm_fra, tol_fra,  multi):
    # peak_index是建一个0-谱图碎片离子最大质量的区间，1mz为一个区间，反映谱图碎片离子分布，之后匹配直接用这个索引看这里面的碎片离子进行匹配
    # multi是放缩谱图质荷比区间大小，默认是1mz
    match_num = 0
    match_inten = 0.0
    match_peak_index, match_ion_name, match_peak2ion = [], [], []
    match_ion_score = 0.0
    list_ppm = []
    # 这里是算误差允许范围，ppm的max_me就是20，max_me是0.02
    if is_ppm_fra:
        max_bias = tol_fra * 1e6
    else:
        max_bias = tol_fra
    # 遍历候选肽段的碎片离子

    last_ion_type = ''
    last_ion_charge = 0
    ion_type_dic = {}
    data_all_ioncontinue_score = []
    one_type_ion_list = []

    for mi, mz in enumerate(input_ResPeptidedata.list_ion):
        # 确定碎片允许的偏差范围
        # 下面这里是判断离子类型和电荷数
        one_mz_flag = input_ResPeptidedata.list_ion_flag[mi]
        ion_charge = 0
        ion_type_index = 0
        for i, item in enumerate(one_mz_flag):
            if item >= '0' and item <= '9':
                ion_type_index = i
                break
        for i, item in enumerate(one_mz_flag):
            if item == '+':
                ion_charge += 1
        ion_type = one_mz_flag[:ion_type_index]

        #  判断离子类型
        if last_ion_type == ion_type:
            if ion_charge == last_ion_charge:
                # 上面两个都表示不是同一类型离子，进到这里面算是同一个类型电荷数的离子
                one_type_ion_list.append(0)  # 先构建这个类型碎片离子连续性列表
            else:
                # 新的类型离子
                data_all_ioncontinue_score.append([last_ion_type, last_ion_charge, one_type_ion_list])
                last_ion_charge = ion_charge
                last_ion_type = ion_type
                one_type_ion_list = [0]
                if ion_type in ion_type_dic.keys():
                    ion_type_dic[ion_type] = ion_type_dic[ion_type] + [len(data_all_ioncontinue_score)]
                else:
                    ion_type_dic[ion_type] = [len(data_all_ioncontinue_score)]
        else:
            # 进到这里也是一个新的类型离子，要对之前的汇总一下
            if last_ion_type == '':
                pass
            else:
                data_all_ioncontinue_score.append([last_ion_type, last_ion_charge, one_type_ion_list])
            if ion_type in ion_type_dic.keys():
                ion_type_dic[ion_type] = ion_type_dic[ion_type] + [len(data_all_ioncontinue_score)]
            else:
                ion_type_dic[ion_type] = [len(data_all_ioncontinue_score)]
            last_ion_charge = ion_charge
            last_ion_type = ion_type
            one_type_ion_list = [0]

        if is_ppm_fra:
            tol_mass = mz * tol_fra
        else:
            tol_mass = tol_fra
        start_mass = mz - tol_mass
        end_mass = mz + tol_mass
        # 对这个偏差范围超出限制的就跳过下一个
        if end_mass <= 0:
            continue
        if start_mass > spe_peaks[-1].mz:
            continue
        if start_mass < 0:
            start_mass = 0.0
        if end_mass > spe_peaks[-1].mz:
            end_mass = spe_peaks[-1].mz
        s_start_ind = peak_index[int(start_mass * multi)]
        # 确定肽段的这个碎片离子所对应谱图离子所在区间
        while s_start_ind < len(spe_peaks) and spe_peaks[s_start_ind].mz < start_mass:
            s_start_ind += 1
        if int(end_mass * multi) >= len(peak_index):
            s_end_ind = peak_index[-1]
        else:
            s_end_ind = peak_index[int(end_mass * multi)]
        while s_end_ind < len(spe_peaks) and spe_peaks[s_end_ind].mz <= end_mass:
            s_end_ind += 1
        if s_start_ind >= s_end_ind:
            continue
        one_mz_max_inten = 0.0
        one_mz_max_inten_index = -1
        # 找这里面强度最高的峰来计分
        for i in range(s_start_ind, s_end_ind):
            if spe_peaks[i].inten > one_mz_max_inten:
                one_mz_max_inten = spe_peaks[i].inten
                one_mz_max_inten_index = i
        # compute Score
        if one_mz_max_inten / spe_max_int < 0.001:
            continue
        if is_ppm_fra:
            one_mz_bias = (spe_peaks[one_mz_max_inten_index].mz - mz) * 1e6 / mz
        else:
            one_mz_bias = spe_peaks[one_mz_max_inten_index].mz - mz
        abs_one_mz_bias = abs(one_mz_bias) / max_bias
        # 获取离子偏差
        if abs_one_mz_bias > max_bias:
            continue
        else:
            match_num += 1
            list_ppm.append(abs_one_mz_bias)
            # 进到这里算匹配上
            # one_ionbias_score = math.sin(1.570796 * (1 - abs_one_mz_bias))
            # one_ionint_score = math.sin(1.570796 * (one_mz_max_inten / spe_max_int))
            one_ionbias_score = math.sin(0.785398 * (1 + 1 - abs_one_mz_bias))
            one_ionint_score = math.sin(0.785398 * (1 + one_mz_max_inten / spe_max_int))
            one_ion_score = math.sin(0.785398 * (1 + one_ionbias_score * one_ionint_score))

            match_ion_name.append([mz, one_ion_score])
            # 记录峰强度占谱图总强度的比例
            if one_mz_max_inten_index not in match_peak_index:
                match_inten += one_mz_max_inten
                match_peak_index.append(one_mz_max_inten_index)
                # one_ion_score = math.sin(0.785398 * (one_ionint_score + one_ionbias_score))
                # one_ion_score = math.sin(1.570796 * one_ionbias_score * one_ionint_score)
                match_ion_score += one_ion_score

            one_type_ion_list[-1] = 1# one_ion_score  # 匹配上了记录连续性的最后一位就由0变成1
    data_all_ioncontinue_score.append([last_ion_type, last_ion_charge, one_type_ion_list])

    # 算b离子连续分数
    continue_data = []
    one_ion_num = len(one_type_ion_list)
    all_type_tmp_list = [0 for i in range(one_ion_num)]
    peptide_coverage_list = [0 for i in range(one_ion_num)]

    for ion_type in ion_type_dic.keys():
        one_type_tmp_list = [0 for i in range(one_ion_num)]
        for ion_type_index in ion_type_dic[ion_type]:
            for ion_num_index in range(one_ion_num):
                if data_all_ioncontinue_score[ion_type_index][2][ion_num_index] == 1:
                    one_type_tmp_list[ion_num_index] = 1
                if ion_type == 'b':
                    if data_all_ioncontinue_score[ion_type_index][2][ion_num_index] == 1:
                        all_type_tmp_list[ion_num_index] = 1
                    if data_all_ioncontinue_score[ion_type_index][2][ion_num_index] > 0.001:
                        peptide_coverage_list[ion_num_index] = 1
                elif ion_type == 'y':
                    if data_all_ioncontinue_score[ion_type_index][2][ion_num_index] == 1:
                        all_type_tmp_list[one_ion_num - ion_num_index - 1] = 1
                    if data_all_ioncontinue_score[ion_type_index][2][ion_num_index] > 0.001:
                        peptide_coverage_list[one_ion_num - ion_num_index - 1] = 1
        continue_data.append(one_type_tmp_list)
    continue_data.append(all_type_tmp_list)
    tmp_one_type_match_score = 0.0
    for item in all_type_tmp_list:
        tmp_one_type_match_score += item
    out_continue_data = []

    TAG_LEN = 2
    for one_continue_data in continue_data:
        tmp_list = []
        # tmp_len_list = []
        tmp_continue_len = 0
        tmp_continue_score = 0.0
        all_contine_len = 0
        for one_continue_ion_score in one_continue_data:
            if one_continue_ion_score > 0.001:
                tmp_continue_len += 1
                tmp_continue_score += one_continue_ion_score
            else:
                if tmp_continue_len <= TAG_LEN:
                    pass
                else:
                    tmp_list.append(tmp_continue_score)
                    # tmp_len_list.append(tmp_continue_len)
                    all_contine_len += tmp_continue_len
                tmp_continue_score = 0.0
                tmp_continue_len = 0
        if tmp_continue_len <= TAG_LEN:
            pass
        else:
            tmp_list.append(tmp_continue_score)
            # tmp_len_list.append(tmp_continue_len)
            all_contine_len += tmp_continue_len

        if len(one_continue_data) == 0 or len(tmp_list) == 0:
            continue_score = 0.0
        else:
            tmp_continue_score = 0.0
            for tmp_index, one_continue_score in enumerate(tmp_list):
                tmp_continue_score += one_continue_score / TAG_LEN
            continue_score = tmp_continue_score
        out_continue_data.append(continue_score)

    tmp_coverage = 0
    for one_peptide_site in peptide_coverage_list:
        tmp_coverage += one_peptide_site
    peptide_coverage = tmp_coverage
    if spe_all_int == 0:
        match_inten_score = 0
    else:
        match_inten_score = match_inten / spe_all_int
    if len(spe_peaks) == 0:
        match_spe_score = 0
    else:
        match_spe_score = len(match_peak_index) / len(spe_peaks)

    sum = 0.0
    for item in out_continue_data:
        sum += item
    out_continue_data.append(sum / len(out_continue_data))
    # match_ion_name
    return CMatchOnePeptideScore(out_continue_data, match_ion_score, match_inten_score, match_spe_score, peptide_coverage, list_ppm, match_inten, match_ion_name)

def op_get_all_ion_match_score(peak_index, spe_peaks, spe_max_int, spe_all_int, ions, is_ppm_fra, tol_fra, multi):
    # peak_index是建一个0-谱图碎片离子最大质量的区间，1mz为一个区间，反映谱图碎片离子分布，之后匹配直接用这个索引看这里面的碎片离子进行匹配
    # multi是放缩谱图质荷比区间大小，默认是1mz
    match_num = 0
    match_inten = 0.0
    match_peak_index, match_ion_name, match_peak2ion = [], [], []
    match_ion_score = 0.0
    list_ppm = []

    # 这里是算误差允许范围，ppm的max_me就是20，max_me是0.02
    if is_ppm_fra:
        max_bias = tol_fra * 1e6
    else:
        max_bias = tol_fra
    # 遍历候选肽段的碎片离子

    for mi, mz in enumerate(ions):
        # 确定碎片允许的偏差范围

        if is_ppm_fra:
            tol_mass = mz * tol_fra
        else:
            tol_mass = tol_fra
        start_mass = mz - tol_mass
        end_mass = mz + tol_mass
        # 对这个偏差范围超出限制的就跳过下一个
        if end_mass <= 0:
            continue
        if start_mass > spe_peaks[-1].mz:
            continue
        if start_mass < 0:
            start_mass = 0.0
        if end_mass > spe_peaks[-1].mz:
            end_mass = spe_peaks[-1].mz
        s_start_ind = peak_index[int(start_mass * multi)]
        # 确定肽段的这个碎片离子所对应谱图离子所在区间
        while s_start_ind < len(spe_peaks) and spe_peaks[s_start_ind].mz < start_mass:
            s_start_ind += 1
        if int(end_mass * multi) >= len(peak_index):
            s_end_ind = peak_index[-1]
        else:
            s_end_ind = peak_index[int(end_mass * multi)]
        while s_end_ind < len(spe_peaks) and spe_peaks[s_end_ind].mz <= end_mass:
            s_end_ind += 1
        if s_start_ind >= s_end_ind:
            continue
        one_mz_max_inten = 0.0
        one_mz_max_inten_index = -1
        # 找这里面强度最高的峰来计分
        for i in range(s_start_ind, s_end_ind):
            if spe_peaks[i].inten > one_mz_max_inten:
                one_mz_max_inten = spe_peaks[i].inten
                one_mz_max_inten_index = i
        # compute Score
        if one_mz_max_inten / spe_max_int < 0.001:
            continue
        if is_ppm_fra:
            one_mz_bias = (spe_peaks[one_mz_max_inten_index].mz - mz) * 1e6 / mz
        else:
            one_mz_bias = spe_peaks[one_mz_max_inten_index].mz - mz
        abs_one_mz_bias = abs(one_mz_bias) / max_bias
        # 获取离子偏差
        if abs_one_mz_bias > max_bias:
            continue
        else:
            match_num += 1
            list_ppm.append(abs_one_mz_bias)
            # 进到这里算匹配上
            # one_ionbias_score = (math.sin(math.pi / 2 * (2 * abs(1 - abs_one_mz_bias) - 1)) + 1) / 2  # (1 - abs_one_me)  math.cos(math.pi * (1 - abs_one_me) / 2)
            # one_ionint_score = (math.sin(math.pi / 2 * (2 * one_mz_max_inten / spe_max_int - 1)) + 1) / 2
            # one_ionbias_score = math.sin(1.570796 * (1 - abs_one_mz_bias))
            # one_ionint_score = math.sin(1.570796 * (one_mz_max_inten / spe_max_int))
            one_ionbias_score = math.sin(0.785398 * (1 + 1 - abs_one_mz_bias))
            one_ionint_score = math.sin(0.785398 * (1 + one_mz_max_inten / spe_max_int))
            # 记录峰强度占谱图总强度的比例
            if one_mz_max_inten_index not in match_peak_index:
                match_inten += one_mz_max_inten
                match_peak_index.append(one_mz_max_inten_index)
                # one_ion_score = math.sin(0.785398 * (one_ionint_score + one_ionbias_score))
                # one_ion_score = math.sin(1.570796 * one_ionbias_score * one_ionint_score)
                one_ion_score = math.sin(0.785398 * (1 + one_ionbias_score * one_ionint_score))
                match_ion_score += one_ion_score
                match_ion_name.append([mz, one_ion_score])
    match_inten_percent = match_inten / spe_all_int
    match_spe_percent = len(match_peak_index) / len(spe_peaks)
    match_inten_score = match_inten_percent#len(match_peak_index) * match_inten_percent
    match_spe_score = match_spe_percent#len(match_peak_index) * match_spe_percent

    # match_ion_name
    return CMatchIonScore(match_ion_score, match_inten_score, match_spe_score, match_inten_percent, match_spe_percent, list_ppm, match_num, match_inten, spe_all_int, match_ion_name)
