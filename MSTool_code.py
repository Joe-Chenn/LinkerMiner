# -*- mode: python ; coding: utf-8 -*-..
# 这个考虑换C版本，python位运算天然慢

def tool_new_get_ion_code(pkl_file_index, ind_file_index, one_ind_num, link_site_combination_num):
    tmp_link_site_combination_num = link_site_combination_num
    tmp_one_ind_num = one_ind_num << 6
    tmp_ind_file_index = ind_file_index << 50
    tmp_pkl_file_index = pkl_file_index << 57
    out = tmp_link_site_combination_num | tmp_one_ind_num | tmp_ind_file_index | tmp_pkl_file_index
    return out


def tool_new_get_ion_data(input_data):
    link_site_combination = input_data & 0x3f
    one_ind_num = input_data >> 6 & 0x7ffffffffff
    ind_file_index = input_data >> 50 & 0x7f
    pkl_file_index = input_data >> 57 & 0x7f
    return pkl_file_index, ind_file_index, one_ind_num, link_site_combination


def tool_get_peptide_info_code(is_target, protein_index, peptide_start, peptide_length):
    tmp_peptide_length = peptide_length
    tmp_peptide_start = peptide_start << 21
    tmp_protein_index = protein_index << 42
    tmp_is_target = is_target << 63
    out = tmp_is_target | tmp_protein_index | tmp_peptide_start | tmp_peptide_length
    return out


def tool_get_peptide_info(input_data):
    peptide_length = input_data & 0x1fffff
    peptide_start = input_data >> 21 & 0x1fffff
    protein_index = input_data >> 42 & 0x1fffff
    is_target = input_data >> 63 & 1
    return is_target, protein_index, peptide_start, peptide_length


def tool_get_fix_mod_code(input_mod_site_list):
    fix_mod_code = 0
    for site, site_data in enumerate(input_mod_site_list):
        if site_data != -1:
            fix_mod_code = fix_mod_code | (1 << site)
    return fix_mod_code


def tool_get_fix_mod_site_list(fix_mod_code, mod_site_list, peptide_length):
    for site in range(peptide_length + 2):
        site_data = fix_mod_code & 1
        if site_data == 0:
            mod_site_list.append(-1)
        else:
            mod_site_list.append(0)
        fix_mod_code = fix_mod_code >> 1

def tool_update_var_mod_code_by_site_data(cur_var_mode_code, fix_site_value, cur_var_mod_num, site, site_data):

    if site_data >= fix_site_value:
        # 因为修饰索引是从0开始的，所以fix_site_value表示固定修饰数目,要取等于
        tmp_var_mod_code = (site << 6) | (site_data)  # 位点在前，修饰索引在后
        if cur_var_mod_num % 5 == 0:
            cur_var_mode_code = cur_var_mode_code | (tmp_var_mod_code << 48)
        else:
            bit_index = (cur_var_mod_num % 5 - 1) * 12
            cur_var_mode_code = cur_var_mode_code | (tmp_var_mod_code << bit_index)
    return cur_var_mode_code

def tool_get_var_mod_code(input_mod_site_list, input_var_mod_code_list, fix_site_value):
    cur_var_mod_num = 0
    cur_var_mode_code = 0
    for site, site_data in enumerate(input_mod_site_list):
        if site_data >= fix_site_value:
            # 因为修饰索引是从0开始的，所以fix_site_value表示固定修饰数目,要取等于
            cur_var_mod_num += 1
            tmp_var_mod_code = (site << 6) | (site_data)  # 位点在前，修饰索引在后
            if cur_var_mod_num % 5 == 0:
                cur_var_mode_code = cur_var_mode_code | (tmp_var_mod_code << 48)
                input_var_mod_code_list.append(cur_var_mode_code)
                cur_var_mode_code = 0
            else:
                bit_index = (cur_var_mod_num % 5 - 1) * 12
                cur_var_mode_code = cur_var_mode_code | (tmp_var_mod_code << bit_index)
    if cur_var_mode_code != 0:
        input_var_mod_code_list.append(cur_var_mode_code)

def tool_get_var_mod_data(input_var_mod_code_list, mod_site_list):

    for var_mod_code in input_var_mod_code_list:
        for i in range(5):
            one_var_mod_code = var_mod_code >> (12 * i) & 0xfff
            if one_var_mod_code == 0:
                break
            var_mod_index = one_var_mod_code & 0x3f
            var_mod_site = one_var_mod_code >> 6 & 0x3f
            mod_site_list[var_mod_site] = var_mod_index

