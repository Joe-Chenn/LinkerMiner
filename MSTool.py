# -*- mode: python ; coding: utf-8 -*-..
import bisect
import copy
import re


def tool_compute_mass_from_element_for_isotope(inputMolecularformula, DIC_ELEMENT_MASS):
    tmp = re.split('\(|\)', inputMolecularformula)
    mass = 0.0
    for i in range(len(tmp)):
        if i % 2 != 0:
            continue
        if tmp[i].strip() == "":
            continue
        mass += DIC_ELEMENT_MASS[tmp[i]] * int(tmp[i + 1])
    return mass


def tool_compute_mass_from_element(inputMolecularformula, DIC_ELEMENT_MASS):
    tmp = re.split('\(|\)', inputMolecularformula)
    mass = 0.0
    for i in range(len(tmp)):
        if i % 2 != 0:
            continue
        if tmp[i].strip() == "":
            continue
        if tmp[i] not in DIC_ELEMENT_MASS:
            # if tmp[i][0]
            continue
        mass += (float(DIC_ELEMENT_MASS[tmp[i]]) * int(tmp[i + 1]))
    return mass


def tool_get_peptide_mass(peptide_seq, start_mass, mods, DIC_AA_MASS):  # 计算修饰后肽段质量
    # DIC_AA_MASS: mass of amino acid (A-Z)
    mass = start_mass
    for p in peptide_seq:
        if p < 'A' or p > 'Z': continue
        mass += DIC_AA_MASS[int(ord(p) - ord('A'))]
    for m in mods:
        mass += m.mass
    return mass


def tool_create_aa_list(sq, mods, aa2mass):
    # create the mass list of peptide sequence
    # e.g., peptide is "ACEDFK with C+57 modification"
    # return [M(A), M(C)+57, M(E), ..., M(K)] in which M(X) is the mass of X
    aa = [0.0 for i in range(len(sq))]
    for m in mods:
        aa[m.site - 1] += m.mass
    for i, p in enumerate(sq):
        if p < 'A' or p > 'Z':
            continue
        aa[i] += aa2mass[int(ord(p) - ord('A'))]
    return aa


def tool_binary_search(lis, key):
    low = 0
    high = len(lis) - 1
    while low < high:
        mid = int((low + high) / 2)
        if key == lis[low]:
            return mid
        if key == lis[high]:
            return high
        if key < lis[mid]:
            high = mid - 1
        elif key > lis[mid]:
            low = mid + 1
        else:
            return high
    return high


def tool_binary_search_index(input_list, find_key, start_pos=None, end_pos=None):
    if start_pos is None:
        start_pos = 0
    if end_pos is None:
        end_pos = len(input_list) - 1

    if len(input_list) == 0:
        return 0, 0

    if input_list[start_pos] > find_key:
        return -1, 0
    elif input_list[end_pos] < find_key:
        return end_pos + 1, end_pos + 2
    elif input_list[start_pos] == find_key:
        return 0, 1
    elif input_list[end_pos] == find_key:
        return end_pos, end_pos + 1
    else:
        while True:
            if (start_pos + end_pos) % 2 == 0:
                mid_pos = int((start_pos + end_pos) / 2)
            else:
                mid_pos = int((start_pos + end_pos - 1) / 2)
            if input_list[mid_pos] > find_key:
                end_pos = mid_pos
            elif input_list[mid_pos] < find_key:
                start_pos = mid_pos
            elif input_list[mid_pos] == find_key:
                start_pos = mid_pos
                end_pos = mid_pos + 1
                break
            if end_pos - start_pos == 1:
                if find_key < input_list[start_pos]:
                    start_pos = -1
                    end_pos = 0
                if find_key > input_list[end_pos]:
                    start_pos = end_pos
                    end_pos = end_pos + 1
                break
        return start_pos, end_pos


def tool_binary_search_left(input_list, find_key):
    start_pos = 0
    end_pos = len(input_list) - 1
    while start_pos < end_pos:
        if (start_pos + end_pos) % 2 == 0:
            mid_pos = int((start_pos + end_pos) / 2)
        else:
            mid_pos = int((start_pos + end_pos - 1) / 2)
        if input_list[mid_pos] < find_key:
            start_pos = mid_pos + 1
        elif input_list[mid_pos] > find_key:
            end_pos = mid_pos
        else:
            end_pos = mid_pos
    return end_pos


def tool_binary_search_right(input_list, find_key):
    start_pos = 0
    end_pos = len(input_list) - 1
    while start_pos < end_pos:
        if (start_pos + end_pos) % 2 == 0:
            mid_pos = int((start_pos + end_pos) / 2)
        else:
            mid_pos = int((start_pos + end_pos - 1) / 2)
        if input_list[mid_pos] < find_key:
            start_pos = mid_pos + 1
        elif input_list[mid_pos] > find_key:
            end_pos = mid_pos
        else:
            start_pos = mid_pos + 1
    return start_pos - 1


def tool_binary_search_keynotsamelist(lis, key):
    low = 0
    high = len(lis) - 1
    if key > lis[high]:
        low = high
    else:
        while low < high:
            mid = int((low + high) / 2)
            if key < lis[mid]:
                high = mid
            else:
                low = mid
            if (high - low) == 1:
                break
    return low


def tool_generate_singlepepdite_by_ion(aa2, ion_charge_max=2, mass_proton=1.00727645224, b_start_mass=0.0,
                                       y_start_mass=18.0105647):
    b_mz = [[] for i in range(ion_charge_max)]
    y_mz = [[] for i in range(ion_charge_max)]

    tmp_mass = 0.0
    for i in range(len(aa2) - 1):
        tmp_mass += aa2[i]
        for j in range(ion_charge_max):
            b_mz[j].append((tmp_mass + b_start_mass + mass_proton * (j + 1)) / (j + 1))

    tmp_mass = 0.0
    for i in range(1, len(aa2))[::-1]:
        tmp_mass += aa2[i]
        for j in range(ion_charge_max):
            y_mz[j].append((tmp_mass + mass_proton * (j + 1) + y_start_mass) / (j + 1))

    out_b_mz = []
    out_y_mz = []
    for i in range(len(b_mz)):
        out_b_mz += b_mz[i]
    for i in range(len(y_mz)):
        out_y_mz += y_mz[i]

    return out_b_mz, out_y_mz


def tool_generate_loop_type_by_ion(aa2, link_site_1, link_site_2, loop_mass, ion_charge_max=2, mono_site=0, mono_mass=0,
                                   mass_proton=1.00727645224, b_start_mass=0.0, y_start_mass=18.0105647):
    b_mz = [[] for i in range(ion_charge_max)]
    y_mz = [[] for i in range(ion_charge_max)]

    tmp_mass = 0.0
    for i in range(len(aa2) - 1):

        if (i + 1) == mono_site:
            tmp_mass += aa2[i] + mono_mass
        else:
            tmp_mass += aa2[i]

        if i >= link_site_1 - 1 and i < link_site_2 - 1:
            if i == link_site_1:
                tmp_mass += loop_mass
        else:
            for j in range(ion_charge_max):
                b_mz[j].append((tmp_mass + mass_proton * (j + 1) + b_start_mass) / (j + 1))

    tmp_mass = 0.0
    for i in range(1, len(aa2))[::-1]:

        if (i + 1) == mono_site:
            tmp_mass += aa2[i] + mono_mass
        else:
            tmp_mass += aa2[i]
        if i >= link_site_1 - 1 and i <= link_site_2 - 1:
            if i == link_site_2:
                tmp_mass += loop_mass
        else:
            for j in range(ion_charge_max):
                y_mz[j].append((tmp_mass + mass_proton * (j + 1) + y_start_mass) / (j + 1))

    out_b = []
    out_y = []
    for item in b_mz:
        out_b = out_b + item
    for item in y_mz:
        out_y = out_y + item

    return out_b, out_y


def tool_generate_only_cross_by_ion(aa2, ion_charge_max=2, mono_site=0, mono_mass=0, mass_proton=1.00727645224,
                                    b_start_mass=0.0, y_start_mass=18.0105647):
    b_mz = [[] for i in range(ion_charge_max)]
    y_mz = [[] for i in range(ion_charge_max)]

    tmp_mass = 0.0
    for i in range(len(aa2) - 1):
        if (i + 1) == mono_site:
            tmp_mass += aa2[i] + mono_mass
        else:
            tmp_mass += aa2[i]
        for j in range(ion_charge_max):
            b_mz[j].append((tmp_mass + mass_proton * (j + 1) + b_start_mass) / (j + 1))

    tmp_mass = 0.0
    for i in range(1, len(aa2))[::-1]:
        if (i + 1) == mono_site:
            tmp_mass += aa2[i] + mono_mass
        else:
            tmp_mass += aa2[i]
        for j in range(ion_charge_max):
            y_mz[j].append((tmp_mass + mass_proton * (j + 1) + y_start_mass) / (j + 1))
    out_b = []
    out_y = []
    for item in b_mz:
        out_b = out_b + item
    for item in y_mz:
        out_y = out_y + item

    return out_b, out_y


def tool_generate_cross_cross_by_ion(aa2, pep_link_site, cross_mass, ion_charge_max=2, mass_proton=1.00727645224,
                                     b_start_mass=0.0, y_start_mass=18.0105647):
    b_mz = [[] for i in range(ion_charge_max)]
    b_flag = [[] for i in range(ion_charge_max)]
    y_mz = [[] for i in range(ion_charge_max)]
    y_flag = [[] for i in range(ion_charge_max)]

    tmp_mass = 0.0
    for i in range(len(aa2) - 1):

        if (i + 1) == pep_link_site[1]:
            tmp_mass += aa2[i] + cross_mass
        else:
            tmp_mass += aa2[i]
        if (i + 1) >= pep_link_site[0] and (i + 1) < pep_link_site[1]:
            continue
        for j in range(ion_charge_max):
            b_mz[j].append((tmp_mass + b_start_mass + mass_proton * (j + 1)) / (j + 1))

    tmp_mass = 0.0
    for i in range(1, len(aa2))[::-1]:
        if (i + 1) == pep_link_site[0]:
            tmp_mass += aa2[i] + cross_mass
        else:
            tmp_mass += aa2[i]
        if (i + 1) <= pep_link_site[1] and (i + 1) > pep_link_site[0]:
            continue
        for j in range(ion_charge_max):
            y_mz[j].append((tmp_mass + mass_proton * (j + 1) + y_start_mass) / (j + 1))

    out_b = []
    out_y = []
    for item in b_mz:
        out_b = out_b + item
    for item in y_mz:
        out_y = out_y + item

    return out_b, out_y


def tool_generate_cross_loop_by_ion(aa2, link_site_1, link_site_2, loop_mass, cross_site, cross_mass, ion_charge_max=2,
                                    mono_site=0, mono_mass=0, mass_proton=1.00727645224, b_start_mass=0.0,
                                    y_start_mass=18.0105647):
    b_mz = [[] for i in range(ion_charge_max)]
    y_mz = [[] for i in range(ion_charge_max)]

    tmp_mass = 0.0
    for i in range(0, len(aa2) - 1):

        if (i + 1) == mono_site:
            tmp_mass += aa2[i] + mono_mass
        else:
            tmp_mass += aa2[i]

        if (i + 1) == cross_site:
            tmp_mass += cross_mass

        if i + 1 >= link_site_1 and i + 1 < link_site_2:
            if i == link_site_1:
                tmp_mass += loop_mass
        else:
            for j in range(ion_charge_max):
                b_mz[j].append((tmp_mass + mass_proton * (j + 1) + b_start_mass) / (j + 1))

    tmp_mass = 0.0
    for i in range(1, len(aa2))[::-1]:

        if (i + 1) == mono_site:
            tmp_mass += aa2[i] + mono_mass
        else:
            tmp_mass += aa2[i]

        if (i + 1) == cross_site:
            tmp_mass += cross_mass

        if i + 1 > link_site_1 and i + 1 <= link_site_2:
            if i + 1 == link_site_2:
                tmp_mass += loop_mass
        else:
            for j in range(ion_charge_max):
                y_mz[j].append((tmp_mass + mass_proton * (j + 1) + y_start_mass) / (j + 1))

    out_b = []
    out_y = []
    for item in b_mz:
        out_b = out_b + item
    for item in y_mz:
        out_y = out_y + item

    return out_b, out_y


def tool_set_one_list(input_list, input_list_same_num, max_bias, out_list, out_same_num):
    if len(input_list) == 0:
        pass
    elif len(input_list) == 1:
        out_list.append(input_list[0])
    else:
        left_data = input_list[0]
        out_list.append(left_data)
        out_same_num.append(input_list_same_num[0])
        for right_index in range(1, len(input_list), 1):
            right_data = input_list[right_index]
            if right_data - left_data < max_bias:
                out_same_num[-1] += input_list_same_num[right_index]
            else:
                out_list.append(right_data)
                left_data = right_data
                out_same_num.append(input_list_same_num[right_index])


def tool_set_two_list(input_list_alpha, input_list_beta, max_bias, out_list, out_same_num, input_list_alpha_num=None,
                      input_list_beta_num=None, out_only_alpha=None, out_only_beta=None):
    # 输入是两个有序的, out_list是把两个list合并之后的结果，out_same_num是对应每个值的重复个数
    if out_only_alpha is None:
        out_only_alpha = []
    if out_only_beta is None:
        out_only_beta = []

    if input_list_alpha_num is None:
        input_list_alpha_num = [1 for i in range(len(input_list_alpha))]
    if input_list_beta_num is None:
        input_list_beta_num = [1 for i in range(len(input_list_beta))]
    if len(input_list_alpha) == 0 and len(input_list_beta) == 0:
        pass
    else:
        alpha_ion_index = 0
        beta_ion_index = 0
        if len(input_list_alpha) == 0:
            last_data = input_list_beta[beta_ion_index]
            beta_ion_index = 1
            out_only_beta.append(input_list_beta[beta_ion_index])
            out_same_num.append(1)
        elif len(input_list_beta) == 0:
            last_data = input_list_alpha[alpha_ion_index]
            alpha_ion_index = 1
            out_only_alpha.append(input_list_alpha[alpha_ion_index])
            out_same_num.append(1)
        else:
            if abs(input_list_alpha[alpha_ion_index] - input_list_beta[beta_ion_index]) <= max_bias:
                last_data = input_list_alpha[alpha_ion_index]
                alpha_ion_index += 1
                beta_ion_index += 1
                out_same_num.append(2)
            else:
                if input_list_alpha[alpha_ion_index] > input_list_beta[beta_ion_index]:
                    last_data = input_list_beta[beta_ion_index]
                    beta_ion_index += 1
                    out_only_beta.append(input_list_beta[beta_ion_index])
                    out_same_num.append(1)
                else:
                    last_data = input_list_alpha[alpha_ion_index]
                    alpha_ion_index += 1
                    out_only_alpha.append(input_list_alpha[alpha_ion_index])
                    out_same_num.append(1)
        out_list.append(last_data)
        while True:
            if alpha_ion_index == len(input_list_alpha) and beta_ion_index == len(input_list_beta):
                break
            if alpha_ion_index == len(input_list_alpha):

                while beta_ion_index < len(input_list_beta):

                    cur_data = input_list_beta[beta_ion_index]

                    if abs(cur_data - last_data) < max_bias:
                        out_same_num[-1] += input_list_beta_num[beta_ion_index]
                    else:
                        out_list.append(cur_data)
                        out_same_num.append(input_list_beta_num[beta_ion_index])
                        out_only_beta.append(input_list_beta[beta_ion_index])

                    beta_ion_index += 1
                    last_data = cur_data

            elif beta_ion_index == len(input_list_beta):

                while alpha_ion_index < len(input_list_alpha):

                    cur_data = input_list_alpha[alpha_ion_index]

                    if abs(cur_data - last_data) < max_bias:
                        out_same_num[-1] += input_list_alpha_num[alpha_ion_index]
                    else:
                        out_list.append(cur_data)
                        out_same_num.append(input_list_alpha_num[alpha_ion_index])
                        out_only_alpha.append(input_list_alpha[alpha_ion_index])

                    alpha_ion_index += 1
                    last_data = cur_data

            else:
                # 先判断两个是否一样
                if abs(input_list_alpha[alpha_ion_index] - input_list_beta[beta_ion_index]) < max_bias:
                    cur_data = (input_list_alpha[alpha_ion_index] + input_list_beta[beta_ion_index]) / 2
                    # 再判断和列表最后一个是否一样
                    if abs(cur_data - last_data) < max_bias:
                        out_same_num[-1] += input_list_alpha_num[alpha_ion_index] + input_list_beta_num[beta_ion_index]
                    else:
                        out_list.append(cur_data)
                        out_same_num.append(input_list_alpha_num[alpha_ion_index] + input_list_beta_num[beta_ion_index])
                    alpha_ion_index += 1
                    beta_ion_index += 1
                else:
                    if input_list_alpha[alpha_ion_index] > input_list_beta[beta_ion_index]:
                        cur_data = input_list_beta[beta_ion_index]
                        if abs(cur_data - last_data) < max_bias:
                            out_same_num[-1] += input_list_beta_num[beta_ion_index]
                        else:
                            out_list.append(cur_data)
                            out_same_num.append(input_list_beta_num[beta_ion_index])
                            out_only_beta.append(input_list_beta[beta_ion_index])
                        beta_ion_index += 1
                    else:
                        cur_data = input_list_alpha[alpha_ion_index]
                        if abs(cur_data - last_data) < max_bias:
                            out_same_num[-1] += input_list_alpha_num[alpha_ion_index]
                        else:
                            out_list.append(cur_data)
                            out_same_num.append(input_list_alpha_num[alpha_ion_index])
                            out_only_alpha.append(input_list_alpha[alpha_ion_index])
                        alpha_ion_index += 1
                last_data = cur_data


def tool_get_N_combine(input_data, out_data, get_num):
    # 这个是产生从M数组里取小于等于N个数的所有组合
    if get_num > len(input_data):
        get_num = len(input_data)
    for one_get_num in range(1, get_num + 1, 1):
        tool_combine_list_N_point(one_get_num, 0, len(input_data), input_data, out_data)


def tool_combine_list_N_point(get_num, start, end, data_list, out):
    # 对一个数组，获取N个点的构成的所有子集集合
    for i in range(start, end - get_num + 1, 1):
        cur_site = []
        tool_combine_fixstart_point(get_num, i, end, data_list, [], cur_site)
        out += cur_site


def tool_combine_fixstart_point(get_num, start, end, data_list, cur_site, out_site):
    # 对一个数组，指定起始位置点，获取N个点组合的子集
    if get_num == 1:
        cur_site.append(data_list[start])
        out_site.append(cur_site)
    else:
        cur_site.append(data_list[start])
        for i in range(start + 1, end):
            new_site = copy.copy(cur_site)
            tool_combine_fixstart_point(get_num - 1, i, end, data_list, new_site, out_site)


def toolCountCharInString(inputStr, inputChar):
    result = 0

    for c in inputStr:
        if c == inputChar:
            result = result + 1

    return result


def toolGetWord(inputString, index, d):
    if inputString[0] != d:
        inputString = d + inputString

    if inputString[-1] != d:
        inputString = inputString + d

    p_d = []

    i = 0
    for c in inputString:

        if c == d:
            p_d.append(i)

        i = i + 1

    result = inputString[p_d[index] + 1:p_d[index + 1]]

    return result


def toolFindNeighborFromSortedList1(inputSortedList, number):
    start, end = 0, len(inputSortedList) - 1

    if number <= inputSortedList[start]:
        return start

    if number >= inputSortedList[end]:
        return end

    neighbor = bisect.bisect_left(inputSortedList, number)
    if abs(inputSortedList[neighbor] - number) < abs(inputSortedList[neighbor - 1] - number):
        return neighbor
    else:
        return neighbor - 1


def toolFindFromListByKey(inputList, key):
    n = len(inputList)

    for i in range(len(inputList)):

        tmpELE = inputList[i]

        if -1 == tmpELE.find(key):
            result = -1
        else:
            result = i
            return result

    return n + 1


def toolFindNeighborListFromSortedList1(inputSortedList, numberList):
    flag = 0
    pointer = 0
    pointerNext = len(numberList) - 1
    list_border = [0, len(inputSortedList) - 1]
    list_result = [0] * len(numberList)

    for _ in range(len(numberList)):
        start, end = list_border
        number = numberList[pointer]

        if number <= inputSortedList[start]:
            neighbor = start

        elif number >= inputSortedList[end]:
            neighbor = end

        else:
            neighbor = bisect.bisect_left(inputSortedList, number, start, end)

            neighbor = neighbor if abs(inputSortedList[neighbor] - number) < \
                                   abs(inputSortedList[neighbor - 1] - number) else neighbor - 1

        list_border[flag] = neighbor
        list_result[pointer] = neighbor
        # 找下一个numberList的点
        next = pointer + 1 if pointer < pointerNext else pointer - 1
        pointer = pointerNext
        pointerNext = next

        flag = 0 if flag == 1 else 1

    return list_result


import math
from typing import List
from tqdm import tqdm
import numpy as np
import time


class timer:
    def __init__(self):
        self.start = time.time()
        self.end = time.time()

    def reset(self):
        self.start = time.time()
        self.end = time.time()

    def elapsed(self):
        self.end = time.time()
        return round(self.end - self.start, 2)


def flatten_and_dump(array: List[List], dump_path: str):
    flattened_list = []
    prefix = []

    for sublist in array:
        prefix.append(len(flattened_list))
        flattened_list.extend(sublist)

    prefix.append(len(flattened_list))

    flattened_list, prefix = np.array(flattened_list, dtype=np.int32), np.array(prefix, dtype=np.int32)
    np.savez(dump_path, flattened_list=flattened_list, prefix=prefix)
    print("Dump file to {}".format(dump_path))


def read_file(fileName: str, tolerance: int):
    data = []
    with open(fileName, 'r') as f:
        one_data = []
        for line in tqdm(f):
            if line.startswith("#"):
                if len(one_data) != 0:
                    data.append(one_data)
                one_data = []
            else:
                one_data.append(round(float(line.strip()) * tolerance))
    dump_path = fileName.split('.')[0] + '.npz'
    flatten_and_dump(data, dump_path)


def get_max_item(array):
    return np.max(array)


def gen_ivf(array, prefix, path):
    max_item = math.ceil(get_max_item(array)) + 1
    # 初始化大小为max_item的ivf数组
    ivf = [[] for _ in range(max_item)]
    for i in tqdm(range(len(prefix) - 1)):
        for j in range(prefix[i], prefix[i + 1]):
            ivf[array[j]].append(i)

    dump_path = path.split('.')[0] + '_dic.npz'
    flatten_and_dump(ivf, dump_path)
    return max_item


def get_search_arrays_idx(a):
    # Step 1: Compute the differences
    counts = np.diff(a, prepend=a[0])

    # Step 2: Compute the total size of the expanded vector
    total = a[-1]

    # Step 3: Generate the expanded index values
    indices = np.zeros(total, dtype=int)
    np.cumsum(counts[:-1], out=indices[1:])

    # Step 4: Scatter the indices based on the counts
    b = np.zeros(total, dtype=int)
    for i in range(len(counts)):
        if counts[i] > 0:
            b[indices[i]:indices[i] + counts[i]] = i

    return b


def expand(array, tolerance: int):
    array = (array * tolerance)
    return array.astype(int)