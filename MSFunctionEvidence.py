import numpy

from MSData import CEvidence
from MSLogging import logGetError
from MSOperator import opGetStartAndEndForProfile
from MSSysterm import VALUE_ILLEGAL
from MSTool import toolFindNeighborListFromSortedList1, toolFindNeighborFromSortedList1


# from scipy.signal import find_peaks
# from switch_plot import evidence_plot
# from scipy.interpolate import spline
# import torch


def toolSumList0(inputList, start, end):
    result = 0.0

    for i in range(start, end + 1):
        result = result + inputList[i]

    return result


def toolFindIndexFromSortedList1(inputList, number):
    return toolFindIndexFromSortedList0(inputList, 0, len(inputList) - 1, number)


def toolFindIndexFromSortedList0(inputList, start, end, number):
    if number < inputList[start]:
        return VALUE_ILLEGAL  # 没找到必须返回一个非法值

    if number > inputList[end]:
        return VALUE_ILLEGAL  # 没找到必须返回一个非法值

    index = soldierBinarySearch(inputList, start, end, number)

    if inputList[index] != number:

        return VALUE_ILLEGAL  # 没找到必须返回一个非法值

    else:

        return index


def soldierBinarySearch(inputList, start, end, number):
    # check
    if end >= len(inputList):
        logGetError("MSTool, MK311: The length of inputList is " + str(len(inputList)) + ", but you want to get " + str(
            end) + "?")

    if start < 0:
        logGetError(
            "MSTool, MK315: start is " + str(start) + "?")

    if number == inputList[end]:  # 最后一个有问题

        return end

    # find
    while start < end:
        mid = (start + end) // 2
        if number < inputList[mid]:
            end = mid
        elif number > inputList[mid]:
            start = mid + 1
        else:
            return mid

    start = start - 1

    return start


class CFunctionEvidence:

    def __init__(self, inputDP):
        self.dp = inputDP

    def __soldierGetPeaksFromPeaks(self, inputList1, inputList2, inputWinMOZ):

        '''
        根据list1去list2里面找最近的，两个list都是排好序的。返回list2的索引
        '''

        indexList = toolFindNeighborListFromSortedList1(inputList2, inputList1)

        result = [VALUE_ILLEGAL if inputWinMOZ < abs((i - inputList2[j]) / i * 1e6)
                  else j for i, j in zip(inputList1, indexList)]

        return result

    # DDALF & DDALL -> GetReference

    def __captainFillEvidenceGetReference(self, inputDataMS1,
                                          inputSeed, iMid, inputWinRT, flagGetStartAndEnd):

        # 确定保留时间的左界和右界
        rtMid = inputDataMS1.INDEX_RT[iMid]
        iLeft = iMid
        basicStride = 2 ** 7

        borderRTLeft = max(0, rtMid - inputWinRT * 60)  # 左界
        stride = basicStride
        # 找左右边界
        while iLeft > 0 and stride > 0:
            if inputDataMS1.INDEX_RT[max(iLeft - stride, 0)] < borderRTLeft or iLeft < stride:
                stride = int(stride / 2)
            else:
                iLeft -= stride

        iRight = iMid
        borderRTRight = min(rtMid + inputWinRT * 60, inputDataMS1.INDEX_RT[-1])  # 右界
        Max_iRight = len(inputDataMS1.INDEX_SCAN) - 1
        stride = basicStride
        while iRight < len(inputDataMS1.INDEX_SCAN) - 1 and stride > 0:

            if inputDataMS1.INDEX_RT[min(iRight + stride, Max_iRight)] > borderRTRight or iRight + stride > Max_iRight:
                stride = int(stride / 2)
            else:
                iRight += stride

        # init
        nScan = iRight - iLeft + 1
        outputEvidence = CEvidence()
        outputEvidence.MATRIX_PROFILE = numpy.zeros(
            [len(inputSeed.DIS_ISO_MOZ_CLC), nScan])  # nPeak行，nScan列, 搞为numpy是为了方便操作，不然面临list的传值传地址的问题，还有取行取列的问题
        outputEvidence.LIST_RET_TIME = [0] * nScan
        outputEvidence.LIST_SCAN = [0] * nScan

        for iScan in range(iLeft, iRight + 1):  # 这个iRight+1其实索引不到
            tmpScan = inputDataMS1.INDEX_SCAN[iScan]

            outputEvidence.LIST_SCAN[iScan - iLeft] = inputDataMS1.INDEX_SCAN[iScan]
            outputEvidence.LIST_RET_TIME[iScan - iLeft] = inputDataMS1.INDEX_RT[iScan]

            # fill
            listPeakMOZ = inputDataMS1.MATRIX_PEAK_MOZ[tmpScan]
            listPeakINT = inputDataMS1.MATRIX_PEAK_INT[tmpScan]

            indexPeak = self.__soldierGetPeaksFromPeaks(inputSeed.DIS_ISO_MOZ_CLC, listPeakMOZ, 20)

            for iPeak in range(len(indexPeak)):
                tmpIndex = indexPeak[iPeak]
                if tmpIndex == VALUE_ILLEGAL:
                    outputEvidence.MATRIX_PROFILE[iPeak, iScan - iLeft] = 0.0
                else:
                    outputEvidence.MATRIX_PROFILE[iPeak, iScan - iLeft] = listPeakINT[tmpIndex]  # iScan-iLeft不用加1

        # start and end
        if flagGetStartAndEnd:  # label free的sample不需要知道起点和终点
            nProfile = len(outputEvidence.MATRIX_PROFILE)
            outputEvidence.LIST_I_START = [0] * nProfile
            outputEvidence.LIST_I_END = [0] * nProfile

            outputEvidence.DIS_ISO_MOZ_EXP = [0] * nProfile
            outputEvidence.DIS_ISO_INT_EXP = [0] * nProfile

            iSeed = toolFindNeighborFromSortedList1(outputEvidence.LIST_SCAN, inputSeed.MID_SCAN)
            for iProfile in range(nProfile):
                listStartEnd = opGetStartAndEndForProfile(outputEvidence.MATRIX_PROFILE[iProfile, :], iSeed,
                                                          0.1,
                                                          # 0.3,
                                                          2)

                outputEvidence.LIST_I_START[iProfile] = listStartEnd[0]
                outputEvidence.LIST_I_END[iProfile] = listStartEnd[1]

                # isotopic distribution
                outputEvidence.DIS_ISO_INT_EXP[iProfile] = toolSumList0(outputEvidence.MATRIX_PROFILE[iProfile, ...],
                                                                        listStartEnd[0], listStartEnd[1])
                outputEvidence.DIS_ISO_MOZ_EXP[iProfile] = 0.0  # 需要根据强度作为权重求平均

            outputEvidence.I_START = outputEvidence.LIST_I_START[inputSeed.INDEX_MONO]
            outputEvidence.I_END = outputEvidence.LIST_I_END[inputSeed.INDEX_MONO]

            outputEvidence.PROFILE_ALL = [0]  # 这个需要进一步更新

            # if self.dp.myCFG.B11_FLAG_TEST_PLOT_XIC == CFG_FLAG_TEST_PLOT_XIC['Yes']:
            #
            #     switch_smooth_for_MIR = False
            #     plt.figure(figsize=(12, 3))
            #     x_RT = copy.deepcopy(outputEvidence.LIST_RET_TIME)
            #     y_profile = copy.deepcopy(outputEvidence.MATRIX_PROFILE)
            #
            #     # profile遍历加曲线
            #     if switch_smooth_for_MIR == False:
            #         # 不进行平滑处理 switch_smooth_for_MIR = False
            #         for i in range(len(y_profile)):
            #             plt.plot(x_RT, y_profile[i], linestyle='-',
            #                      label="{}".format(str(round(inputSeed.DIS_ISO_MOZ_CLC[i], 3))))
            #             plt.legend(loc="upper right")
            #     else:
            #         # 进行平滑处理 switch_smooth_for_MIR = True
            #         for i in range(len(y_profile)):
            #             x_smooth = numpy.linspace(min(x_RT), max(x_RT), 500)
            #             y_smooth = make_interp_spline(x_RT, y_profile[i])(x_smooth)
            #             for item in range(len(y_smooth)):  # 拒绝负值
            #                 if y_smooth[item] < 10:
            #                     y_smooth[item] = 0.0
            #             plt.plot(x_smooth, y_smooth, linestyle='-',
            #                      label="{}".format(str(round(inputSeed.DIS_ISO_MOZ_CLC[i], 3))))
            #             # plt.legend(loc="best")
            #             plt.legend(loc="upper right")
            #
            #     # 用于在画图限制标记的最高处
            #     matrix_int_T = numpy.transpose(outputEvidence.MATRIX_PROFILE).tolist()  # 转置 即每一个元素列表包含同一保留时间下的不同同位素的强度
            #     list_max = [max(i) for i in matrix_int_T]  # 选取每一个保留时间下的强度最高值
            #     index_rt_start = outputEvidence.LIST_I_START[inputSeed.INDEX_MONO]
            #     index_rt_end = outputEvidence.LIST_I_END[inputSeed.INDEX_MONO]
            #     rt_start = x_RT[index_rt_start]
            #     rt_end = x_RT[index_rt_end]
            #     inten_max = int(max(list_max))  # 目前用整个时间窗口内的最高强度值, 不用[start: end]内的最高强度值
            #
            #
            #     '''
            #     # 以下这部分是取超过相似性阈值的保留时间内的最高强度值
            #     # if index_rt_start != index_rt_end:
            #     #     inten_max = int(max(list_max[index_rt_start: index_rt_end]))
            #     # else:
            #     #     inten_max = int(list_max[index_rt_start])
            #     '''
            #
            #     # print("inten_max:{}".format(inten_max))
            #     # print("SCORE_IS_PEPTIDE:{}".format(outputEvidence.SCORE_IS_PEPTIDE))
            #
            #     # 如果再报错, 应该是这里出问题, 目前原因不明确 2022/03/08
            #     plt.plot([rt_start for i in range(int(0.9 * inten_max)) if i % 1000 == 0], [i for i in range(int(0.9 * inten_max)) if i % 1000 == 0], linestyle='--', color='#BEBEBE')
            #     plt.plot([rt_end for i in range(int(0.9 * inten_max)) if i % 1000 == 0], [i for i in range(int(0.9 * inten_max)) if i % 1000 == 0], linestyle='--', color='#BEBEBE')
            #     mid_rt = inputDataMS1.INDEX_RT[iMid]
            #     plt.plot([mid_rt for i in range(int(0.9 * inten_max)) if i % 1000 == 0], [i for i in range(int(0.9 * inten_max)) if i % 1000 == 0], linestyle='--', color='red')
            #     # plt.plot([rt_start for i in range(int(0.9 * inten_max)) if i % 10000 == 0], [i for i in range(int(0.9 * inten_max)) if i % 10000 == 0], linestyle='--', color='#BEBEBE')
            #     # plt.plot([rt_end for i in range(int(0.9 * inten_max)) if i % 10000 == 0], [i for i in range(int(0.9 * inten_max)) if i % 10000 == 0], linestyle='--', color='#BEBEBE')
            #
            #     # 添加transfer前后的糖肽组成label
            #     # index_ID_glycom = self.dp.myID.PSM_LINE_T[0].split("\t").index("GlycanComposition")  # 因为只有测试的时候才用得到, 所以放到这里了
            #     # ori_glycom = self.dp.myID.PSM_LINE_C[iPSM].split("\t")[index_ID_glycom]
            #     # new_glycom = self.__weaponGly2GlycoCom(self.dp.myID.PSM6_GLC[iPSM])
            #     # text_glycom_info = ori_glycom + "->" + new_glycom
            #
            #     # print("Transfer_Info: {}".format(text_glycom_info))
            #
            #     # 添加肽段、电荷数信息
            #     text_pep_charge_info = "RAW: " + str(self.dp.myID.PSM1_RAW_NAME[iPSM]) + "    " + \
            #                            "SCAN: " + str(self.dp.myID.PSM2_SCAN_ID[iPSM]) + "    " + \
            #                            "PEP: " + self.dp.myID.PSM4_SEQ[iPSM] + "    " + \
            #                            "CHARGE: " + str(self.dp.myID.PSM9_CHARGE[iPSM]) + "    " + \
            #                            "RT: " + str(round(rt_start, 3)) + "-" + str(round(rt_end, 3)) + "s    " + \
            #                            str(round(rt_start/60, 3)) + "-" + str(round(rt_end/60, 3)) + "min." + \
            #                            "GLC: " + self.__weaponGly2GlycoCom(self.dp.myID.PSM6_GLC[iPSM])
            #
            #     temp_int = 0.0
            #     indexMono = inputSeed.INDEX_MONO
            #     distribution = inputSeed.DIS_ISO_INT_CLC
            #     for iScan in range(outputEvidence.I_START, outputEvidence.I_END + 1):  # 最后一个还是需要加1的
            #         temp_int += outputEvidence.MATRIX_PROFILE[indexMono, iScan] / distribution[indexMono] * toolSumList1(distribution)  # 根据mono峰强度, 计算所有同位素峰的总强度
            #     farray = [str(int(x)) for x in inputSeed.DIS_ISO_INT_CLC]
            #     if str(self.dp.myID.PSM2_SCAN_ID[iPSM]) == '10961':
            #         print(1)
            #     text_pep_charge_info = "[SPECTRUM INFORMATION]\n\n" + \
            #                            "[RAW]: " + str(self.dp.myID.PSM1_RAW_NAME[iPSM]) + "\n" + \
            #                            "[SCAN]: " + str(self.dp.myID.PSM2_SCAN_ID[iPSM]) + "\n" + \
            #                            "[PEP]: " + self.dp.myID.PSM4_SEQ[iPSM] + "\n" + \
            #                            "[E]: " + str(self.dp.myID.PSM9_CHARGE[iPSM]) + "\n" + \
            #                            "[MONO]: " + str(indexMono) + "\n" + \
            #                            "[RT]: " + str(round(rt_start, 3)) + "-" + str(round(rt_end, 3)) + "s  " + str(round(rt_start / 60, 3)) + "-" + str(round(rt_end / 60, 3)) + "min.\n" + \
            #                            "[GLC]: " + self.__weaponGly2GlycoCom(self.dp.myID.PSM6_GLC[iPSM]) + "\n" + \
            #                            "[INT]: " + str("{:.2e}".format(temp_int)) + '\n' + '[CLC]: ' + str(','.join(farray))
            #     # print("[CLC_INT]" + str(' '.join(inputSeed.DIS_ISO_INT_CLC)))
            #
            #     # print()
            #
            #     # 添加text
            #     all_text = text_pep_charge_info
            #     text_x = outputEvidence.LIST_RET_TIME[0]
            #     text_y = inten_max * 0.95
            #     plt.text(text_x,
            #              text_y,
            #              all_text,
            #              fontsize=10,
            #              verticalalignment="top",
            #              horizontalalignment="left"
            #              )
            #
            #     # 输出和展示色谱曲线
            #     path_out_fig = self.dp.myCFG.E1_PATH_EXPORT + "fig"
            #     if not os.path.exists(path_out_fig):
            #         os.mkdir(path_out_fig)
            #
            #     plt.ticklabel_format(axis="y", style='sci', scilimits=(0, 0), useMathText=True)
            #     plt.rcParams['font.sans-serif'] = "Arial"
            #
            #     text_pep_charge_info = str(self.dp.myID.PSM1_RAW_NAME[iPSM]) + "_" + str(self.dp.myID.PSM2_SCAN_ID[iPSM]) + "_" + self.dp.myID.PSM4_SEQ[iPSM] + "_" + str(self.dp.myID.PSM9_CHARGE[iPSM]) + "_" + str(round(self.dp.myID.PSM3_RT[iPSM], 3))
            #     path_out_fig_1 = path_out_fig + r"\PSM{}_reference_{}.png".format(iPSM, text_pep_charge_info)
            #     # path_out_fig_2 = path_out_fig + r"\{}.pdf".format(text_pep_charge_info)
            #     plt.savefig(path_out_fig_1, bbox_inches ="tight")
            #     # plt.savefig(path_out_fig_2)
            #     # plt.show()
            #     plt.close()
        else:

            nProfile = len(outputEvidence.MATRIX_PROFILE)
            outputEvidence.LIST_I_START = [0] * nProfile
            outputEvidence.LIST_I_END = [0] * nProfile

            outputEvidence.I_START = VALUE_ILLEGAL
            outputEvidence.I_END = VALUE_ILLEGAL

            outputEvidence.DIS_ISO_MOZ_EXP = [0]
            outputEvidence.DIS_ISO_INT_EXP = [0]
            outputEvidence.PROFILE_ALL = [0]  # 这个需要进一步更新

        # output
        return outputEvidence

    def fillEvidenceGetReferenceForDDALF(self, inputDataMS1, inputSeed):  # DDA LF的reference用的

        iMid = toolFindIndexFromSortedList1(inputDataMS1.INDEX_SCAN, inputSeed.MID_SCAN)  # 理论上搞seed的时候，已经搞成ms1的了

        # print(inputDataMS1.INDEX_SCAN[-1])

        if inputDataMS1.INDEX_SCAN[iMid] != inputSeed.MID_SCAN:
            logGetError("MSFunctionEvidence, MK142: Can not find scan " + str(
                inputSeed.MID_SCAN) + " in inputDataMS1.INDEX_SCAN. inputDataMS1.INDEX_SCAN[iMid] is " + str(
                inputDataMS1.INDEX_SCAN[iMid]))
        # TODO 这种情况什么时候会发生? 是在前面使用了保留时间寻找一级谱的时候?

        winRT = 2

        # fill evidence

        return self.__captainFillEvidenceGetReference(inputDataMS1, inputSeed, iMid, winRT, True)  # DDALF
