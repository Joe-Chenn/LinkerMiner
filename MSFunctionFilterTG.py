import os
import pickle
import time
from typing import List

from sklearn.preprocessing import MinMaxScaler

from MSData import CFileMS1, CDataPack, CINI, CSeed, CEvidence
from MSDataResult import CWriteResultOnlyCross
from MSEmass import CEmass
from MSFunction import CFunctionPickle
from MSFunctionComposition import CFunctionComposition
from MSFunctionEvidence import CFunctionEvidence
from MSLogging import logGetError, logToUser
from MSSysterm import FILE_NAME_PSM_RESULT
from MSTool import toolFindNeighborFromSortedList1, toolCountCharInString

import numpy as np
from sklearn.metrics.pairwise import cosine_similarity, euclidean_distances


class CFunctionReadMS1:
    def __init__(self, inputDP: CDataPack):
        self.dp = inputDP

    def read(self, pathMS1):
        # check
        path_pkl = pathMS1 + ".pkl"
        if os.access(path_pkl, os.F_OK):
            pass
        else:
            dataMS1 = CFileMS1.newInstance()
            # open
            with open(pathMS1, 'r') as f:
                i_MS1 = -1
                for line in f.readlines():
                    len_line = len(line)
                    if len_line > 1:
                        if line.startswith("S	"):
                            i_MS1 = i_MS1 + 1
                            tmpScan = int(self.toolGetWord(line, 1, '	'))
                            dataMS1.INDEX_SCAN.append(tmpScan)
                            dataMS1.MATRIX_PEAK_MOZ[tmpScan] = []  # 这两行不能少
                            dataMS1.MATRIX_PEAK_INT[tmpScan] = []
                        elif line.startswith("H	"):
                            pass
                        elif line.startswith("I	"):
                            if line.startswith("I	IonInjectionTime"):
                                t = self.toolGetWord(line, 2, '	')
                                dataMS1.LIST_ION_INJECTION_TIME[tmpScan] = float(t)
                            elif line.startswith("I	RetTime	"):
                                t = self.toolGetWord(line, 2, '	')
                                dataMS1.INDEX_RT.append(float(t))
                        else:
                            dataMS1.MATRIX_PEAK_MOZ[tmpScan].append(float(self.toolGetWord(line, 0, ' ')))
                            dataMS1.MATRIX_PEAK_INT[tmpScan].append(float(self.toolGetWord(line, 1, ' ')))

            # write pkl
            fid_pkl = open(path_pkl, 'wb')
            pickle.dump(dataMS1, fid_pkl)
            fid_pkl.close()

    @staticmethod
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


class CFunctionFilterTG:
    def __init__(self, inputDP: CDataPack):
        self.dp = inputDP

    def filter(self, search_result: List[CWriteResultOnlyCross]):
        nPSM = len(search_result)

        nMS1 = len(self.dp.LIST_PATH_MS1)
        filtered_result = []
        for iMS1 in range(nMS1):
            pathMS1 = self.dp.LIST_PATH_MS1[iMS1]
            nameRawMS1 = self.toolGetNameFromPath(pathMS1)

            # load
            time0 = time.perf_counter()

            pathPKL = pathMS1 + ".pkl"
            pklFile = open(pathPKL, 'rb')
            dataMS1 = pickle.load(pklFile)
            pklFile.close()

            time1 = time.perf_counter()
            print("loadMS1Data: {}".format(time1 - time0))

            for index, one_res in enumerate(search_result):
                if one_res.spectrum_title.split('.')[0] != nameRawMS1:
                    continue
                canBeSingle, singleSeq, singleMod = self.getSingle(one_res.peptide_sq, one_res.peptide_modification,
                                                                   one_res.peptide_type, one_res.peptide_protein)
                if canBeSingle == False:
                    continue

                ## 1. 计算单肽&交联的同位素分布
                psmComposition = self.__getStrComposition(one_res.peptide_sq,
                                                          one_res.peptide_modification,
                                                          one_res.linker,
                                                          self.dp.myINI)

                psmCompositionSingle = self.__getStrComposition(singleSeq, singleMod, "", self.dp.myINI)

                # 得到种子
                tmpSeed = self.__captainFillSeed(int(one_res.scan), int(one_res.spectrum_charge), psmComposition,
                                                 dataMS1.INDEX_SCAN, dataMS1.INDEX_RT)
                singleTmpSeed = self.__captainFillSeed(int(one_res.scan), int(one_res.spectrum_charge), psmCompositionSingle, dataMS1.INDEX_SCAN,
                                                       dataMS1.INDEX_SCAN)

                # 得到evidence
                functionEvidence = CFunctionEvidence(self.dp)

                # if self.dp.myID.PSM2_SCAN_ID[iPSM] == 6091 or self.dp.myID.PSM2_SCAN_ID[iPSM] == 6276:
                # 	print(111)

                tmpEvidence = functionEvidence.fillEvidenceGetReferenceForDDALF(dataMS1,
                                                                                tmpSeed)  # __captainGetReference
                singleTmpEvidence = functionEvidence.fillEvidenceGetReferenceForDDALF(dataMS1,
                                                                                  singleTmpSeed)


                # formatted_arr = [f"{x:.2f}" for x in tmpSeed.DIS_ISO_MOZ_CLC]
                # print("cross:{}\tmass:[{}]".format(one_res.peptide_sq, ', '.join(formatted_arr)))
                crossSimi = self.__calc_similarity(tmpEvidence, tmpSeed)
                # formatted_arr = [f"{x:.2f}" for x in singleTmpSeed.DIS_ISO_MOZ_CLC]
                # print("single:{}\tmass:[{}]".format(singleSeq, ', '.join(formatted_arr)))
                singleSimi = self.__calc_similarity(singleTmpEvidence, singleTmpSeed)


                if crossSimi > singleSimi:
                    filtered_result.append(one_res)

        result_folder = self.dp.myCFG.A5_PATH_RESULT_EXPORT
        result_folder_pkl = result_folder + '\\tmp\\'
        write_result_path = result_folder + "\\" + FILE_NAME_PSM_RESULT[2][1] + '.txt'
        write_result_pkl_path = result_folder_pkl + FILE_NAME_PSM_RESULT[2][1] + '.pkl'
        self.__captain_write_onlycross_PSM(filtered_result, write_result_pkl_path, write_result_path)

        return filtered_result

    def getSingle(self, inputSeq: str, inputMod: str, inputPepType: str, inputPro: str):
        if inputPepType != "cross":
            return False, "", ""
        else:
            peps = inputSeq.split('-')
            pro = inputPro.split('\\')
            mods = inputMod.split(';')

            outputSeq = ""
            outputMod = ""
            flag = False

            for pros in pro[:-1]:
                pros = pros.split('-')

                pep_index_1 = int(peps[0].split('(')[1][:-1])
                pep_index_2 = int(peps[1].split('(')[1][:-1])
                pro_index_1 = int(pros[0].split('(')[1][:-1])
                pro_index_2 = int(pros[1].split('(')[1][:-1])

                if pro_index_1 < pro_index_2:
                    length = len(peps[0].split('(')[0]) - pep_index_1 + pep_index_2
                    # 把str=A(xx)-B(yy)变成str=AB的形式
                    if length == pro_index_2 - pro_index_1:
                        flag = True
                        outputSeq += peps[0].split("(")[0]
                        outputSeq += peps[1].split("(")[0]

                        for mod in mods:
                            if mod == "-" or mod == "":
                                continue
                            index, modificationName = mod.split(":")
                            index = int(index)
                            if index < 0:
                                index = index * -1
                                index -= 1
                                index += len(peps[0].split('(')[0])
                            outputMod += "{}:{};".format(index, modificationName)
                        break
                else:
                    length = len(peps[1].split('(')[0]) - pep_index_2 + pep_index_1
                    if length == pro_index_1 - pro_index_2:
                        flag = True
                        outputSeq += peps[1].split("(")[0]
                        outputSeq += peps[0].split("(")[0]

                        for mod in mods:
                            if mod == "-" or mod == "":
                                continue
                            index, modificationName = mod.split(":")
                            index = int(index)
                            if index < 0:
                                index = index * -1
                            else:
                                index += len(peps[1].split('(')[0])
                            outputMod += "{}:{};".format(index, modificationName)
                        break

            return flag, outputSeq, outputMod

    def __getStrComposition(self, inputSeq, inputMod, inputLIK: str, inputINI: CINI):
        result = 'H(2)O(1)'  # 写死的

        for aa in inputSeq:
            if inputINI.DIC_AA_COMPOSITION.__contains__(aa):
                result = result + inputINI.DIC_AA_COMPOSITION[aa]
            else:
                continue

        # modification
        nMOD = inputMod.split(";")  # 标准形式为6:Carbamidomethyl[C];14:Carbamidomethyl[C];

        for iMOD in nMOD:
            if iMOD == "-" or iMOD == "":
                continue
            nameMOD = iMOD.split(':')[1]
            result = result + inputINI.DIC_MOD[nameMOD].composition


        try:
            if not inputLIK == "":
                result = result + inputINI.DIC_LINKER[inputLIK].composition
            if inputSeq.find('-') != -1:
                result += 'H(2)O(1)'
        except:
            print('error in getStrComposition, unknown linker: [{}]'.format(inputLIK))

        return result

    def toolGetNameFromPath(self, path: str):

        lenStr = len(path)
        iStart = 0
        iEnd = -1

        for i in range(lenStr):

            j = lenStr - 1 - i

            if path[j] == '.':
                iEnd = j  # 亲测必须这么写，不用减一

            if path[j] == '\\' or path[j] == '/':
                iStart = j + 1
                break

        return path[iStart:iEnd]

    def __captainFillSeed(self, scanNum, charge, compositionPSM, listMS1Scan, listMS1RT):  # 和LL最大的不同是，这里不用考虑标记信息

        mySeed = CSeed()

        indexMS1 = toolFindNeighborFromSortedList1(listMS1Scan, scanNum)
        mySeed.MID_SCAN = listMS1Scan[indexMS1]
        mySeed.MID_RT = listMS1RT[indexMS1]

        # 元素组成
        fComposition = CFunctionComposition(self.dp)
        mySeed.DICT_COMPOSITION = fComposition.getDictComposition(compositionPSM)

        # 调emass，得到质量分布；
        emass = CEmass(self.dp)
        # try:
        tmpIsoDis = emass.getCalculatedIsotopicPeaks(mySeed.DICT_COMPOSITION, charge)
        # except:
        #     logGetError("Wrong element!")

        mySeed.DIS_ISO_INT_CLC = tmpIsoDis[0]
        mySeed.DIS_ISO_MOZ_CLC = tmpIsoDis[1]

        # 更新下mono
        tmpMOZMono = emass.getCalculatedMonoMZ(mySeed.DICT_COMPOSITION, charge)
        mySeed.INDEX_MONO = toolFindNeighborFromSortedList1(mySeed.DIS_ISO_MOZ_CLC, tmpMOZMono)

        return mySeed

    def __calc_similarity(self, evidence: CEvidence, seed: CSeed):
        clc_tensor = np.array(seed.DIS_ISO_INT_CLC)

        def normalize_to_100(arr):
            max_val = np.max(arr)
            if max_val == 0:  # 防止除以0的情况
                return arr
            return (arr / max_val) * 100.0

        exp_tensor = normalize_to_100(np.array(evidence.DIS_ISO_INT_EXP))


        if len(clc_tensor) != len(exp_tensor):
            return 0
        similarity = cosine_similarity(clc_tensor.reshape(1, -1), exp_tensor.reshape(1, -1))
        # print("clc:{}, ecp:{}, simi:{}".format(np.array2string(clc_tensor, precision=0, separator=', '),
        #                                        np.array2string(exp_tensor, precision=0, separator=', '),
        #                                        similarity))

        return similarity[0][0]

    def __captain_write_onlycross_PSM(self, all_match_result, write_result_pkl_path, write_result_path,
                                      separate_FDR=1, clean_decoy=0):

        functionPickle = CFunctionPickle()
        functionPickle.write_data_to_pkl(all_match_result, write_result_pkl_path)

        index_all_match_result = 0
        with open(write_result_path, 'w', encoding='utf-8') as f:
            f.write('Order\tTitle\tCharge\tPrecursor_Mass\tPeptide\tPeptide_Type'
                    '\tLinker\tPeptide_Mass\tModifications\tScore\tPrecursor_Mass_Error(Da)'
                    '\tPrecursor_Mass_Error(ppm)\tProteins\tProtein_Type\tAlpha_Matched'
                    '\tBeta_Matched\tMatch Ion\tLink type\tInter FDR\tIntra FDR\tFDR\tRerank Score\tMobility\n')
            for one_match_result in all_match_result:
                if clean_decoy == 1 and one_match_result.protein_type_math != 1:
                    pass
                else:
                    f.write(str(index_all_match_result + 1))
                    f.write('\t')
                    f.write(one_match_result.spectrum_title)  # Title
                    f.write('\t')
                    f.write(one_match_result.spectrum_charge)  # Charge
                    f.write('\t')
                    f.write(one_match_result.spectrum_mass)  # Precursor_mass
                    f.write('\t')
                    f.write(one_match_result.peptide_sq)
                    f.write('\t')
                    f.write(one_match_result.peptide_type)  # Peptide_type
                    f.write('\t')
                    f.write(one_match_result.linker)  # Linker
                    f.write('\t')
                    f.write(one_match_result.peptide_mass)  # Peptide_Mass
                    f.write('\t')
                    f.write(one_match_result.peptide_modification)  # Modification
                    f.write('\t')
                    f.write('%.4f' % one_match_result.peptide_score)  # Score
                    f.write('\t')
                    f.write(one_match_result.precursor_mass_error_Da)  # Precursor_Mass_Error(Da)
                    f.write('\t')
                    f.write(one_match_result.precursor_mass_error_ppm)  # Precursor_Mass_Error(ppm)
                    f.write('\t')
                    f.write(one_match_result.peptide_protein)
                    f.write('\t')
                    f.write(one_match_result.protein_type)  # Protein_Type
                    f.write('\t')
                    f.write("1")
                    # f.write("%.4f"%one_match_result.feature.match_alpha_score)  # Alpha_Matched
                    f.write('\t')
                    f.write("1")
                    # f.write("%.4f"%one_match_result.feature.match_beta_score)  # Beta_Matched
                    f.write('\t')
                    f.write("1")
                    # f.write("%.4f"%one_match_result.feature.match_ion_score)
                    f.write('\t')
                    if one_match_result.inter_protein:
                        f.write('Inter-protein')
                    else:
                        f.write('Intra-protein')
                    f.write('\t')
                    if one_match_result.inter_protein:
                        f.write('%.4f' % one_match_result.inter_FDR)
                    else:
                        f.write('')
                    f.write('\t')
                    if one_match_result.inter_protein:
                        f.write('')
                    else:
                        f.write('%.4f' % one_match_result.intra_FDR)
                    f.write('\t')
                    f.write('%.4f' % one_match_result.FDR)
                    f.write('\t')
                    f.write('%.4f' % one_match_result.rerank_score)
                    f.write('\t')
                    if one_match_result.mobility is None:
                        f.write('none')
                    else:
                        f.write('%.4f'%one_match_result.mobility)
                    f.write('\n')
                    index_all_match_result += 1
