# -*- mode: python ; coding: utf-8 -*-..
import os.path

from MSFunction import CFunctionPickle, CFunctionFDR
from MSSysterm import FILE_NAME_PSM_RESULT, VALUE_MAX_SCAN

import operator

class CFunctionOutOnePeptidePSM:

    FDR_threshold = 1
    is_rerank = 0
    clean_decoy = 0

    def __init__(self, inputDP):

        self.dp = inputDP

    def update_class_variable(self, FDR_threshold, out_type, clean_decoy):

        self.FDR_threshold = FDR_threshold
        self.is_rerank = out_type
        self.clean_decoy = clean_decoy

    def outputSinglePeptidePSM(self, all_match_result, write_pkl_path, write_path, write_separate=True):

        write_single_result_path = write_path + FILE_NAME_PSM_RESULT[1][1] + '.txt'
        write_single_result_pkl_path = write_pkl_path + FILE_NAME_PSM_RESULT[1][1] + '.pkl'

        if write_separate:

            write_mono_result_path = write_path + "result-mono_link" + '.txt'
            write_mono_result_pkl_path = write_pkl_path + "result-mono_link" + '.pkl'

            list_mono = []
            list_single = []
            for one_match_result in all_match_result:
                if one_match_result.momo_link:
                    list_mono.append(one_match_result)
                else:
                    list_single.append(one_match_result)
            # is_rerank:0是PSM，1是rerank后的PSM
            if self.is_rerank:
                cmpfun = operator.attrgetter('rerank_score')
                list_single.sort(key=cmpfun, reverse=True)
                cmpfun = operator.attrgetter('rerank_score')
                list_mono.sort(key=cmpfun, reverse=True)

            else:
                cmpfun = operator.attrgetter('peptide_score')
                list_single.sort(key=cmpfun, reverse=True)
                cmpfun = operator.attrgetter('peptide_score')
                list_mono.sort(key=cmpfun, reverse=True)
            #  获取FDR
            functionFDR = CFunctionFDR(self.dp)
            functionFDR.calculate_single_result_FDR(list_single)
            functionFDR.calculate_single_result_FDR(list_mono)
            # 写检索结果
            self.__captain_write_singlepeptide_PSM(list_single, write_single_result_pkl_path, write_single_result_path)
            self.__captain_write_monopeptide_PSM(list_mono, write_mono_result_pkl_path, write_mono_result_path)

        else:

            # is_rerank:0是PSM，1是rerank后的PSM
            if self.is_rerank:
                cmpfun = operator.attrgetter('rerank_score')
                all_match_result.sort(key=cmpfun, reverse=True)

            else:
                cmpfun = operator.attrgetter('peptide_score')
                all_match_result.sort(key=cmpfun, reverse=True)
            #  获取FDR
            functionFDR = CFunctionFDR(self.dp)
            functionFDR.calculate_single_result_FDR(all_match_result)
            # 写检索结果
            self.__captain_write_mono_single_PSM(all_match_result, write_single_result_pkl_path, write_single_result_path)

    def __captain_write_mono_single_PSM(self, all_match_result, write_result_pkl_path, write_result_path):

        functionPickle = CFunctionPickle()
        functionPickle.write_data_to_pkl(all_match_result, write_result_pkl_path)

        index_all_match_result = 0
        with open(write_result_path, 'w') as f:
            f.write('Order\tTitle\tCharge\tPrecursor_Mass\tPeptide\tPeptide_Type\tLinker'
                    '\tPeptide_Mass\tModifications\tScore\tPrecursor_Mass_Error(Da)'
                    '\tPrecursor_Mass_Error(ppm)\tProteins\tProtein_Type\tPeptide_Matched'
                    '\tFDR\tRerank Score\tMobility\n')
            for one_match_result in all_match_result:
                if one_match_result.FDR > self.FDR_threshold:
                    break
                if self.clean_decoy == 1 and one_match_result.protein_type_math != 1:
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
                    if one_match_result.momo_link:
                        f.write(one_match_result.peptide_mono)
                    else:
                        f.write('Common')
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
                    f.write("%.4f"%one_match_result.feature.match_score)  # Matched
                    f.write('\t')
                    f.write('%.4f' % one_match_result.FDR)
                    f.write('\t')
                    f.write('%.4f' % one_match_result.rerank_score)
                    f.write('\t')
                    if one_match_result.mobility is None:
                        f.write('none')
                    else:
                        f.write('%.4f' % one_match_result.mobility)
                    f.write('\n')
                    index_all_match_result += 1

    def __captain_write_singlepeptide_PSM(self, all_match_result, write_result_pkl_path, write_result_path):

        functionPickle = CFunctionPickle()
        functionPickle.write_data_to_pkl(all_match_result, write_result_pkl_path)

        index_all_match_result = 0
        with open(write_result_path, 'w') as f:
            f.write('Order\tTitle\tCharge\tPrecursor_Mass\tPeptide\tPeptide_Type'
                    '\tPeptide_Mass\tModifications\tScore\tPrecursor_Mass_Error(Da)'
                    '\tPrecursor_Mass_Error(ppm)\tProteins\tProtein_Type\tPeptide_Matched'
                    '\tFDR\tRerank Score\tMobility\n')
            for one_match_result in all_match_result:
                if one_match_result.FDR > self.FDR_threshold:
                    break
                if self.clean_decoy == 1 and one_match_result.protein_type_math != 1:
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
                    f.write("%.4f"%one_match_result.feature.match_score)  # Matched
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

    def __captain_write_monopeptide_PSM(self, all_match_result, write_result_pkl_path, write_result_path):

        functionPickle = CFunctionPickle()
        functionPickle.write_data_to_pkl(all_match_result, write_result_pkl_path)

        index_all_match_result = 0
        with open(write_result_path, 'w') as f:
            f.write('Order\tTitle\tCharge\tPrecursor_Mass\tPeptide\tPeptide_Type\tLinker'
                    '\tPeptide_Mass\tModifications\tScore\tPrecursor_Mass_Error(Da)'
                    '\tPrecursor_Mass_Error(ppm)\tProteins\tProtein_Type\tPeptide_Matched'
                    '\tFDR\tRerank Score\tMobility\n')
            for one_match_result in all_match_result:
                if one_match_result.FDR > self.FDR_threshold:
                    break
                if self.clean_decoy == 1 and one_match_result.protein_type_math != 1:
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
                    f.write(one_match_result.peptide_mono)
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
                    f.write("%.4f"%one_match_result.feature.match_score)  # Matched
                    f.write('\t')
                    f.write('%.4f' % one_match_result.FDR)
                    f.write('\t')
                    f.write('%.4f' % one_match_result.rerank_score)
                    if one_match_result.mobility is None:
                        f.write('none')
                    else:
                        f.write('%.4f'%one_match_result.mobility)
                    f.write('\n')
                    index_all_match_result += 1

    def outputLooplinkPSM(self, all_match_result, write_pkl_path, write_path):

        write_loop_link_result_path = write_path + FILE_NAME_PSM_RESULT[1][2] + '.txt'
        write_loop_link_result_pkl_path = write_pkl_path + FILE_NAME_PSM_RESULT[1][2] + '.pkl'

        # is_rerank:0是PSM，1是rerank后的PSM
        if self.is_rerank:
            cmpfun = operator.attrgetter('rerank_score')
            all_match_result.sort(key=cmpfun, reverse=True)

        else:
            cmpfun = operator.attrgetter('peptide_score')
            all_match_result.sort(key=cmpfun, reverse=True)

        #  获取FDR
        functionFDR = CFunctionFDR(self.dp)
        functionFDR.calculate_single_result_FDR(all_match_result)
        # 写检索结果
        self.__captain_write_loop_link_PSM(all_match_result, write_loop_link_result_pkl_path, write_loop_link_result_path)

    def __captain_write_loop_link_PSM(self, all_match_result, write_result_pkl_path, write_result_path):

        functionPickle = CFunctionPickle()
        functionPickle.write_data_to_pkl(all_match_result, write_result_pkl_path)

        index_all_match_result = 0
        with open(write_result_path, 'w') as f:
            f.write('Order\tTitle\tCharge\tPrecursor_Mass\tPeptide\tPeptide_Type\tLinker'
                    '\tPeptide_Mass\tModifications\tScore\tPrecursor_Mass_Error(Da)'
                    '\tPrecursor_Mass_Error(ppm)\tProteins\tProtein_Type\tPeptide_Matched'
                    '\tFDR\tRerank Score\tMobility\n')
            for one_match_result in all_match_result:
                if one_match_result.FDR > self.FDR_threshold:
                    break
                if self.clean_decoy == 1 and one_match_result.protein_type_math != 1:
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
                    f.write(one_match_result.peptide_loop)
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
                    f.write("%.4f"%one_match_result.feature.match_score)
                    f.write('\t')
                    f.write('%.4f' % one_match_result.FDR)
                    f.write('\t')
                    f.write('%.4f' % one_match_result.rerank_score)
                    if one_match_result.mobility is None:
                        f.write('none')
                    else:
                        f.write('%.4f'%one_match_result.mobility)
                    f.write('\n')
                    index_all_match_result += 1


class CFunctionOutTwoPeptidePSM:

    FDR_threshold = 1
    is_rerank = 0
    clean_decoy = 0

    def __init__(self, inputDP):
        self.dp = inputDP

    def update_class_variable(self, FDR_threshold, is_rerank, clean_decoy):

        self.FDR_threshold = FDR_threshold
        self.is_rerank = is_rerank
        self.clean_decoy = clean_decoy

    def outputOnlyCrossPSM(self, all_match_result, write_result_pkl_path, write_result_path, separate_FDR=True):
        # outtype:0是PSM，1是rerank后的PSM
        if self.is_rerank == 0:
            cmpfun = operator.attrgetter('peptide_score')
            all_match_result.sort(key=cmpfun, reverse=True)
        else:
            cmpfun = operator.attrgetter('rerank_score')
            all_match_result.sort(key=cmpfun, reverse=True)
        #  获取FDR
        functionFDR = CFunctionFDR(self.dp)
        functionFDR.calculate_cross_result_FDR(all_match_result)
        # 写检索结果
        self.__captain_write_onlycross_PSM(all_match_result, write_result_pkl_path, write_result_path, separate_FDR)

    def __captain_write_onlycross_PSM(self, all_match_result, write_result_pkl_path, write_result_path, separate_FDR=1):

        functionPickle = CFunctionPickle()
        functionPickle.write_data_to_pkl(all_match_result, write_result_pkl_path)

        index_all_match_result = 0
        with open(write_result_path, 'w') as f:
            f.write('Order\tTitle\tCharge\tPrecursor_Mass\tPeptide\tPeptide_Type'
                    '\tLinker\tPeptide_Mass\tModifications\tScore\tPrecursor_Mass_Error(Da)'
                    '\tPrecursor_Mass_Error(ppm)\tProteins\tProtein_Type\tAlpha_Matched'
                    '\tBeta_Matched\tMatch Ion\tLink type\tInter FDR\tIntra FDR\tFDR\tRerank Score\tMobility\n')
            for one_match_result in all_match_result:
                if separate_FDR:
                    if one_match_result.inter_FDR > self.FDR_threshold and one_match_result.intra_FDR > self.FDR_threshold:
                        break
                    if one_match_result.inter_FDR > self.FDR_threshold:
                        continue
                    if one_match_result.intra_FDR > self.FDR_threshold:
                        continue
                else:
                    if one_match_result.FDR > self.FDR_threshold:
                        break
                if self.clean_decoy == 1 and one_match_result.protein_type_math != 1:
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
                    f.write("%.4f"%one_match_result.feature.match_alpha_score)  # Alpha_Matched
                    f.write('\t')
                    f.write("%.4f"%one_match_result.feature.match_beta_score)  # Beta_Matched
                    f.write('\t')
                    f.write("%.4f"%one_match_result.feature.match_ion_score)
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

class CFunctionReadReasult():

    def __init__(self, inputDP):
        self.dp = inputDP

    def get_id_spectrum(self, FDR_threshold, record_spectrum_name, record_spectrum_scan):

        psm_result = self.read_psm_result()

        for one_result in psm_result:
            if one_result.FDR < FDR_threshold and one_result.protein_type_math == 1:
                title_name = one_result.spectrum_title.split('.')[0]
                scan = int(one_result.spectrum_title.split('.')[2])
                if scan == 16451:
                    print('aa')
                if title_name in record_spectrum_name:
                    record_spectrum_scan[record_spectrum_name.index(title_name)][scan].append(one_result.spectrum_title)
                else:
                    record_spectrum_name.append(title_name)
                    record_spectrum_scan.append([[] for i in range(VALUE_MAX_SCAN)])
                    record_spectrum_scan[-1][scan].append(one_result.spectrum_title)
            else:
                break
        del psm_result

    def read_psm_result(self):

        if self.dp.myCFG.D14_TYPE_PEPTIDE < len(FILE_NAME_PSM_RESULT):

            if self.dp.myCFG.D15_TYPE_LINK < len(FILE_NAME_PSM_RESULT[self.dp.myCFG.D14_TYPE_PEPTIDE]):

                functionPickle = CFunctionPickle()
                path_psm_result = self.dp.myCFG.A5_PATH_RESULT_EXPORT + '\\tmp\\res-target-' + FILE_NAME_PSM_RESULT[self.dp.myCFG.D14_TYPE_PEPTIDE][self.dp.myCFG.D15_TYPE_LINK] + '.pkl'
                psm_result = functionPickle.load_pkl_to_data(path_psm_result)

            else:

                psm_result = []

        else:

            psm_result = []

        return psm_result