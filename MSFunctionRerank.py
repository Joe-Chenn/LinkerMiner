# -*- mode: python ; coding: utf-8 -*-..
from MSFunction import CFunctionFDR, CFunctionPickle
from MSFunctionResult import CFunctionReadReasult, CFunctionOutOnePeptidePSM, CFunctionOutTwoPeptidePSM
from MSSysterm import FILE_NAME_PSM_RESULT, FILE_NAME_RERANK
from MSSysterm import PATH_LIBSVM, SCALE_EXE, TRAIN_EXE, PREDICT_EXE
from MSLogging import logGetError

import os
import operator
import matplotlib.pyplot as plt


class CFunctionRerank():

    def __init__(self):

        pass

    def delete_old_file(self, test_scale_file, train_scale_file, model_file, out_file):

        # 删掉旧的文件
        if os.path.exists(test_scale_file):
            os.remove(test_scale_file)
        if os.path.exists(train_scale_file):
            os.remove(train_scale_file)
        if os.path.exists(model_file):
            os.remove(model_file)
        if os.path.exists(out_file):
            os.remove(out_file)

    def get_data_scale(self, test_file, test_scale_file):
        # 执行的libsvm代码
        cmd_code = 'cd ' + os.getcwd() + '\\' + PATH_LIBSVM + ' && '
        test_scale_code = cmd_code + SCALE_EXE + ' ' + test_file + ' > ' + test_scale_file
        # 归一化、
        os.system(test_scale_code)

    def get_train_dataset_separate_FDR(self, search_result, test_scale_file, train_scale_file, target_FDR_threshold=1.0, decoy_FDR_threshold=1.0):
        # 输入正例和负例的数目直接从所有的数据中拿
        with open(test_scale_file, 'r') as f:
            test_scale_data = f.readlines()
        target_num = 0
        decoy_num = 0
        # 这是挑用来训练的结果的特征写入训练文件中
        with open(train_scale_file, 'w') as f:
            for result_index, one_result in enumerate(search_result):
                if one_result.protein_type_math is None:
                    break

                if one_result.inter_protein:
                    # inter-protein
                    if one_result.protein_type_math == 1:
                        # 正例
                        if one_result.inter_FDR <= target_FDR_threshold:
                            f.write(test_scale_data[result_index])
                            target_num += 1
                    elif one_result.protein_type_math == 0:
                        # 负例
                        if one_result.inter_FDR <= decoy_FDR_threshold:
                            f.write(test_scale_data[result_index])
                            decoy_num += 1
                else:
                    # intra-protein
                    if one_result.protein_type_math == 1:
                        # 正例
                        if one_result.intra_FDR <= target_FDR_threshold:
                            f.write(test_scale_data[result_index])
                            target_num += 1
                    elif one_result.protein_type_math == 0:
                        # 负例
                        if one_result.intra_FDR <= decoy_FDR_threshold:
                            f.write(test_scale_data[result_index])
                            decoy_num += 1
        print("[Info] Target num = " + str(target_num) + " Decoy_num = " +str(decoy_num))
        return target_num, decoy_num

    def get_train_dataset_global_FDR(self, search_result, test_scale_file, train_scale_file, target_FDR_threshold=1.0, decoy_FDR_threshold=1.0):
        # 输入正例和负例的数目直接从所有的数据中拿
        with open(test_scale_file, 'r') as f:
            test_scale_data = f.readlines()
        target_num = 0
        decoy_num = 0
        # 这是挑用来训练的结果的特征写入训练文件中
        with open(train_scale_file, 'w') as f:
            for result_index, one_result in enumerate(search_result):
                if one_result.protein_type_math is None:
                    break
                if one_result.protein_type_math == 1:
                    # 正例
                    if one_result.FDR <= target_FDR_threshold:
                        f.write(test_scale_data[result_index])
                        target_num += 1
                else:
                    # 负例
                    if one_result.FDR <= decoy_FDR_threshold:
                        f.write(test_scale_data[result_index])
                        decoy_num += 1
        print("[Info] Target num : %d"%target_num)
        print("[Info] Decoy num : %d" % decoy_num)

        return target_num, decoy_num

    def run_libsvm(self, train_scale_file, test_scale_file, model_file, out_file, C):
        # 执行的libsvm代码
        cmd_code = 'cd ' + os.getcwd() + '\\' + PATH_LIBSVM + ' && '
        train_code = cmd_code + TRAIN_EXE + ' -b 1 -h 0  -c ' + str(C) + ' ' + train_scale_file + ' ' + model_file
        predict_code = cmd_code + PREDICT_EXE + ' -b 1 ' + test_scale_file + ' ' + model_file + ' ' + out_file

        #训练、预测
        os.system(train_code)
        os.system(predict_code)

    def update_rerank_score(self, out_file, search_result):

        # 输入数据进行SVM训练后的打分结果，获取分数
        score_list = []
        pos_flag = 1
        with open(out_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line == "":
                    continue
                if line.startswith("labels"):
                    tmp = line.split(' ')
                    if tmp[2] == '1':
                        pos_flag = 2
                    continue
                tmp = line.split(' ')
                score_list.append(float(tmp[pos_flag]))

        # 分数从高到底获取
        for sort_index in range(len(score_list)):
            search_result[sort_index].rerank_score = score_list[sort_index]  # list_sort_index[sort_index]是分高到低的在原结果里的索引号
            search_result[sort_index].feature.rerank_score = score_list[sort_index]
        # 按rerank分数排好序的结果列表
        cmpfun = operator.attrgetter('rerank_score')
        search_result.sort(key=cmpfun, reverse=True)


class CFunctionRerankOnePeptide():

    target_FDR_threshold = 1
    decoy_FDR_threshold = 1

    def __init__(self, inputDP):

        self.dp = inputDP

    def updata_class_variable(self, target_FDR_threshold, decoy_FDR_threshold):

        self.target_FDR_threshold = target_FDR_threshold
        self.decoy_FDR_threshold = decoy_FDR_threshold

    def rerank_singlepeptide(self):
        # 获取匹配的详细信息，第一次是从检索结果来
        funcitonPickle = CFunctionPickle()
        psm_result_file = self.dp.myCFG.A5_PATH_RESULT_EXPORT + "\\tmp\\" + FILE_NAME_PSM_RESULT[1][1] + '.pkl'
        search_result = funcitonPickle.load_pkl_to_data(psm_result_file)

        functionOutOnePeptidePSM = CFunctionOutOnePeptidePSM(self.dp)
        # target_result_pkl_path = self.dp.myCFG.A5_PATH_RESULT_EXPORT + "\\tmp\\target-"
        # target_result_path = self.dp.myCFG.A5_PATH_RESULT_EXPORT + "target-"
        # functionOutOnePeptidePSM.update_class_variable(1, 0, 1)
        # functionOutOnePeptidePSM.outputSinglePeptidePSM(search_result, target_result_pkl_path, target_result_path)

        rerank_folder = self.dp.myCFG.A5_PATH_RESULT_EXPORT + "\\tmp\\" + '\\' + FILE_NAME_RERANK[1][1]
        if not os.path.exists(rerank_folder):
            os.makedirs(rerank_folder)

        functionRerank = CFunctionRerank()
        functionFDR = CFunctionFDR(self.dp)

        for i in range(self.dp.myCFG.single_rerank_time):  # 进行5次
            # 上次写好的结果pkl文件
            # 训练集的特征值文件，后面还要归一化处理
            # 测试集的特征值文件，后面还要归一化处理
            # 训练的模型文件
            # 训练好的模型预测的结果
            # 新的写好的训练好的预测结果的pkl文件

            test_file = rerank_folder + "\\" + "test_" + FILE_NAME_RERANK[1][1] + '-' + str(i) + '.txt'
            # 特征归一化后的两个数据集的特征值文件
            train_scale_file = rerank_folder + "\\" + "train_scale_" + FILE_NAME_RERANK[1][1] + '-' + str(i) + '.txt'
            test_scale_file = rerank_folder + "\\" + "test_scale" + FILE_NAME_RERANK[1][1] + '-' + str(i) + '.txt'
            model_file = rerank_folder + "\\" + "model_" + FILE_NAME_RERANK[1][1] + '-' + str(i) + '.model'
            out_file = rerank_folder + "\\" + "out_" + FILE_NAME_RERANK[1][1] + '-' + str(i) + '.txt'

            functionRerank.delete_old_file(test_scale_file, train_scale_file, model_file, out_file)
            # 这里是为了看特征分布
            self.__captain_show_singlepeptide_feature_distribution(search_result, rerank_folder, 1, 1, i)

            # 这里是取正例和负例的结果作为训练内容
            self.__captain_get_singlepeptide_test_dataset(search_result, test_file)
            functionRerank.get_data_scale(test_file, test_scale_file)

            target_num, decoy_num = functionRerank.get_train_dataset_global_FDR(search_result, test_scale_file, train_scale_file, self.target_FDR_threshold, self.decoy_FDR_threshold)
            if target_num == 0 or decoy_num == 0:
                break
            else:
                if target_num > decoy_num:
                    C = target_num / decoy_num
                else:
                    C = decoy_num / target_num
            functionRerank.run_libsvm(train_scale_file, test_scale_file, model_file, out_file, C)

            functionRerank.update_rerank_score(out_file, search_result)

            functionFDR.calculate_single_result_FDR(search_result)

        rerank_result_path = self.dp.myCFG.A5_PATH_RESULT_EXPORT + "\\tmp\\res-"
        rerank_result_pkl_path = self.dp.myCFG.A5_PATH_RESULT_EXPORT + "\\tmp\\res-"
        functionOutOnePeptidePSM.update_class_variable(1, 1, 0)
        functionOutOnePeptidePSM.outputSinglePeptidePSM(search_result, rerank_result_pkl_path, rerank_result_path)

        rerank_target_result_path = self.dp.myCFG.A5_PATH_RESULT_EXPORT + "\\res-target-"
        rerank_target_result_pkl_path = self.dp.myCFG.A5_PATH_RESULT_EXPORT + "\\tmp\\res-target-"
        functionOutOnePeptidePSM.update_class_variable(self.dp.myCFG.E1_FDR_PSM, 1, 1)
        functionOutOnePeptidePSM.outputSinglePeptidePSM(search_result, rerank_target_result_pkl_path, rerank_target_result_path)

    def __captain_show_singlepeptide_feature_distribution(self, search_result, rerank_folder, target_FDR_threshold=1.0, decoy_FDR_threshold=1.0, input_i=None):

        target_data_rerank_score = []
        target_data_match_score = []

        target_data_match_error_sum = []
        target_data_match_error_average = []
        target_data_match_error_var = []

        target_data_match_ion_num = []
        target_data_match_ion_intensity = []

        target_data_match_ion_num_percent = []
        target_data_match_ion_intensity_percent = []

        target_data_match_spe_ion_percent = []
        target_data_match_spe_inten_percent = []

        target_data_continue_peptide_score = []
        target_data_precursor_bias = []

        target_data_delta_score = []

        target_data_pParseNum = []

        # =======================================================
        decoy_data_rerank_score = []
        decoy_data_match_score = []

        decoy_data_match_error_sum = []
        decoy_data_match_error_average = []
        decoy_data_match_error_var = []

        decoy_data_match_ion_num = []
        decoy_data_match_ion_intensity = []

        decoy_data_match_ion_num_percent = []
        decoy_data_match_ion_intensity_percent = []

        decoy_data_match_spe_ion_percent = []
        decoy_data_match_spe_inten_percent = []

        decoy_data_continue_peptide_score = []
        decoy_data_precursor_bias = []

        decoy_data_delta_score = []

        decoy_data_pParseNum = []

        for n, one_search_result in enumerate(search_result):
            '''
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
            '''
            if one_search_result.FDR > target_FDR_threshold and one_search_result.protein_type_math == 1:
                continue
            if one_search_result.protein_type_math == 1 and one_search_result.FDR < target_FDR_threshold:
                # 分数类的特征
                target_data_rerank_score.append(one_search_result.feature.rerank_score)
                target_data_match_score.append(one_search_result.feature.match_score)
                # 匹配碎片离子偏差类特征
                target_data_match_error_sum.append(one_search_result.feature.match_error_sum)
                target_data_match_error_average.append(one_search_result.feature.match_error_average)
                target_data_match_error_var.append(one_search_result.feature.match_error_var)

                ## 匹配上数目和总强度类特征
                target_data_match_ion_num.append(one_search_result.feature.match_ion_num)
                target_data_match_ion_intensity.append(one_search_result.feature.match_ion_intensity)
                # 匹配上离子的百分比特征
                target_data_match_ion_num_percent.append(one_search_result.feature.match_ion_num_percent)
                target_data_match_ion_intensity_percent.append(one_search_result.feature.match_ion_intensity_percent)
                # 匹配对于谱图的百分比特征
                target_data_match_spe_ion_percent.append(one_search_result.feature.match_spe_ion_percent)
                target_data_match_spe_inten_percent.append(one_search_result.feature.match_spe_inten_percent)
                # 连续性的特征
                target_data_continue_peptide_score.append(one_search_result.feature.continue_peptide_score)
                # 其他和分数以及匹配情况不太相关的特征
                target_data_precursor_bias.append(one_search_result.feature.precursor_bias)
                target_data_delta_score.append(one_search_result.feature.delta_score)
                target_data_pParseNum.append(one_search_result.feature.pParseNum)

            elif (one_search_result.protein_type_math == 0 or one_search_result.protein_type_math == -1) and one_search_result.FDR < decoy_FDR_threshold:
                # 分数类的特征
                decoy_data_rerank_score.append(one_search_result.feature.rerank_score)
                target_data_match_score.append(one_search_result.feature.match_score)
                # 匹配碎片离子偏差类特征
                decoy_data_match_error_sum.append(one_search_result.feature.match_error_sum)
                decoy_data_match_error_average.append(one_search_result.feature.match_error_average)
                decoy_data_match_error_var.append(one_search_result.feature.match_error_var)

                ## 匹配上数目和总强度类特征
                decoy_data_match_ion_num.append(one_search_result.feature.match_ion_num)
                decoy_data_match_ion_intensity.append(one_search_result.feature.match_ion_intensity)
                # 匹配上离子的百分比特征
                decoy_data_match_ion_num_percent.append(one_search_result.feature.match_ion_num_percent)
                decoy_data_match_ion_intensity_percent.append(one_search_result.feature.match_ion_intensity_percent)
                # 匹配对于谱图的百分比特征
                decoy_data_match_spe_ion_percent.append(one_search_result.feature.match_spe_ion_percent)
                decoy_data_match_spe_inten_percent.append(one_search_result.feature.match_spe_inten_percent)
                # 连续性的特征
                decoy_data_continue_peptide_score.append(one_search_result.feature.continue_peptide_score)
                # 其他和分数以及匹配情况不太相关的特征
                decoy_data_precursor_bias.append(one_search_result.feature.precursor_bias)
                decoy_data_delta_score.append(one_search_result.feature.delta_score)
                decoy_data_pParseNum.append(one_search_result.feature.pParseNum)
            else:
                pass

        f, axs = plt.subplots(5, 3, figsize=(60, 20))

        axs[0][0].hist([target_data_rerank_score, decoy_data_rerank_score], label=['Target', 'Decoy'], bins=100,
                       stacked=True)
        axs[0][0].title.set_text('rerank_score')

        axs[0][1].hist([target_data_match_score, decoy_data_match_score], label=['Target', 'Decoy'],
                       bins=100, stacked=True)
        axs[0][1].title.set_text('match_score')

        # ==============================
        axs[0][2].hist([target_data_match_error_sum, decoy_data_match_error_sum], label=['Target', 'Decoy'],
                       bins=100, stacked=True)
        axs[0][2].title.set_text('match_error_sum')

        axs[1][0].hist([target_data_match_error_average, decoy_data_match_error_average],
                       label=['Target', 'Decoy'], bins=100, stacked=True,)
        axs[1][0].title.set_text('match_error_average')

        axs[1][1].hist([target_data_match_error_var, decoy_data_match_error_var],
                       label=['Target', 'Decoy'], bins=100, stacked=True)
        axs[1][1].title.set_text('match_error_var')

        axs[1][2].hist([target_data_match_ion_num, decoy_data_match_ion_num],
                       label=['Target', 'Decoy'], bins=100, stacked=True)
        axs[1][2].title.set_text('match_ion_num')
        # ==============================
        axs[2][0].hist([target_data_match_ion_intensity, decoy_data_match_ion_intensity],
                       label=['Target', 'Decoy'], bins=100, stacked=True)
        axs[2][0].title.set_text('match_ion_intensity')

        axs[2][1].hist([target_data_match_ion_num_percent, decoy_data_match_ion_num_percent], label=['Target', 'Decoy'],
                       bins=100, stacked=True)
        axs[2][1].title.set_text('match_ion_num_percent')

        axs[2][2].hist([target_data_match_ion_intensity_percent, decoy_data_match_ion_intensity_percent], label=['Target', 'Decoy'], bins=100,
                       stacked=True)
        axs[2][2].title.set_text('match_ion_intensity_percent')

        # ==============================
        axs[3][0].hist([target_data_match_spe_ion_percent, decoy_data_match_spe_ion_percent], label=['Target', 'Decoy'],
                       bins=100, stacked=True)
        axs[3][0].title.set_text('match_spe_ion_percent')

        axs[3][1].hist([target_data_match_spe_inten_percent, decoy_data_match_spe_inten_percent], label=['Target', 'Decoy'],
                       bins=100, stacked=True)
        axs[3][1].title.set_text('match_spe_inten_percent')

        axs[3][2].hist([target_data_continue_peptide_score, decoy_data_continue_peptide_score],
                       label=['Target', 'Decoy'], bins=100, stacked=True)
        axs[3][2].title.set_text('data_delta_score')

        axs[4][0].hist([target_data_precursor_bias, decoy_data_precursor_bias], label=['Target', 'Decoy'],
                       bins=100, stacked=True)
        axs[4][0].title.set_text('precursor_bias')

        axs[4][1].hist([target_data_delta_score, decoy_data_delta_score], label=['Target', 'Decoy'],
                       bins=100, stacked=True)
        axs[4][1].title.set_text('delta_score')

        axs[4][2].hist([target_data_pParseNum, decoy_data_pParseNum], label=['Target', 'Decoy'],
                       bins=100, stacked=True)
        axs[4][2].title.set_text('pParseNum')

        axs[0][1].legend()
        axs[0][1].legend()
        axs[0][1].legend()
        axs[1][0].legend()
        axs[1][1].legend()
        axs[1][2].legend()
        axs[2][0].legend()
        axs[2][1].legend()
        axs[2][2].legend()
        axs[3][0].legend()
        axs[3][1].legend()
        axs[3][2].legend()
        axs[4][0].legend()
        axs[4][1].legend()
        axs[4][2].legend()

        if input_i is None:
            plt.savefig(rerank_folder + "\\" + 'rerank_figure.png', format='png')
        else:
            plt.savefig(rerank_folder + "\\" + 'rerank_figure-' + str(input_i) + '.png', format='png')

    def __captain_get_singlepeptide_test_dataset(self, search_result, test_file):
        # 这是把所有数据都归一化写入到文件中
        with open(test_file, 'w') as f:
            for one_search_result in search_result:
                if len(one_search_result.protein_type) == 0:
                    break
                if one_search_result.protein_type_math == 1:
                    flag = '1'
                else:
                    flag = '-1'
                # write_str += (str(flag) +
                #                 " 1:" + str(one_search_result.feature.rerank_score) +
                #                 " 2:" + str(one_search_result.feature.match_score) +
                #
                #                 " 3:" + str(one_search_result.feature.match_error_sum) +
                #                 " 4:" + str(one_search_result.feature.match_error_average) +
                #                 " 5:" + str(one_search_result.feature.match_error_var) +
                #
                #                 " 6:" + str(one_search_result.feature.match_ion_num) +
                #                 " 7:" + str(one_search_result.feature.match_ion_intensity) +
                #
                #                 " 8:" + str(one_search_result.feature.match_ion_num_percent) +
                #                 " 9:" + str(one_search_result.feature.match_ion_intensity_percent) +
                #
                #                 " 10:" + str(one_search_result.feature.match_spe_ion_percent) +
                #                 " 11:" + str(one_search_result.feature.match_spe_inten_percent) +
                #
                #                 " 12:" + str(one_search_result.feature.continue_peptide_score) +
                #                 " 13:" + str(one_search_result.feature.precursor_bias) +
                #
                #                 " 14:" + str(one_search_result.feature.delta_score) +
                #                 " 15:" + str(one_search_result.feature.pParseNum) +
                #
                #               "\n")
                write_str = (str(flag) +
                                 " 1:" + str(one_search_result.feature.rerank_score) +
                                 " 2:" + str(one_search_result.feature.match_score) +

                                 # " 3:" + str(one_search_result.feature.match_error_sum) +
                                 # " 4:" + str(one_search_result.feature.match_error_average) +

                                 " 5:" + str(one_search_result.feature.match_ion_num) +

                                 " 6:" + str(one_search_result.feature.continue_peptide_score) +
                                 " 7:" + str(one_search_result.feature.precursor_bias) +

                                 " 8:" + str(one_search_result.feature.delta_score) +
                              "\n")
                f.write(write_str)

    def rerank_looplink(self):
        # 获取匹配的详细信息，第一次是从检索结果来
        funcitonPickle = CFunctionPickle()
        psm_result_file = self.dp.myCFG.A5_PATH_RESULT_EXPORT + "\\tmp\\" + FILE_NAME_PSM_RESULT[1][2] + '.pkl'
        search_result = funcitonPickle.load_pkl_to_data(psm_result_file)

        functionOutOnePeptidePSM = CFunctionOutOnePeptidePSM(self.dp)
        # target_result_path = self.dp.myCFG.A5_PATH_RESULT_EXPORT + "target-"
        # target_result_pkl_path = self.dp.myCFG.A5_PATH_RESULT_EXPORT + "\\tmp\\target-"
        # functionOutOnePeptidePSM.update_class_variable(1, 0, 1)
        # functionOutOnePeptidePSM.outputLooplinkPSM(search_result, target_result_pkl_path, target_result_path)

        rerank_folder = self.dp.myCFG.A5_PATH_RESULT_EXPORT + '\\tmp\\' + FILE_NAME_RERANK[1][2]
        if not os.path.exists(rerank_folder):
            os.makedirs(rerank_folder)

        functionRerank = CFunctionRerank()
        functionFDR = CFunctionFDR(self.dp)
        for i in range(self.dp.myCFG.rerank_time):  # 进行5次
            # 上次写好的结果pkl文件
            # 训练集的特征值文件，后面还要归一化处理
            # 测试集的特征值文件，后面还要归一化处理
            # 训练的模型文件
            # 训练好的模型预测的结果
            # 新的写好的训练好的预测结果的pkl文件

            test_file = rerank_folder + "\\" + "test_" + FILE_NAME_RERANK[1][2] + '-' + str(i) + '.txt'
            # 特征归一化后的两个数据集的特征值文件
            train_scale_file = rerank_folder + "\\" + "train_scale_" + FILE_NAME_RERANK[1][2] + '-' + str(i) + '.txt'
            test_scale_file = rerank_folder + "\\" + "test_scale" + FILE_NAME_RERANK[1][2] + '-' + str(i) + '.txt'
            model_file = rerank_folder + "\\" + "model_" + FILE_NAME_RERANK[1][2] + '-' + str(i) + '.model'
            out_file = rerank_folder + "\\" + "out_" + FILE_NAME_RERANK[1][2] + '-' + str(i) + '.txt'

            functionRerank.delete_old_file(test_scale_file, train_scale_file, model_file, out_file)
            # 这里是为了看特征分布
            self.__captain_show_looplink_feature_distribution(search_result, rerank_folder, 1, 1, i)

            # 这里是取正例和负例的结果作为训练内容
            self.__captain_get_looplink_test_dataset(search_result, test_file)
            functionRerank.get_data_scale(test_file, test_scale_file)

            target_num, decoy_num = functionRerank.get_train_dataset_global_FDR(search_result, test_scale_file, train_scale_file, self.target_FDR_threshold, self.decoy_FDR_threshold)
            if target_num == 0 or decoy_num == 0:
                break
            else:
                C = target_num / decoy_num
            functionRerank.run_libsvm(train_scale_file, test_scale_file, model_file, out_file, C)

            functionRerank.update_rerank_score(out_file, search_result)

            functionFDR.calculate_single_result_FDR(search_result)

        rerank_result_path = self.dp.myCFG.A5_PATH_RESULT_EXPORT + "\\tmp\\res-"
        rerank_result_pkl_path = self.dp.myCFG.A5_PATH_RESULT_EXPORT + "\\tmp\\res-"
        functionOutOnePeptidePSM.update_class_variable(self.dp.myCFG.E1_FDR_PSM, 1, 0)
        functionOutOnePeptidePSM.outputLooplinkPSM(search_result, rerank_result_pkl_path, rerank_result_path)

        rerank_target_result_path = self.dp.myCFG.A5_PATH_RESULT_EXPORT + "\\res-target-"
        rerank_target_result_pkl_path = self.dp.myCFG.A5_PATH_RESULT_EXPORT + "\\tmp\\res-target-"
        functionOutOnePeptidePSM.update_class_variable(self.dp.myCFG.E1_FDR_PSM, 1, 1)
        functionOutOnePeptidePSM.outputLooplinkPSM(search_result, rerank_target_result_pkl_path, rerank_target_result_path)

    def __captain_show_looplink_feature_distribution(self, search_result, rerank_folder, target_FDR_threshold=1.0, decoy_FDR_threshold=1.0, input_i=None):

        target_data_rerank_score = []
        target_data_match_score = []

        target_data_match_error_sum = []
        target_data_match_error_average = []
        target_data_match_error_var = []

        target_data_match_ion_num = []
        target_data_match_ion_intensity = []

        target_data_match_ion_num_percent = []
        target_data_match_ion_intensity_percent = []

        target_data_match_spe_ion_percent = []
        target_data_match_spe_inten_percent = []

        target_data_continue_peptide_score = []
        target_data_precursor_bias = []

        target_data_delta_score = []

        target_data_pParseNum = []

        # =======================================================
        decoy_data_rerank_score = []
        decoy_data_match_score = []

        decoy_data_match_error_sum = []
        decoy_data_match_error_average = []
        decoy_data_match_error_var = []

        decoy_data_match_ion_num = []
        decoy_data_match_ion_intensity = []

        decoy_data_match_ion_num_percent = []
        decoy_data_match_ion_intensity_percent = []

        decoy_data_match_spe_ion_percent = []
        decoy_data_match_spe_inten_percent = []

        decoy_data_continue_peptide_score = []
        decoy_data_precursor_bias = []

        decoy_data_delta_score = []

        decoy_data_pParseNum = []

        for n, one_search_result in enumerate(search_result):
            '''
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
            '''
            if one_search_result.FDR > target_FDR_threshold and one_search_result.protein_type_math == 1:
                continue
            if one_search_result.protein_type_math == 1 and one_search_result.FDR < target_FDR_threshold:
                # 分数类的特征
                target_data_rerank_score.append(one_search_result.feature.rerank_score)
                target_data_match_score.append(one_search_result.feature.match_score)
                # 匹配碎片离子偏差类特征
                target_data_match_error_sum.append(one_search_result.feature.match_error_sum)
                target_data_match_error_average.append(one_search_result.feature.match_error_average)
                target_data_match_error_var.append(one_search_result.feature.match_error_var)

                ## 匹配上数目和总强度类特征
                target_data_match_ion_num.append(one_search_result.feature.match_ion_num)
                target_data_match_ion_intensity.append(one_search_result.feature.match_ion_intensity)
                # 匹配上离子的百分比特征
                target_data_match_ion_num_percent.append(one_search_result.feature.match_ion_num_percent)
                target_data_match_ion_intensity_percent.append(one_search_result.feature.match_ion_intensity_percent)
                # 匹配对于谱图的百分比特征
                target_data_match_spe_ion_percent.append(one_search_result.feature.match_spe_ion_percent)
                target_data_match_spe_inten_percent.append(one_search_result.feature.match_spe_inten_percent)
                # 连续性的特征
                target_data_continue_peptide_score.append(one_search_result.feature.continue_peptide_score)
                # 其他和分数以及匹配情况不太相关的特征
                target_data_precursor_bias.append(one_search_result.feature.precursor_bias)
                target_data_delta_score.append(one_search_result.feature.delta_score)
                target_data_pParseNum.append(one_search_result.feature.pParseNum)

            elif (
                    one_search_result.protein_type_math == 0 or one_search_result.protein_type_math == -1) and one_search_result.FDR < decoy_FDR_threshold:
                # 分数类的特征
                decoy_data_rerank_score.append(one_search_result.feature.rerank_score)
                target_data_match_score.append(one_search_result.feature.match_score)
                # 匹配碎片离子偏差类特征
                decoy_data_match_error_sum.append(one_search_result.feature.match_error_sum)
                decoy_data_match_error_average.append(one_search_result.feature.match_error_average)
                decoy_data_match_error_var.append(one_search_result.feature.match_error_var)

                ## 匹配上数目和总强度类特征
                decoy_data_match_ion_num.append(one_search_result.feature.match_ion_num)
                decoy_data_match_ion_intensity.append(one_search_result.feature.match_ion_intensity)
                # 匹配上离子的百分比特征
                decoy_data_match_ion_num_percent.append(one_search_result.feature.match_ion_num_percent)
                decoy_data_match_ion_intensity_percent.append(one_search_result.feature.match_ion_intensity_percent)
                # 匹配对于谱图的百分比特征
                decoy_data_match_spe_ion_percent.append(one_search_result.feature.match_spe_ion_percent)
                decoy_data_match_spe_inten_percent.append(one_search_result.feature.match_spe_inten_percent)
                # 连续性的特征
                decoy_data_continue_peptide_score.append(one_search_result.feature.continue_peptide_score)
                # 其他和分数以及匹配情况不太相关的特征
                decoy_data_precursor_bias.append(one_search_result.feature.precursor_bias)
                decoy_data_delta_score.append(one_search_result.feature.delta_score)
                decoy_data_pParseNum.append(one_search_result.feature.pParseNum)
            else:
                pass

        f, axs = plt.subplots(5, 3, figsize=(20, 60))

        axs[0][0].hist([target_data_rerank_score, decoy_data_rerank_score], label=['Target', 'Decoy'], bins=100,
                       stacked=True)
        axs[0][0].title.set_text('rerank_score')

        axs[0][1].hist([target_data_match_score, decoy_data_match_score], label=['Target', 'Decoy'],
                       bins=100, stacked=True)
        axs[0][1].title.set_text('match_score')

        # ==============================
        axs[0][2].hist([target_data_match_error_sum, decoy_data_match_error_sum], label=['Target', 'Decoy'],
                       bins=100, stacked=True)
        axs[0][2].title.set_text('match_error_sum')

        axs[1][0].hist([target_data_match_error_average, decoy_data_match_error_average],
                       label=['Target', 'Decoy'], bins=100, stacked=True, )
        axs[1][0].title.set_text('match_error_average')

        axs[1][1].hist([target_data_match_error_var, decoy_data_match_error_var],
                       label=['Target', 'Decoy'], bins=100, stacked=True)
        axs[1][1].title.set_text('match_error_var')

        axs[1][2].hist([target_data_match_ion_num, decoy_data_match_ion_num],
                       label=['Target', 'Decoy'], bins=100, stacked=True)
        axs[1][2].title.set_text('match_ion_num')
        # ==============================
        axs[2][0].hist([target_data_match_ion_intensity, decoy_data_match_ion_intensity],
                       label=['Target', 'Decoy'], bins=100, stacked=True)
        axs[2][0].title.set_text('match_ion_intensity')

        axs[2][1].hist([target_data_match_ion_num_percent, decoy_data_match_ion_num_percent], label=['Target', 'Decoy'],
                       bins=100, stacked=True)
        axs[2][1].title.set_text('match_ion_num_percent')

        axs[2][2].hist([target_data_match_ion_intensity_percent, decoy_data_match_ion_intensity_percent],
                       label=['Target', 'Decoy'], bins=100,
                       stacked=True)
        axs[2][2].title.set_text('match_ion_intensity_percent')

        # ==============================
        axs[3][0].hist([target_data_match_spe_ion_percent, decoy_data_match_spe_ion_percent], label=['Target', 'Decoy'],
                       bins=100, stacked=True)
        axs[3][0].title.set_text('match_spe_ion_percent')

        axs[3][1].hist([target_data_match_spe_inten_percent, decoy_data_match_spe_inten_percent],
                       label=['Target', 'Decoy'],
                       bins=100, stacked=True)
        axs[3][1].title.set_text('match_spe_inten_percent')

        axs[3][2].hist([target_data_continue_peptide_score, decoy_data_continue_peptide_score],
                       label=['Target', 'Decoy'], bins=100, stacked=True)
        axs[3][2].title.set_text('data_delta_score')

        axs[4][0].hist([target_data_precursor_bias, decoy_data_precursor_bias], label=['Target', 'Decoy'],
                       bins=100, stacked=True)
        axs[4][0].title.set_text('precursor_bias')

        axs[4][1].hist([target_data_delta_score, decoy_data_delta_score], label=['Target', 'Decoy'],
                       bins=100, stacked=True)
        axs[4][1].title.set_text('delta_score')

        axs[4][2].hist([target_data_pParseNum, decoy_data_pParseNum], label=['Target', 'Decoy'],
                       bins=100, stacked=True)
        axs[4][2].title.set_text('pParseNum')

        axs[0][1].legend()
        axs[0][1].legend()
        axs[0][1].legend()
        axs[1][0].legend()
        axs[1][1].legend()
        axs[1][2].legend()
        axs[2][0].legend()
        axs[2][1].legend()
        axs[2][2].legend()
        axs[3][0].legend()
        axs[3][1].legend()
        axs[3][2].legend()
        axs[4][0].legend()
        axs[4][1].legend()
        axs[4][2].legend()

        if input_i is None:
            plt.savefig(rerank_folder + "\\" + 'rerank_figure.png', format='png')
        else:
            plt.savefig(rerank_folder + "\\" + 'rerank_figure-' + str(input_i) + '.png', format='png')

    def __captain_get_looplink_test_dataset(self, search_result, test_file):
        # 这是把所有数据都归一化写入到文件中
        with open(test_file, 'w') as f:
            for one_search_result in search_result:
                if len(one_search_result.protein_type) == 0:
                    break
                if one_search_result.protein_type_math == 1:
                    flag = '1'
                else:
                    flag = '-1'
                # write_str += (str(flag) +
                #               " 1:" + str(one_search_result.feature.rerank_score) +
                #               " 2:" + str(one_search_result.feature.match_score) +
                #
                #               " 3:" + str(one_search_result.feature.match_error_average) +
                #               " 4:" + str(one_search_result.feature.match_error_var) +
                #
                #               " 5:" + str(one_search_result.feature.match_ion_num_percent) +
                #               " 6:" + str(one_search_result.feature.match_ion_intensity_percent) +
                #
                #               " 7:" + str(one_search_result.feature.match_ion_num) +
                #               " 8:" + str(one_search_result.feature.match_ion_intensity) +
                #
                #               " 9:" + str(one_search_result.feature.continue_peptide_score) +
                #
                #               " 10:" + str(one_search_result.feature.spectrum_intensity) +
                #               " 11:" + str(one_search_result.feature.precursor_bias) +
                #               " 12:" + str(one_search_result.feature.delta_score) +
                #               " 13:" + str(one_search_result.feature.pParseNum) +
                #               "\n")
                write_str = (str(flag) +
                              " 1:" + str(one_search_result.feature.rerank_score) +
                              " 2:" + str(one_search_result.feature.match_score) +
                              # " 3:" + str(one_search_result.feature.match_ion_num) +
                              " 4:" + str(one_search_result.feature.match_ion_num_percent) +
                              " 5:" + str(one_search_result.feature.match_ion_intensity_percent) +
                              # " 6:" + str(one_search_result.feature.precursor_bias) +
                              " 7:" + str(one_search_result.feature.delta_score) +
                              " 8:" + str(one_search_result.feature.continue_peptide_score) +
                              " 9:" + str(one_search_result.feature.pParseNum) +
                              "\n")
                f.write(write_str)

class CFunctionRerankTwoPeptide():

    target_FDR_threshold = 1
    decoy_FDR_threshold = 1

    def __init__(self, inputDP):

        self.dp = inputDP

    def updata_class_variable(self, target_FDR_threshold, decoy_FDR_threshold):

        self.target_FDR_threshold = target_FDR_threshold
        self.decoy_FDR_threshold = decoy_FDR_threshold

    def rerank_onlycross(self):

        # 获取匹配的详细信息，第一次是从检索结果来
        psm_result_file = self.dp.myCFG.A5_PATH_RESULT_EXPORT + "\\tmp\\" + FILE_NAME_PSM_RESULT[2][1] + '.pkl'
        funcitonPickle = CFunctionPickle()
        search_result = funcitonPickle.load_pkl_to_data(psm_result_file)

        functionFDR = CFunctionFDR(self.dp)
        functionFDR.calculate_cross_result_FDR(search_result)

        functionOutTwoPeptidePSM = CFunctionOutTwoPeptidePSM(self.dp)
        # target_result_path = self.dp.myCFG.A5_PATH_RESULT_EXPORT + "target-" + FILE_NAME_PSM_RESULT[2][1] + '.txt'
        # target_pkl_path = self.dp.myCFG.A5_PATH_RESULT_EXPORT + '\\tmp\\' + "target-" + FILE_NAME_PSM_RESULT[2][1] + '.pkl'
        # functionOutTwoPeptidePSM.update_class_variable(self.dp.myCFG.E1_FDR_PSM, 0, 1)# FDR阈值 排序方式 清除反
        # functionOutTwoPeptidePSM.outputOnlyCrossPSM(search_result, target_pkl_path, target_result_path, self.dp.myCFG.separate_FDR)

        functionRerank = CFunctionRerank()
        rerank_folder = self.dp.myCFG.A5_PATH_RESULT_EXPORT + "\\tmp\\" + FILE_NAME_RERANK[2][1]
        if not os.path.exists(rerank_folder):
            os.makedirs(rerank_folder)

        for i in range(self.dp.myCFG.rerank_time):  # 进行5次
            # 上次写好的结果pkl文件
            # 训练集的特征值文件，后面还要归一化处理
            # 测试集的特征值文件，后面还要归一化处理
            # 训练的模型文件
            # 训练好的模型预测的结果
            # 新的写好的训练好的预测结果的pkl文件

            test_file = rerank_folder + "\\" + "test_" + FILE_NAME_RERANK[2][1] + '-' + str(i) + '.txt'
            # 特征归一化后的两个数据集的特征值文件
            train_scale_file = rerank_folder + "\\" + "train_scale_" + FILE_NAME_RERANK[2][1] + '-' + str(i) + '.txt'
            test_scale_file = rerank_folder + "\\" + "test_scale_" + FILE_NAME_RERANK[2][1] + '-' + str(i) + '.txt'
            model_file = rerank_folder + "\\" + "model_" + FILE_NAME_RERANK[2][1] + '-' + str(i) + '.model'
            out_file = rerank_folder + "\\" + "out_" + FILE_NAME_RERANK[2][1] + '-' + str(i) + '.txt'

            functionRerank.delete_old_file(test_scale_file, train_scale_file, model_file, out_file)
            # 这里是为了看特征分布
            self.__captain_show_crosslink_feature_distribution(search_result, rerank_folder, 1, 1, i)
            self.__captain_get_crosslink_test_dataset(search_result, test_file)

            functionRerank.get_data_scale(test_file, test_scale_file)
            if self.dp.myCFG.E2_FDR_SEPARATE:
                target_num, decoy_num = functionRerank.get_train_dataset_separate_FDR(search_result, test_scale_file, train_scale_file, self.target_FDR_threshold, self.decoy_FDR_threshold)
            else:
                target_num, decoy_num = functionRerank.get_train_dataset_global_FDR(search_result, test_scale_file, train_scale_file, self.target_FDR_threshold, self.decoy_FDR_threshold)

            if target_num == 0 or decoy_num == 0:
                break
            else:
                if target_num > decoy_num:
                    C = decoy_num / (target_num + decoy_num)
                else:
                    C = target_num / (target_num + decoy_num)

            # 这里是取正例和负例的结果作为训练内容
            functionRerank.run_libsvm(train_scale_file, test_scale_file, model_file, out_file, C)

            functionRerank.update_rerank_score(out_file, search_result)
            functionFDR.calculate_cross_result_FDR(search_result)

            # write_result_pkl_path = rerank_folder + "\\check_" + str(i) + '.pkl'
            # write_result_path = rerank_folder + "\\check_" + str(i) + '.test'
            # functionOutTwoPeptidePSM.update_class_variable(1, 1, 0)  # FDR阈值 排序方式 清除反库
            # functionOutTwoPeptidePSM.outputOnlyCrossPSM(search_result, write_result_pkl_path, write_result_path, self.dp.myCFG.separate_FDR)

        rerank_result_pkl_path = self.dp.myCFG.A5_PATH_RESULT_EXPORT + "\\tmp\\" + "\\res-" + FILE_NAME_PSM_RESULT[2][1] + '.pkl'
        rerank_result_path = self.dp.myCFG.A5_PATH_RESULT_EXPORT + "\\tmp\\res-" + FILE_NAME_PSM_RESULT[2][1] + '.txt'

        functionFDR.calculate_cross_result_FDR(search_result)
        functionOutTwoPeptidePSM.update_class_variable(1, 1, 0)# FDR阈值 排序方式 清除反库
        functionOutTwoPeptidePSM.outputOnlyCrossPSM(search_result, rerank_result_pkl_path, rerank_result_path, self.dp.myCFG.separate_FDR)

        rerank_target_result_pkl_path = self.dp.myCFG.A5_PATH_RESULT_EXPORT + "\\tmp\\" + "\\res-target-" + FILE_NAME_PSM_RESULT[2][1] + '.pkl'
        rerank_target_result_path = self.dp.myCFG.A5_PATH_RESULT_EXPORT + "\\res-target-" + FILE_NAME_PSM_RESULT[2][1] + '.txt'
        functionOutTwoPeptidePSM.update_class_variable(self.dp.myCFG.E1_FDR_PSM, 1, 1)# FDR阈值 排序方式 清除反库
        functionOutTwoPeptidePSM.outputOnlyCrossPSM(search_result, rerank_target_result_pkl_path, rerank_target_result_path, self.dp.myCFG.separate_FDR)

    def __captain_show_crosslink_feature_distribution(self, search_result, rerank_folder, target_FDR_threshold=1.0, decoy_FDR_threshold=1.0, input_i=None):

        target_data_rerank_score = []
        target_data_match_score = []
        target_data_match_alpha_score = []
        target_data_match_beta_score = []

        target_data_match_error_sum = []
        target_data_match_error_average = []
        target_data_match_error_var = []

        target_data_match_ion_score = []
        target_data_match_ion_num = []
        target_data_match_ion_intensity = []

        target_data_match_ion_num_percent = []
        target_data_match_ion_intensity_percent = []

        target_data_match_spe_ion_percent = []
        target_data_match_spe_inten_percent = []

        target_data_continue_alpha_score = []
        target_data_continue_beta_score = []

        target_data_spectrum_average_intensity = []

        target_data_alpha_pep_len = []
        target_data_beta_pep_len = []

        target_data_precursor_bias = []
        target_data_delta_score = []
        target_data_pParseNum = []

        target_data_crosslink_delta_score = []

        target_data_mobility = []

        # =======================================================
        decoy_data_rerank_score = []
        decoy_data_match_score = []
        decoy_data_match_alpha_score = []
        decoy_data_match_beta_score = []

        decoy_data_match_error_sum = []
        decoy_data_match_error_average = []
        decoy_data_match_error_var = []

        decoy_data_match_ion_score = []
        decoy_data_match_ion_num = []
        decoy_data_match_ion_intensity = []

        decoy_data_match_ion_num_percent = []
        decoy_data_match_ion_intensity_percent = []

        decoy_data_match_spe_ion_percent = []
        decoy_data_match_spe_inten_percent = []

        decoy_data_continue_alpha_score = []
        decoy_data_continue_beta_score = []

        decoy_data_spectrum_average_intensity = []

        decoy_data_alpha_pep_len = []
        decoy_data_beta_pep_len = []

        decoy_data_precursor_bias = []
        decoy_data_delta_score = []
        decoy_data_pParseNum = []

        decoy_data_crosslink_delta_score = []

        decoy_data_mobility = []

        # =========
        td_data_rerank_score = []
        td_data_match_score = []
        td_data_match_alpha_score = []
        td_data_match_beta_score = []

        td_data_match_error_sum = []
        td_data_match_error_average = []
        td_data_match_error_var = []

        td_data_match_ion_score = []
        td_data_match_ion_num = []
        td_data_match_ion_intensity = []

        td_data_match_ion_num_percent = []
        td_data_match_ion_intensity_percent = []

        td_data_match_spe_ion_percent = []
        td_data_match_spe_inten_percent = []

        td_data_continue_alpha_score = []
        td_data_continue_beta_score = []

        td_data_spectrum_average_intensity = []

        td_data_alpha_pep_len = []
        td_data_beta_pep_len = []

        td_data_precursor_bias = []
        td_data_delta_score = []
        td_data_pParseNum = []

        td_data_crosslink_delta_score = []

        td_data_mobility = []

        for n, one_search_result in enumerate(search_result):
            '''
            # 分数类的特征
            self.rerank_score = rerank_score
            self.match_score = match_score
            self.match_alpha_score = match_alpha_score
            self.match_beta_score = match_beta_score
            # 匹配碎片离子偏差类特征
            self.match_error_average = mathc_error_average
            self.match_error_var = match_error_var
            # 匹配上离子的百分比特征
            self.match_ion_num_percent = match_ion_num_percent
            self.match_ion_intensity_percent = match_ion_intensity_percent
            # 匹配上数目和总强度类特征
            self.match_ion_num = match_ion_num
            self.match_ion_intensity = match_ion_intensity
            # 连续性的特征
            self.continue_peptide_score = continue_peptide_score
            self.continue_beta_score = continue_beta_score
            # 其他和分数以及匹配情况不太相关的特征
            self.delta_score = delta_score
            self.FDR = FDR
            self.pParseNum = pParseNum
            self.precursor_bias = precursor_bias
            '''
            if one_search_result.FDR > target_FDR_threshold and one_search_result.protein_type_math == 1:
                continue
            if one_search_result.protein_type_math == 1 and one_search_result.FDR < target_FDR_threshold:
                # 分数类的特征
                target_data_rerank_score.append(one_search_result.feature.rerank_score)
                target_data_match_score.append(one_search_result.feature.match_score)
                target_data_match_alpha_score.append(one_search_result.feature.match_alpha_score)
                target_data_match_beta_score.append(one_search_result.feature.match_beta_score)
                # 匹配碎片离子偏差类特征
                target_data_match_error_sum.append(one_search_result.feature.match_error_sum)
                target_data_match_error_average.append(one_search_result.feature.match_error_average)
                target_data_match_error_var.append(one_search_result.feature.match_error_var)
                # 匹配上数目和总强度类特征
                target_data_match_ion_score.append(one_search_result.feature.match_ion_score)
                target_data_match_ion_num.append(one_search_result.feature.match_ion_num)
                target_data_match_ion_intensity.append(one_search_result.feature.match_ion_intensity)
                # 匹配上离子的百分比特征
                target_data_match_ion_num_percent.append(one_search_result.feature.match_ion_num_percent)
                target_data_match_ion_intensity_percent.append(one_search_result.feature.match_ion_intensity_percent)
                # 匹配对于谱图的百分比特征
                target_data_match_spe_ion_percent.append(one_search_result.feature.match_spe_ion_percent)
                target_data_match_spe_inten_percent.append(one_search_result.feature.match_spe_inten_percent)
                # 连续性的特征
                target_data_continue_alpha_score.append(one_search_result.feature.continue_alpha_score)
                target_data_continue_beta_score.append(one_search_result.feature.continue_beta_score)
                # 谱图的特征
                target_data_spectrum_average_intensity.append(one_search_result.feature.spectrum_average_intensity)
                # 肽段的特征
                target_data_alpha_pep_len.append(one_search_result.feature.alpha_pep_len)
                target_data_beta_pep_len.append(one_search_result.feature.beta_pep_len)
                # 和肽段母离子相关的
                target_data_precursor_bias.append(one_search_result.feature.precursor_bias)
                target_data_delta_score.append(one_search_result.feature.delta_score)
                target_data_pParseNum.append(one_search_result.feature.pParseNum)
                # 两条肽段差别的分数
                target_data_crosslink_delta_score.append(one_search_result.feature.crosslink_delta_score)
                target_data_mobility.append(one_search_result.feature.mobility)
            elif (one_search_result.protein_type_math == 0) and one_search_result.FDR < decoy_FDR_threshold:
                # 分数类的特征
                decoy_data_rerank_score.append(one_search_result.feature.rerank_score)
                decoy_data_match_score.append(one_search_result.feature.match_score)
                decoy_data_match_alpha_score.append(one_search_result.feature.match_alpha_score)
                decoy_data_match_beta_score.append(one_search_result.feature.match_beta_score)
                # 匹配碎片离子偏差类特征
                decoy_data_match_error_sum.append(one_search_result.feature.match_error_sum)
                decoy_data_match_error_average.append(one_search_result.feature.match_error_average)
                decoy_data_match_error_var.append(one_search_result.feature.match_error_var)
                # 匹配上数目和总强度类特征
                decoy_data_match_ion_score.append(one_search_result.feature.match_ion_score)
                decoy_data_match_ion_num.append(one_search_result.feature.match_ion_num)
                decoy_data_match_ion_intensity.append(one_search_result.feature.match_ion_intensity)
                # 匹配上离子的百分比特征
                decoy_data_match_ion_num_percent.append(one_search_result.feature.match_ion_num_percent)
                decoy_data_match_ion_intensity_percent.append(one_search_result.feature.match_ion_intensity_percent)
                # 匹配对于谱图的百分比特征
                decoy_data_match_spe_ion_percent.append(one_search_result.feature.match_spe_ion_percent)
                decoy_data_match_spe_inten_percent.append(one_search_result.feature.match_spe_inten_percent)
                # 连续性的特征
                decoy_data_continue_alpha_score.append(one_search_result.feature.continue_alpha_score)
                decoy_data_continue_beta_score.append(one_search_result.feature.continue_beta_score)
                # 谱图的特征
                decoy_data_spectrum_average_intensity.append(one_search_result.feature.spectrum_average_intensity)
                # 肽段的特征
                decoy_data_alpha_pep_len.append(one_search_result.feature.alpha_pep_len)
                decoy_data_beta_pep_len.append(one_search_result.feature.beta_pep_len)
                # 和肽段母离子相关的
                decoy_data_precursor_bias.append(one_search_result.feature.precursor_bias)
                decoy_data_delta_score.append(one_search_result.feature.delta_score)
                decoy_data_pParseNum.append(one_search_result.feature.pParseNum)
                # 两条肽段差别的分数
                decoy_data_crosslink_delta_score.append(one_search_result.feature.crosslink_delta_score)
                decoy_data_mobility.append(one_search_result.feature.mobility)
            elif one_search_result.protein_type_math == -1:
                # 分数类的特征
                td_data_rerank_score.append(one_search_result.feature.rerank_score)
                td_data_match_score.append(one_search_result.feature.match_score)
                td_data_match_alpha_score.append(one_search_result.feature.match_alpha_score)
                td_data_match_beta_score.append(one_search_result.feature.match_beta_score)
                # 匹配碎片离子偏差类特征
                td_data_match_error_sum.append(one_search_result.feature.match_error_sum)
                td_data_match_error_average.append(one_search_result.feature.match_error_average)
                td_data_match_error_var.append(one_search_result.feature.match_error_var)
                # 匹配上数目和总强度类特征
                td_data_match_ion_score.append(one_search_result.feature.match_ion_score)
                td_data_match_ion_num.append(one_search_result.feature.match_ion_num)
                td_data_match_ion_intensity.append(one_search_result.feature.match_ion_intensity)
                # 匹配上离子的百分比特征
                td_data_match_ion_num_percent.append(one_search_result.feature.match_ion_num_percent)
                td_data_match_ion_intensity_percent.append(one_search_result.feature.match_ion_intensity_percent)
                # 匹配对于谱图的百分比特征
                td_data_match_spe_ion_percent.append(one_search_result.feature.match_spe_ion_percent)
                td_data_match_spe_inten_percent.append(one_search_result.feature.match_spe_inten_percent)
                # 连续性的特征
                td_data_continue_alpha_score.append(one_search_result.feature.continue_alpha_score)
                td_data_continue_beta_score.append(one_search_result.feature.continue_beta_score)
                # 谱图的特征
                td_data_spectrum_average_intensity.append(one_search_result.feature.spectrum_average_intensity)
                # 肽段的特征
                td_data_alpha_pep_len.append(one_search_result.feature.alpha_pep_len)
                td_data_beta_pep_len.append(one_search_result.feature.beta_pep_len)
                # 和肽段母离子相关的
                td_data_precursor_bias.append(one_search_result.feature.precursor_bias)
                td_data_delta_score.append(one_search_result.feature.delta_score)
                td_data_pParseNum.append(one_search_result.feature.pParseNum)
                # 两条肽段差别的分数
                td_data_crosslink_delta_score.append(one_search_result.feature.crosslink_delta_score)
                td_data_mobility.append(one_search_result.feature.mobility)
        f, axs = plt.subplots(8, 3, figsize=(40, 20))
        # 分数类的特征
        axs[0][0].hist([target_data_rerank_score, decoy_data_rerank_score, td_data_rerank_score], label=['Target', 'Decoy', 'td'], bins=100,
                       stacked=True)
        axs[0][0].title.set_text('rerank_score')

        axs[0][1].hist([target_data_match_score, decoy_data_match_score, td_data_match_score], label=['Target', 'Decoy', 'td'],
                       bins=100, stacked=True)
        axs[0][1].title.set_text('match_score')

        axs[0][2].hist([target_data_match_alpha_score, decoy_data_match_alpha_score, td_data_match_alpha_score], label=['Target', 'Decoy', 'td'],
                       bins=100, stacked=True)
        axs[0][2].title.set_text('match_alpha_score')

        axs[1][0].hist([target_data_match_beta_score, decoy_data_match_beta_score, td_data_match_beta_score], label=['Target', 'Decoy', 'td'],
                       bins=100, stacked=True)
        axs[1][0].title.set_text('match_beta_score')
        # 匹配碎片离子偏差类特征
        axs[1][1].hist([target_data_match_error_average, decoy_data_match_error_average, td_data_match_error_average], label=['Target', 'Decoy', 'td'],
                       bins=100, stacked=True)
        axs[1][1].title.set_text('match_error_average')
        # ==============================
        axs[1][2].hist([target_data_match_error_var, decoy_data_match_error_var, td_data_match_error_var],
                       label=['Target', 'Decoy', 'td'], bins=100, stacked=True,)
        axs[1][2].title.set_text('match_error_var')
        # 匹配上数目和总强度类特征
        axs[2][0].hist([target_data_match_ion_num, decoy_data_match_ion_num, td_data_match_ion_num],
                       label=['Target', 'Decoy', 'td'], bins=100, stacked=True)
        axs[2][0].title.set_text('match_ion_num')

        axs[2][1].hist([target_data_match_ion_intensity, decoy_data_match_ion_intensity, td_data_match_ion_intensity],
                       label=['Target', 'Decoy', 'td'], bins=100, stacked=True)
        axs[2][1].title.set_text('match_ion_intensity')
        # 匹配上离子的百分比特征
        axs[2][2].hist([target_data_match_ion_num_percent, decoy_data_match_ion_num_percent, td_data_match_ion_num_percent],
                       label=['Target', 'Decoy', 'td'], bins=100, stacked=True)
        axs[2][2].title.set_text('match_ion_num_percent')

        axs[3][0].hist([target_data_match_ion_intensity_percent, decoy_data_match_ion_intensity_percent, td_data_match_ion_intensity_percent], label=['Target', 'Decoy', 'td'],
                       bins=100, stacked=True)
        axs[3][0].title.set_text('match_ion_intensity_percent')
        # ==============================
        # 连续性的特征
        axs[3][1].hist([target_data_continue_alpha_score, decoy_data_continue_alpha_score, td_data_continue_alpha_score], label=['Target', 'Decoy', 'td'], bins=100,
                       stacked=True)
        axs[3][1].title.set_text('continue_alpha_score')

        axs[3][2].hist([target_data_continue_beta_score, decoy_data_continue_beta_score, td_data_continue_beta_score], label=['Target', 'Decoy', 'td'], bins=100,
                       stacked=True)
        axs[3][2].title.set_text('continue_beta_score')
        # 其他和分数以及匹配情况不太相关的特
        # 谱图的特征
        axs[4][0].hist([target_data_spectrum_average_intensity, decoy_data_spectrum_average_intensity, td_data_spectrum_average_intensity],
                       label=['Target', 'Decoy', 'td'], bins=100, stacked=True)
        axs[4][0].title.set_text('spectrum_average_intensity')
        # 肽段的特征
        axs[4][1].hist([target_data_alpha_pep_len, decoy_data_alpha_pep_len, td_data_alpha_pep_len],
                       label=['Target', 'Decoy', 'td'], bins=100, stacked=True)
        axs[4][1].title.set_text('alpha_pep_len')

        axs[4][2].hist([target_data_beta_pep_len, decoy_data_beta_pep_len, td_data_beta_pep_len],
                       label=['Target', 'Decoy', 'td'], bins=100, stacked=True)
        axs[4][2].title.set_text('beta_pep_len')
        # ==============================
        # 和肽段母离子相关的
        axs[5][0].hist([target_data_precursor_bias, decoy_data_precursor_bias, td_data_precursor_bias], label=['Target', 'Decoy', 'td'],
                       bins=100, stacked=True)
        axs[5][0].title.set_text('precursor_bias')

        axs[5][1].hist([target_data_delta_score, decoy_data_delta_score, td_data_delta_score], label=['Target', 'Decoy', 'td'],
                       bins=100, stacked=True)
        axs[5][1].title.set_text('delta_score')

        # 两条肽段差别的分数
        axs[5][2].hist([target_data_crosslink_delta_score, decoy_data_crosslink_delta_score, td_data_crosslink_delta_score], label=['Target', 'Decoy', 'td'],
                       bins=100, stacked=True)
        axs[5][2].title.set_text('crosslink_delta_score')


        axs[6][0].hist([target_data_pParseNum, decoy_data_pParseNum, td_data_pParseNum], label=['Target', 'Decoy', 'td'],
                       bins=100, stacked=True)
        axs[6][0].title.set_text('pParseNum')

        axs[6][1].hist([target_data_match_spe_ion_percent, decoy_data_match_spe_ion_percent, td_data_match_spe_ion_percent], label=['Target', 'Decoy', 'td'],
                       bins=100, stacked=True)
        axs[6][1].title.set_text('match_match_spe_ion_percent')

        axs[6][2].hist([target_data_match_spe_inten_percent, decoy_data_match_spe_inten_percent, td_data_match_spe_inten_percent], label=['Target', 'Decoy', 'td'],
                       bins=100, stacked=True)
        axs[6][2].title.set_text('match_spe_inten_percent')

        axs[7][0].hist([target_data_match_error_sum, decoy_data_match_error_sum, td_data_match_error_sum], label=['Target', 'Decoy', 'td'],
                       bins=100, stacked=True)
        axs[7][0].title.set_text('match_error_sum')

        axs[7][1].hist([target_data_match_ion_score, decoy_data_match_ion_score, td_data_match_ion_score], label=['Target', 'Decoy', 'td'],
                       bins=100, stacked=True)
        axs[7][1].title.set_text('match_ion_score')

        axs[7][2].hist([target_data_mobility, decoy_data_mobility, td_data_mobility], label=['Target', 'Decoy', 'td'],
                       bins=100, stacked=True)
        if self.dp.myCFG.E3_MOBILITY_RERANK:
            axs[7][2].title.set_text('mobility')
        else:
            axs[7][2].title.set_text('mobility_closed')
        # ==============================

        axs[0][0].legend()
        axs[0][1].legend()
        axs[0][2].legend()

        axs[1][0].legend()
        axs[1][1].legend()
        axs[1][2].legend()

        axs[2][0].legend()
        axs[2][1].legend()
        axs[2][2].legend()

        axs[3][0].legend()
        axs[3][1].legend()
        axs[3][2].legend()

        axs[4][0].legend()
        axs[4][1].legend()
        axs[4][2].legend()

        axs[5][0].legend()
        axs[5][1].legend()
        axs[5][2].legend()

        axs[6][0].legend()
        axs[6][1].legend()
        axs[6][2].legend()

        axs[7][0].legend()
        axs[7][1].legend()
        axs[7][2].legend()

        if input_i is None:
            plt.savefig(rerank_folder + "\\" + 'rerank_figure.png', format='png')
        else:
            plt.savefig(rerank_folder + "\\" + 'rerank_figure-' + str(input_i) + '.png', format='png')

    def __captain_get_crosslink_test_dataset(self, search_result, test_file):
        # 这是把所有数据都归一化写入到文件中
        with open(test_file, 'w') as f:
            for one_search_index, one_search_result in enumerate(search_result):
                if len(one_search_result.protein_type) == 0:
                    break
                if one_search_result.protein_type_math == 1:
                    flag = '1'
                elif one_search_result.protein_type_math == 0:
                    flag = '-1'
                elif one_search_result.protein_type_math == -1:
                    flag = '1'
                # write_str += (str(flag) +
                #               # 分数类的特征
                #               " 1:" + str(one_search_result.feature.rerank_score) +
                #               " 2:" + str(one_search_result.feature.match_score) +
                #               " 3:" + str(one_search_result.feature.match_alpha_score) +
                #               " 4:" + str(one_search_result.feature.match_beta_score) +
                #               # 匹配碎片离子偏差类特征
                #               " 5:" + str(one_search_result.feature.match_error_sum) +
                #               " 6:" + str(one_search_result.feature.match_error_average) +
                #               " 7:" + str(one_search_result.feature.match_error_var) +
                #               # 匹配上数目和总强度类特征
                #               " 8:" + str(one_search_result.feature.match_ion_score) +
                #               " 9:" + str(one_search_result.feature.match_ion_num_percent) +
                #               " 10:" + str(one_search_result.feature.match_ion_intensity_percent) +
                #               # 匹配上离子的百分比特征
                #               " 11:" + str(one_search_result.feature.match_ion_num) +
                #               " 12:" + str(one_search_result.feature.match_ion_intensity) +
                #               # 匹配对于谱图的百分比特征
                #               " 13:" + str(one_search_result.feature.match_spe_ion_percent) +
                #               " 14:" + str(one_search_result.feature.match_spe_inten_percent) +
                #               # 连续性的特征
                #               " 15:" + str(one_search_result.feature.continue_alpha_score) +
                #               " 16:" + str(one_search_result.feature.continue_beta_score) +
                #               # 谱图的特征
                #               " 17:" + str(one_search_result.feature.spectrum_average_intensity) +
                #               # 肽段的特征
                #               " 18:" + str(one_search_result.feature.alpha_pep_len) +
                #               " 19:" + str(one_search_result.feature.beta_pep_len) +
                #               # 和肽段母离子相关的
                #               " 20:" + str(one_search_result.feature.precursor_bias) +
                #               " 21:" + str(one_search_result.feature.delta_score) +
                #               " 22:" + str(one_search_result.feature.pParseNum) +
                #               # 两条肽段差别的分数
                #               " 23:" + str(one_search_result.feature.crosslink_delta_score) +
                #               "\n")
                if self.dp.myCFG.E3_MOBILITY_RERANK:
                    write_str = (str(flag) +
                               " 1:" + str(one_search_result.feature.match_score) +
                               # " 3:" + str(one_search_result.feature.match_alpha_score) +
                               # " 4:" + str(one_search_result.feature.match_beta_score) +

                               # " 5:" + str(one_search_result.feature.match_error_sum) +
                               " 6:" + str(one_search_result.feature.match_error_average) +
                               " 7:" + str(one_search_result.feature.match_error_var) +
                               " 8:" + str(one_search_result.feature.match_ion_score) +
                               " 9:" + str(one_search_result.feature.match_ion_num_percent) +
                               # " 10:" + str(one_search_result.feature.match_ion_intensity_percent) +
                               " 13:" + str(one_search_result.feature.match_spe_ion_percent) +
                               " 14:" + str(one_search_result.feature.match_spe_inten_percent) +
                               # " 15:" + str(one_search_result.feature.continue_alpha_score) +
                               " 16:" + str(one_search_result.feature.continue_beta_score) +
                               # " 18:" + str(one_search_result.feature.alpha_pep_len) +
                               # " 19:" + str(one_search_result.feature.beta_pep_len) +
                               # " 20:" + str(one_search_result.feature.precursor_bias) +
                               " 21:" + str(one_search_result.feature.delta_score) +
                               # " 22:" + str(one_search_result.feature.pParseNum) +
                               # " 23:" + str(one_search_result.feature.crosslink_delta_score) +
                               " 24:" + str(one_search_result.feature.rerank_score) +
                               " 25:" + str(one_search_result.feature.mobility) +
                               "\n")
                else:
                    write_str = (str(flag) +
                               " 1:" + str(one_search_result.feature.match_score) +

                               " 3:" + str(one_search_result.feature.match_alpha_score) +
                               " 4:" + str(one_search_result.feature.match_beta_score) +

                               # " 5:" + str(one_search_result.feature.match_error_sum) +
                               " 6:" + str(one_search_result.feature.match_error_average) +
                               " 7:" + str(one_search_result.feature.match_error_var) +
                               # " 8:" + str(one_search_result.feature.match_ion_score) +
                               " 9:" + str(one_search_result.feature.match_ion_num_percent) +
                               # " 10:" + str(one_search_result.feature.match_ion_intensity_percent) +
                               " 13:" + str(one_search_result.feature.match_spe_ion_percent) +
                               " 14:" + str(one_search_result.feature.match_spe_inten_percent) +
                               # " 15:" + str(one_search_result.feature.continue_alpha_score) +
                               " 16:" + str(one_search_result.feature.continue_beta_score) +
                               # " 18:" + str(one_search_result.feature.alpha_pep_len) +
                               # " 19:" + str(one_search_result.feature.beta_pep_len) +
                               # " 20:" + str(one_search_result.feature.precursor_bias) +
                               " 21:" + str(one_search_result.feature.delta_score) +
                               # " 22:" + str(one_search_result.feature.pParseNum) +
                               # " 23:" + str(one_search_result.feature.crosslink_delta_score) +
                               " 24:" + str(one_search_result.feature.rerank_score) +
                               "\n")
                f.write(write_str)