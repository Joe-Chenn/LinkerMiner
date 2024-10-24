# -*- mode: python ; coding: utf-8 -*-..
from MSTask import CTaskCheck, CTaskFillConfig, CTaskFilterTG
from MSTask import CTaskCreatePeptideIndex, CTaskCreateIonIndex
from MSTask import CTaskPreProcessSpectrum, CTaskSearch, CTaskRerank, CTaskGetIdPSM
from MSLogging import logToUser, logGetError
from MSLogging import INFO_TO_USER_Flow1
from MSSysterm import FOLDER_INDEX

import time


class CFlow1:

    # 只搜一种形式，通过D14_TYPE_PEPTIDE和D15_TYPE_LINK控制

    def __init__(self, inputDP):
        self.dp = inputDP

    def run(self):

        logToUser(INFO_TO_USER_Flow1[0])
        taskCheck = CTaskCheck(self.dp)
        taskCheck.work()

        logToUser(INFO_TO_USER_Flow1[1])
        taskLoadFile = CTaskFillConfig(self.dp)
        taskLoadFile.work()

        time_1 = time.time()

        taskPreProcessSpectrum = CTaskPreProcessSpectrum(self.dp)
        taskPreProcessSpectrum.work()

        time_2 = time.time()
        if self.dp.myCFG.S1_SKIP_CREATE_PEPTIDE_INDEX:
            taskCreatePeptideIndex = CTaskCreatePeptideIndex(self.dp)
            taskCreatePeptideIndex.work()

        time_3 = time.time()
        if self.dp.myCFG.S2_SKIP_CREATE_ION_INDEX:
            taskCreateIonIndex = CTaskCreateIonIndex(self.dp)
            taskCreateIonIndex.work()
        time_4 = time.time()
        if self.dp.myCFG.S3_SKIP_SEARCH:
            taskSearch = CTaskSearch(self.dp)
            taskSearch.work()

        time_5 = time.time()
        if self.dp.myCFG.S4_SKIP_RERANK:
            taskRerank = CTaskRerank(self.dp)
            taskRerank.work()

        time_6 = time.time()

        print("[Info] Load spectra using time : %.4f s." % (time_2 - time_1))
        print("[Info] Create " + str(FOLDER_INDEX[self.dp.myCFG.D14_TYPE_PEPTIDE][
                                         self.dp.myCFG.D15_TYPE_LINK]) + " peptide index time : %.4f s." % (
                          time_3 - time_2))
        print("[Info] Create " + str(
            FOLDER_INDEX[self.dp.myCFG.D14_TYPE_PEPTIDE][self.dp.myCFG.D15_TYPE_LINK]) + " ion index time : %.4f s." % (
                          time_4 - time_3))
        print("[Info] Search " + str(
            FOLDER_INDEX[self.dp.myCFG.D14_TYPE_PEPTIDE][self.dp.myCFG.D15_TYPE_LINK]) + " time : %.4f s." % (
                          time_5 - time_4))
        print("[Info] Rerank " + str(
            FOLDER_INDEX[self.dp.myCFG.D14_TYPE_PEPTIDE][self.dp.myCFG.D15_TYPE_LINK]) + " time : %.4f s." % (
                          time_6 - time_5))


class CFlow2:

    # 只搜单肽，cross和loop

    def __init__(self, inputDP):
        self.dp = inputDP

    def run(self):
        id_psm_name = []
        id_psm_scan = []

        taskCheck = CTaskCheck(self.dp)
        taskCheck.work()
        taskLoadFile = CTaskFillConfig(self.dp)
        taskLoadFile.work()
        # logToUser(INFO_TO_USER_Flow1[2])
        time_1 = time.time()

        taskPreProcessSpectrum = CTaskPreProcessSpectrum(self.dp)
        taskPreProcessSpectrum.work()

        time_2 = time.time()

        if self.dp.myCFG.S1_SKIP_CREATE_PEPTIDE_INDEX:
            taskCreatePeptideIndex = CTaskCreatePeptideIndex(self.dp)
            taskCreatePeptideIndex.work()

        time_3 = time.time()

        taskCreateIonIndex = CTaskCreateIonIndex(self.dp)

        # self.dp.myCFG.D14_TYPE_PEPTIDE = 1
        # self.dp.myCFG.D15_TYPE_LINK = 1
        #
        # if self.dp.myCFG.S2_SKIP_CREATE_ION_INDEX:
        #     taskCreateIonIndex.work()
        #
        # time_4 = time.time()

        self.dp.myCFG.D14_TYPE_PEPTIDE = 1
        self.dp.myCFG.D15_TYPE_LINK = 2
        if self.dp.myCFG.S2_SKIP_CREATE_ION_INDEX:
            taskCreateIonIndex.work()

        time_5 = time.time()

        self.dp.myCFG.D14_TYPE_PEPTIDE = 2
        self.dp.myCFG.D15_TYPE_LINK = 1
        if self.dp.myCFG.S2_SKIP_CREATE_ION_INDEX:
            taskCreateIonIndex.work()

        time_6 = time.time()

        taskSearch = CTaskSearch(self.dp)
        taskRerank = CTaskRerank(self.dp)
        taskGetIdPSM = CTaskGetIdPSM(self.dp)

        time_9 = time.time()

        # self.dp.myCFG.D14_TYPE_PEPTIDE = 1
        # self.dp.myCFG.D15_TYPE_LINK = 2
        # if self.dp.myCFG.S3_SKIP_SEARCH:
        #     taskSearch.work(id_psm_name, id_psm_scan)
        #
        # time_10 = time.time()
        # if self.dp.myCFG.S4_SKIP_RERANK:
        #     taskRerank.work()

        time_11 = time.time()

        taskGetIdPSM.work(id_psm_name, id_psm_scan)

        time_12 = time.time()

        self.dp.myCFG.D14_TYPE_PEPTIDE = 2
        self.dp.myCFG.D15_TYPE_LINK = 1
        if self.dp.myCFG.S3_SKIP_SEARCH:
            taskSearch.work(id_psm_name, id_psm_scan)

        taskFilterTG = CTaskFilterTG(self.dp)
        taskFilterTG.work()

        time_13 = time.time()
        if self.dp.myCFG.S4_SKIP_RERANK:
            taskRerank.work()

        time_14 = time.time()

        print("[Info] Load spectra using time : %.4f s." % (time_2 - time_1))
        print("[Info] Create peptide index time : %.4f s." % (time_3 - time_2))
        # print("[Info] Create single peptide peptide index time : %.4f s."%(time_4 - time_3))
        print("[Info] Create loop-link peptide index time : %.4f s." % (time_5 - time_3))
        print("[Info] Creating cross-link peptide peptide index time : %.4f s." % (time_6 - time_5))
        print("[Info] Creating loop-link ion index time : %.4f s." % (time_5 - time_3))
        print("[Info] Creating cross-link ion index time : %.4f s." % (time_6 - time_5))
        # print("[Info] Searching single peptide time : %.4f s."%(time_9 - time_6))
        # print("[Info] Reranking single peptide time : %.4f s."%(time_8 - time_6))
        # print("[Info] Reading single peptide result time : %.4f s"%(time_9 - time_8))
        print("[Info] Searching loop-link time : %.4f s." % (time_10 - time_9))
        print("[Info] Reranking loop-link time : %.4f s." % (time_11 - time_10))
        print("[Info] Reading loop-link result time : %.4f s" % (time_12 - time_11))
        print("[Info] Searching cross-link time : %.4f s." % (time_13 - time_12))
        print("[Info] Reranking cross-link time : %.4f s." % (time_14 - time_13))

class CFlow3:
    def __init__(self, inputDP):
        self.dp = inputDP

    def run(self):
        id_psm_name = []
        id_psm_scan = []

        taskCheck = CTaskCheck(self.dp)
        taskCheck.work()
        taskLoadFile = CTaskFillConfig(self.dp)
        taskLoadFile.work()
        # logToUser(INFO_TO_USER_Flow1[2])
        time_1 = time.time()

        taskPreProcessSpectrum = CTaskPreProcessSpectrum(self.dp)
        taskPreProcessSpectrum.work()

        time_2 = time.time()

        if self.dp.myCFG.S1_SKIP_CREATE_PEPTIDE_INDEX:
            taskCreatePeptideIndex = CTaskCreatePeptideIndex(self.dp)
            taskCreatePeptideIndex.work()

        time_3 = time.time()

        taskCreateIonIndex = CTaskCreateIonIndex(self.dp)

        # self.dp.myCFG.D14_TYPE_PEPTIDE = 1
        # self.dp.myCFG.D15_TYPE_LINK = 1
        #
        # if self.dp.myCFG.S2_SKIP_CREATE_ION_INDEX:
        #     taskCreateIonIndex.work()
        #
        # time_4 = time.time()

        self.dp.myCFG.D14_TYPE_PEPTIDE = 1
        self.dp.myCFG.D15_TYPE_LINK = 2
        if self.dp.myCFG.S2_SKIP_CREATE_ION_INDEX:
            taskCreateIonIndex.work()

        time_5 = time.time()

        self.dp.myCFG.D14_TYPE_PEPTIDE = 2
        self.dp.myCFG.D15_TYPE_LINK = 1
        if self.dp.myCFG.S2_SKIP_CREATE_ION_INDEX:
            taskCreateIonIndex.work()

        time_6 = time.time()

        taskSearch = CTaskSearch(self.dp)
        taskRerank = CTaskRerank(self.dp)
        taskGetIdPSM = CTaskGetIdPSM(self.dp)

        time_9 = time.time()

        # self.dp.myCFG.D14_TYPE_PEPTIDE = 1
        # self.dp.myCFG.D15_TYPE_LINK = 2
        # if self.dp.myCFG.S3_SKIP_SEARCH:
        #     taskSearch.work(id_psm_name, id_psm_scan)

        # time_10 = time.time()
        # if self.dp.myCFG.S4_SKIP_RERANK:
        #     taskRerank.work()
        #
        # time_11 = time.time()
        #
        # taskGetIdPSM.work(id_psm_name, id_psm_scan)

        time_12 = time.time()

        self.dp.myCFG.D14_TYPE_PEPTIDE = 2
        self.dp.myCFG.D15_TYPE_LINK = 1
        if self.dp.myCFG.S3_SKIP_SEARCH:
            taskSearch.work(id_psm_name, id_psm_scan)

        taskFilterTG = CTaskFilterTG(self.dp)
        taskFilterTG.work()

        time_13 = time.time()
        if self.dp.myCFG.S4_SKIP_RERANK:
            taskRerank.work()

        time_14 = time.time()

        print("[Info] Load spectra using time : %.4f s." % (time_2 - time_1))
        print("[Info] Create peptide index time : %.4f s." % (time_3 - time_2))
        # print("[Info] Create single peptide peptide index time : %.4f s."%(time_4 - time_3))
        print("[Info] Create loop-link peptide index time : %.4f s." % (time_5 - time_3))
        print("[Info] Creating cross-link peptide peptide index time : %.4f s." % (time_6 - time_5))
        print("[Info] Creating loop-link ion index time : %.4f s." % (time_5 - time_3))
        print("[Info] Creating cross-link ion index time : %.4f s." % (time_6 - time_5))
        # print("[Info] Searching single peptide time : %.4f s."%(time_9 - time_6))
        # print("[Info] Reranking single peptide time : %.4f s."%(time_8 - time_6))
        # print("[Info] Reading single peptide result time : %.4f s"%(time_9 - time_8))
        print("[Info] Searching loop-link time : %.4f s." % (time_10 - time_9))
        print("[Info] Reranking loop-link time : %.4f s." % (time_11 - time_10))
        print("[Info] Reading loop-link result time : %.4f s" % (time_12 - time_11))
        print("[Info] Searching cross-link time : %.4f s." % (time_13 - time_12))
        print("[Info] Reranking cross-link time : %.4f s." % (time_14 - time_13))
