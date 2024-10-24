# -*- mode: python ; coding: utf-8 -*-..

from MSFlow import CFlow1, CFlow2, CFlow3
from MSFunction import CFunctionConfigIO
from MSSysterm import EXPIRATION_TIME, FILE_NAME_CONFIG
from MSData import CDataPack, CConfig
from MSLogging import INFO_TO_USER_Staff
from MSLogging import logGetError

import datetime
import time

class CStaff:

    def __init__(self):
        self.dp = CDataPack()

    def start(self, argv):
        self.__checktime()
        self.__runflow(argv)

    def __checktime(self):
        dateNow = datetime.datetime.now()
        dateDead = datetime.datetime(EXPIRATION_TIME['Year'], EXPIRATION_TIME['Month'], EXPIRATION_TIME['Day'], 23, 59)
        n_day = (dateDead - dateNow).days
        if n_day < 0:
            print(INFO_TO_USER_Staff[0])
            exit(-1)
        elif n_day < 7:
            print(INFO_TO_USER_Staff[1])
        else:
            print(INFO_TO_USER_Staff[2])
            print(dateNow)

    def __runflow(self, argv):
        if len(argv) == 1:

            start_time = time.time()

            config = CConfig()
            functionConfigIO = CFunctionConfigIO()

            functionConfigIO.config2file(FILE_NAME_CONFIG, config)

            end_time = time.time()
            print("[Info] Task finish in {}s".format(end_time - start_time))

        elif len(argv) == 2:

            start_time = time.time()
            print("Config file", argv)

            functionConfigIO = CFunctionConfigIO()
            functionConfigIO.file2config(argv[1], self.dp.myCFG)

            if self.dp.myCFG.D13_TYPE_FLOW == 1:

                flow1 = CFlow1(self.dp)
                flow1.run()

            elif self.dp.myCFG.D13_TYPE_FLOW == 2:

                flow2 = CFlow2(self.dp)
                flow2.run()

            elif self.dp.myCFG.D13_TYPE_FLOW == 3:
                flow3 = CFlow3(self.dp)
                flow3.run()

            else:

                pass

            end_time = time.time()
            print("[Info] Task finish in {} s . ".format(end_time - start_time))
            print("[Info] Local time {} ".format(time.asctime(time.localtime(time.time()))))

