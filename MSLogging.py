# -*- mode: python ; coding: utf-8 -*-
import logging
import sys
from MSSysterm import LOG_LEVEL

myLogPath = 'eLink.log'


# 也许向某个日志文件里写东西
def logToUser(strInfo):
    # if os.access(myLogPath, os.W_OK):  # 当文件被excel打开时，这个东东没啥用

    try:
        print(strInfo)
        f_w = open(myLogPath, 'a', encoding='utf8')
        f_w.write(strInfo + '\n')
        f_w.close()
    except IOError:
        print("eLink.log is opened! Please close it and run the program again!")
        sys.exit(0)


def logGetError(info):
    print(info)

    logging.basicConfig(filename='eLink.log',
                        filemode='a',
                        format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s'
                        )
    logging.error(info)
    sys.exit(0)


def logGetWarning(info):
    print(info)

    logging.basicConfig(filename='eLink.log',
                        filemode='a',
                        format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s'
                        )
    logging.warning(info)


INFO_TO_USER_Staff = (
    '\n[eLink] eLink is expired! Please send e-mail to  for the new version.',
    '\n[eLink] Warning! The current license will expired in 7 days. Please send e-mail to for the new version.',
    '\n[eLink] Starting...',
    '\n[eLink] Finished!',
    '\n[eLink] Writing config file...',)

INFO_TO_USER_Flow1 = (
    '\n[eLink] [Info] Flow finish in ',
    ''
)

SOFTWARE_NAME = "eLink"
INFO_TO_USER_FunctionComposition = (

    '\n[{}] <Function Composition> Wrong format: '.format(SOFTWARE_NAME),
    '\n[{}] <Function Composition> Make sure the format of this modifications is correct!'.format(SOFTWARE_NAME),
    '\n[{}] <Function Composition> Make sure the format of this glyco is correct!'.format(SOFTWARE_NAME),
    '\n[{}] <Function Composition> Make sure the format of this linker is correct!'.format(SOFTWARE_NAME),
)

from enum import Enum
import datetime


class LogLevel(Enum):
    ERROR = 0
    WARNING = 1
    INFO = 2
    DEBUG = 3


class Logger:

    @staticmethod
    def __get_time():
        return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    @staticmethod
    def info(message):
        if LOG_LEVEL >= LogLevel.INFO.value:
            print(f"{Logger.__get_time()} [INFO] {message}")

    @staticmethod
    def warning(message):
        if LOG_LEVEL >= LogLevel.WARNING.value:
            print(f"{Logger.__get_time()} [WARNING] {message}")

    @staticmethod
    def error(message):
        if LOG_LEVEL >= LogLevel.ERROR.value:
            print(f"{Logger.__get_time()} [ERROR] {message}")

    @staticmethod
    def debug(message):
        if LOG_LEVEL >= LogLevel.DEBUG.value:
            print(f"{Logger.__get_time()} [DEBUG] {message}")
