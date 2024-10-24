# -*- mode: python ; coding: utf-8 -*-..
import sys
import multiprocessing
import psutil
from MSStaff import CStaff
import faulthandler
if __name__ == "__main__":
    memory_data = psutil.virtual_memory()
    used_memory = float(memory_data.used) / 1024 / 1024 / 1024
    print("start_used_memory : %.4f GB" % (used_memory))
    multiprocessing.freeze_support()
    staff = CStaff()
    staff.start(sys.argv)
    sys.exit(10000)
