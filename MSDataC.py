# -*- mode: python ; coding: utf-8 -*-
from ctypes import *

class ION_INDEX_PEP(Structure):
    _fields_ = [("pep_index_num", c_long)]

class ION_INDEX_MZ(Structure):
    _fields_ = [("pep_num", c_long)]

class ION_INDEX_PEP_LOCK(Structure):
    _fields_ = [("pep_index_num", c_long)]

class ION_INDEX_MZ_LOCK(Structure):
    _fields_ = [("pep_num", c_long)]
