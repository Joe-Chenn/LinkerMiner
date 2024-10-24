import ctypes


class ION_INDEX_MZ(ctypes.Structure):
    _fields_ = [("pep_num", ctypes.c_long),
                ("pep", ctypes.POINTER(ctypes.c_long)),
                ("pep_size", ctypes.c_long)]
