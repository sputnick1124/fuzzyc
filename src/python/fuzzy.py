from ctypes import *

class Rule(Structure):
    _fields_ = [
            ('input', POINTER(POINTER(c_double))),
            ('num_in', c_int),
            ('output', POINTER(POINTER(c_double))),
            ('num_out', c_int)
            ]

class Fis(Structure):
    _fields_ = [
            ('num_rule', c_int),
            ('rule_list', POINTER(POINTER(Rule)))
            ]

def fis_create(params, 
fuzzylib = CDLL('./fuzzy.so')

_fis_create = fuzzylib.fis_create
_fis_create.argtypes = [POINTER(c_double), c_int, c_int, c_int, POINTER(POINTER(c_int)), POINTER(c_int),
        POINTER(c_int)]
_fis_create.restype = POINTER(Fis)

