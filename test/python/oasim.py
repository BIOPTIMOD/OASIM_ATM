#!/usr/bin/python3

import ctypes as ct
import numpy as np


class oasim_lib:
    def __init__(self, lib_path, config_path, lat, lon):
        self._lib = ct.cdll.LoadLibrary(lib_path)
        self._init_lib = self._lib.py_init_lib
        self._init_lib.restype = ct.c_void_p
        self._init_lib.argtypes = [ct.c_int,
                                   ct.c_int,
                                   ct.c_char_p,
                                   np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                   np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS")]

        self._finalize_lib = self._lib.py_finalize_lib
        self._finalize_lib.restype = None
        self._finalize_lib.argtypes = [ct.c_void_p]

        error = ct.c_bool()
        if lat.shape != lon.shape or len(lat.shape) != 1:
            raise ValueError("The lat and lon mesh are incompatible")

        self._ptr = self._init_lib(len(config_path),
                                   lat.shape[0],
                                   config_path.encode("utf-8"),
                                   lat,
                                   lon)

        if not self._ptr:
            raise FileNotFoundError("Corrupted file: {}".format(config_path))
    

    def __del__(self):
        if self._ptr:
            self._finalize_lib(self._ptr)
        


class calc_unit:
    def __init__(self, p_size, lib):
        self._lib = lib
        self._init_calc = self._lib._lib.py_init_calc
        self._init_calc.restype = ct.c_void_p
        self._init_calc.argtypes = [ct.c_int,
                                    ct.c_void_p]

        self._monrad = self._lib._lib.py_monrad
        self._monrad.restype = ct.c_bool
        self._monrad.argtypes = [ct.c_void_p,
                                 ct.c_int,
                                 ct.c_int,
                                 np.ctypeslib.ndpointer(ct.c_int, flags="F_CONTIGUOUS"),
                                 ct.c_int,
                                 ct.c_int,
                                 ct.c_double,
                                 ct.c_double,
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS"),
                                 np.ctypeslib.ndpointer(ct.c_double, flags="F_CONTIGUOUS")]
                                 

        self._finalize_calc = self._lib._lib.py_finalize_calc
        self._finalize_calc.restype = None
        self._finalize_calc.argtypes = [ct.c_void_p]
        
        self._ptr = self._init_calc(p_size, self._lib._ptr)


    def monrad(self, points, iyr, iday, sec_b, sec_e,
               sp, msl, ws10, tco3, t2m, d2m, tcc, tclw, cdrem,
               taua, asymp, ssalb):
        p_size, rows = taua.shape
        edout = np.zeros((p_size, rows), dtype=np.float64, order="F")
        esout = np.zeros((p_size, rows), dtype=np.float64, order="F")

        error = self._monrad(self._ptr, p_size, rows, 
                             np.asfortranarray(points, dtype=np.int32),
                             iyr, iday, sec_b, sec_e,
                             np.asfortranarray(sp, dtype=np.float64), 
                             np.asfortranarray(msl, dtype=np.float64), 
                             np.asfortranarray(ws10, dtype=np.float64), 
                             np.asfortranarray(tco3, dtype=np.float64), 
                             np.asfortranarray(t2m, dtype=np.float64), 
                             np.asfortranarray(d2m, dtype=np.float64), 
                             np.asfortranarray(tcc, dtype=np.float64),
                             np.asfortranarray(tclw, dtype=np.float64), 
                             np.asfortranarray(cdrem, dtype=np.float64),
                             np.asfortranarray(taua, dtype=np.float64), 
                             np.asfortranarray(asymp, dtype=np.float64), 
                             np.asfortranarray(ssalb, dtype=np.float64), 
                             edout, esout)
        
        if error:
            raise ValueError("Error in monrad computation.")

        return edout, esout


    def __del__(self):
        if self._ptr:
            self._finalize_calc(self._ptr)


lat = np.array([42.9375], dtype=np.float64, order='F')
lon = np.array([7.2916665], dtype=np.float64, order='F')
olib = oasim_lib("../../OASIMlib/liboasim-py.so", "config.yaml", lat, lon)
cunit = calc_unit(1, olib)

points = np.array([1])
iyr = 2019
iday = 1
sec_b = 28800.0
sec_e = 36000.0
sp = np.array([102377.13])
msl = np.array([102410.06])
ws10 = np.array([10.0])
tco3 = np.array([0.0065])
t2m = np.array([286.96487])
d2m = np.array([284.10672])
tcc = np.array([12.03943])
tclw = np.array([0.05])
cdrem = np.array([10.042508])
taua = np.array([[0.12861817,0.11693237,0.11303711,0.10914185,0.10524658,0.10135132,
    0.09745605,0.0934308,0.08906767,0.08519087,0.08188277,0.079188,
    0.07697173,0.07506739,0.07332103,0.07159076,0.06982587,0.06804233,
    0.0645105,0.05960191,0.05457829,0.05089346,0.0479044,0.04534442,
    0.04311129,0.04106918,0.03914813,0.0373441,0.03567282,0.03332167,
    0.02901247,0.01838189,0.00623265]])
asymp = np.array([[0.78112185,0.76577055,0.76065344,0.7555364,0.75041926,0.74530214,
    0.7401851,0.7352586,0.7308154,0.72576404,0.71973574,0.7124992,
    0.7045832,0.69671303,0.68963724,0.68404067,0.6799671,0.6770385,
    0.6733678,0.6695107,0.65929943,0.64360243,0.62599105,0.6091321,
    0.59477234,0.58193725,0.5697906,0.5578721,0.54593724,0.5280586,
    0.49310532,0.41066745,0.3164527]])
ssalb = np.array([[0.9653854,0.9685926,0.96966165,0.9707307,0.9717998,0.97286886,
    0.9739379,0.9750599,0.97628933,0.97734255,0.97821194,0.9788913,
    0.9793933,0.97973263,0.9799186,0.9799566,0.9798687,0.9796886,
    0.9791781,0.9784721,0.9785546,0.9787444,0.97678363,0.96918947,
    0.95303714,0.9326887,0.91252154,0.897008,0.88944095,0.88677096,
    0.89569014,0.89257544,0.88901585]])    

edout, esout = cunit.monrad(points, iyr, iday, sec_b, sec_e,
                            sp, msl, ws10, tco3, t2m, d2m, tcc, tclw, cdrem,
                            taua, asymp, ssalb)

print(edout)
print(esout)

del cunit
del olib
