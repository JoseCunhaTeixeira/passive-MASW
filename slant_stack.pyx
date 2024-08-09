# distutils: language=c++
# distutils: libraries=fftw3

from libcpp.vector cimport vector
from libcpp.complex cimport complex
cimport cython

cdef extern from "slant_stack_src.cpp":
    vector[vector[double]] makeFV_src(vector[vector[double]] XT, double si, vector[double] offset, double vmin, double vmax, double dv, double fmax)

cdef class FVResult:
    cdef vector[vector[double]] FV
    cdef vector[double] vs
    cdef vector[double] fs

def slant_stack(vector[vector[double]] XT, double si, offset, double vmin, double vmax, double dv, double fmax):
    cdef vector[vector[double]] XT_vec = XT
    cdef vector[double] offset_vec = offset
    cdef FVResult result = FVResult()
    result.FV = makeFV_src(XT_vec, si, offset_vec, vmin, vmax, dv, fmax)
    x = vmin
    for i in range(int((float(vmax) + float(dv) - float(vmin)) // float(dv))):
        result.vs.push_back(x)
        x = x + dv
    result.fs = [i / (si * len(XT[0])) for i in range(len(XT[0])) if i / (si * len(XT[0])) <= fmax]
    return result.FV, result.vs, result.fs
