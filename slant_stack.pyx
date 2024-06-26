# cython: language_level=3
# distutils: language=c++
# cython: c_string_type=unicode, c_string_encoding=utf8

# Import necessary Cython declarations
from libcpp.complex cimport complex
from libcpp.vector cimport vector

ctypedef complex[double] Complex
ctypedef vector[Complex] CArray

cdef extern from "slant_stack_src.cpp":
    CArray fft_src(CArray x_t)


# Define Python wrapper functions
def slant_stack(CArray x_t):
    fft_src(x_t)
    print(x_t)
    return x_t