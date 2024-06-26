import subprocess
import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize



# if not os.path.exists("./lib/"):
#    os.makedirs("./lib/")



extensions = [
    Extension("slant_stack", 
                ["slant_stack.pyx"],
                language="c++",
                extra_compile_args=["-O3"],
                extra_link_args=["-O3"]),
]

setup(
    ext_modules = cythonize(extensions)
)



# subprocess.run(['mv', 'FFT.cpp', './lib/'])