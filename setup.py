from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

fftw_include_dir = '/usr/local/include'
lib_dir = '/usr/local/lib'

extensions = [
    Extension(name="slant_stack", 
                sources=["slant_stack.pyx"],
                language="c++",
                include_dirs=[fftw_include_dir],
                library_dirs=[lib_dir],
                libraries=["fftw3", "omp"],
                extra_compile_args=["-O3"],
                extra_link_args=["-O3"]),
]

setup(
    name="slant_stack",
    ext_modules=cythonize(extensions),
)



# subprocess.run(['mv', 'FFT.cpp', './lib/'])