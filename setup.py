from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

ext_modules=[Extension("uks_aello",["uks_aello.pyx"],libraries=["m"],extra_compile_args=["-ffast-math"])]
setup(name='uks_aello',ext_modules = cythonize('uks_aello.pyx',language_level=3))

