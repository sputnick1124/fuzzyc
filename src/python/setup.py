from distutils.core import setup, Extension
import sys

py_version = 'PYTHON' + str(sys.version_info.major)
print(py_version)

module1 = Extension('fuzzy',
                    include_dirs = ['.'],
                    define_macros = [(py_version, '1')],
                    libraries = ['fuzzy'],
                    library_dirs = ['.'],
                    sources = ['fuzzymodule.c'])

setup (name = 'FuzzyModule',
       version = '0.1',
       ext_modules = [module1])
