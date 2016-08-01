from distutils.core import setup, Extension
import sys

if sys.version_info[0] == 3:
	print("Using Python 3")
	macros = [('PYTHON3', '1')]
else:
	macros = [('PYTHON2', '1')]

module1 = Extension('fuzzy',
                    sources = ['src/fuzzymodule.c', 'src/fuzzy.c'],
					include_dirs = ['include'],
					define_macros = macros,
					extra_compile_args=["-O4"])

setup (name = 'fuzzy',
       version = '1.0',
       description = 'Fuzzy package in C',
       ext_modules = [module1])

print(macros[0][0])
