from distutils.core import setup, Extension

module1 = Extension('fuzzy',
                    sources = ['fuzzymodule.c', 'fuzzy.c'])

setup (name = 'fuzzy',
       version = '1.0',
       description = 'Fuzzy package in C',
       ext_modules = [module1])
