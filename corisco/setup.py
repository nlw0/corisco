from setuptools import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [
    Extension('corisco.corisco_aux',
              ['corisco/corisco_aux.pyx', 'corisco/corisco_aux_external.c'],
              include_dirs=['/usr/lib/python2.7/site-packages/numpy/core/include/'],
              # extra_objects=['corisco_aux_external.o'],
              extra_compile_args=['-O3', '-fPIC', '-msse2', '-mfpmath=sse'],
              libraries=['m', 'blas'])
    ]

setup(
    name='corisco',
    version='0.1dev',
    packages=['corisco'],
    entry_points={
        'console_scripts': [
            'corisco = corisco.tools.corisco_estimate_orientation:main',
        ]
    },

    tests_require=['Attest'],
    test_loader='attest:Loader',
    test_suite='tests.collection',

    install_requires=[
        'filtersqp'
    ],

    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
)

