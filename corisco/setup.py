from setuptools import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    name="corisco",
    version="0.1dev",
    packages=["corisco"],
    entry_points={
        "console_scripts": [
            "corisco = corisco.tools.corisco_estimate_orientation:main",
        ]
    },

    tests_require=['Attest'],
    test_loader='attest:Loader',
    test_suite='tests.collection',

    install_requires=[
        'filtersqp'
    ],

    cmdclass = {'build_ext': build_ext},
    ext_modules=[Extension("corisco.corisco_aux",
                           ["corisco/corisco_aux.pyx"],
                           include_dirs=['/usr/lib/python2.7/site-packages/numpy/core/include/'],
                           extra_objects=['corisco/corisco_aux_external.o'],
                           libraries=['m','blas'])],

)

