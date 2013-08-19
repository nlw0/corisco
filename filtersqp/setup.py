import distutils
from setuptools import setup

import numpy

setup(
    name="filtersqp",
    version="0.1dev",
    packages=["filtersqp"],
    entry_points={
        "console_scripts": [
            "sqp_plot_root = filtersqp.tools.sqp_plot_root:main",
            "sqp_plot_filter = filtersqp.tools.sqp_plot_filter:main",
        ]
    },

    tests_require=['Attest'],
    test_loader='attest:Loader',
    test_suite='tests.collection',

    ext_modules=[
        distutils.extension.Extension(
            "filtersqp.trust_rootfind",
            sources=["filtersqp/trust_rootfind_module.cc"],
        ),
    ],

    install_requires=[
        'numpy'
    ],

    include_dirs=[numpy.get_include()],
)
