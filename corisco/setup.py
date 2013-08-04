from setuptools import setup

setup(
    name="corisco",
    version="0.1dev",
    packages=["corisco"],
    entry_points={
        "console_scripts": [
            "corisco_xxx = tools.corisco_xxx:main",
        ]
    },

    tests_require=['Attest'],
    test_loader='attest:Loader',
    test_suite='tests.collection',

    install_requires=[
        'filtersqp'
    ],
)
