from setuptools import setup

setup(
    name="corisco",
    version="0.1dev",
    packages=["corisco"],

    tests_require=['Attest'],
    test_loader='attest:Loader',
    test_suite='tests.collection',
)
