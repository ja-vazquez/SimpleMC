#!/usr/bin/env python
# encoding: utf-8

from distutils.core import setup

with open("README.md", "r") as fh:
    desc = fh.read()

required = ["numpy", "scipy", "sklearn", "matplotlib", "numdifftools", "mpi4py"]

setup(
    name="simplemc",
    version='0.9.8',
    author='JA Vazquez, I Gomez-Vargas, A Slosar',
    author_email="jvazquez@icf.unam.mx",
    url="https://github.com/igomezv/simplemc_tests",
    license="GPLv3",
    description="Cosmological parameter estimation with SimpleMC",
    long_description=desc,
    setup_requires=['setuptools_scm'],
    install_requires=required,
    include_package_data=True,
    package_data={'simplemc': ['data/*.txt']},
    keywords=["SimpleMC",
              "parameter estimation",
              "machine learning",
              "dark energy",
              "cosmology",
              "MCMC"],
    classifiers=[
        "Development Status :: 0.9.8 - Beta",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        'License :: OSI Approved :: GPLv3 License',
        "Programming Language :: Python",
        'Natural Language :: English',
    ],
    python_requires='>=3.6',
)
