#!/usr/bin/env python
# encoding: utf-8

from setuptools import setup


with open("README.md", "r") as fh:
    desc = fh.read()

required = ["numpy", "scipy", "sklearn", "matplotlib", "numdifftools", "mpi4py"]

setup(
    name="simplemc",
    version='0.9.8',
    author='JA Vazquez, I Gomez-Vargas, A Slosar',
    author_email="jvazquez@icf.unam.mx",
    url="https://github.com/ja-vazquez/SimpleMC",
    license="GPLv3",
    description="Cosmological parameter estimation with SimpleMC",
    long_description=desc,
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
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Natural Language :: English',
    ],

)
