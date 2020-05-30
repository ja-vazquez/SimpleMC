#!/usr/bin/env python
# encoding: utf-8

import sys
import os
from setuptools import find_packages
os.path = ["simplemc"] + os.path
try:
    from setuptools import setup, Extension
    setup, Extension
except ImportError:
    from distutils.core import setup, Extension
    setup, Extension

from setuptools.command.test import test as TestCommand

class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True
    def run_tests(self):
        #import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.test_args)
        sys.exit(errno)


desc = open("README.rst").read()
required = ["numpy", "emcee", "scipy", "nestle", "keras", "tensorflow", "matplotlib", 
            "corner", "getdist"]
test_requires = ["mock"]

PACKAGE_PATH = os.path.abspath(os.path.join(__file__, os.pardir))

setup(
    name="SimpleMC",
    version='2.0.0',
    author='JA Vazquez, I Gomez-Vargas, A Slosar',
    author_email="jvazquez@icf.unam.mx",
    url="https://github.com/ja-vazquez/SimpleMC",
    license="GPLv3",
    # packages=find_packages(PACKAGE_PATH, "Tests"),
    description="Cosmological parameter estimation with SimpleMC",
    long_description=desc,
    install_requires=required,
    test_requires=test_requires,
    package_data={"": ["LICENSE"],
                  'SimpleMC': ['data/*.dat']},
    include_package_data=True,
    keywords=["SimpleMC",
              "parameter estimation",
              "cosmology",
              "MCMC"],
    cmdclass = {'test': PyTest},
    classifiers=[
        "Development Status :: 2 - Beta",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
    ],
)