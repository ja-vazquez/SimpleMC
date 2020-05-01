#!/usr/bin/env python

from DriverMC import DriverMC

import pathlib

fileConfig = '%s/baseConfig.ini'%pathlib.Path().absolute()

D = DriverMC(fileConfig)

#D.plotter()

