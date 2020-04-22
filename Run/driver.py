#!/usr/bin/env python
from DriverMC import DriverMC
#from PlotterMC import PlotterMC
import pathlib

fileConfig = '%s/baseConfig.ini'%pathlib.Path().absolute()

D = DriverMC(fileConfig)

#D.plotter()
