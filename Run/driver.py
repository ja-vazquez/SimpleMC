#!/usr/bin/env python

from DriverMC import DriverMC
<<<<<<< HEAD

fileConfig = "baseConfig.ini"
=======
#from PlotterMC import PlotterMC
import pathlib

fileConfig = '%s/baseConfig.ini'%pathlib.Path().absolute()
>>>>>>> 33667434ec3900d2c4380db33a366dd782e1a47b

D = DriverMC(fileConfig)

#D.plotter()

