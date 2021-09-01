from simplemc.DriverMC import DriverMC
from os import remove

"read all setting from .ini file"
inifile = "baseConfigGUI.ini"

analyzer = DriverMC(iniFile=inifile)
analyzer.executer()
#analyzer.postprocess()

remove("baseConfigGUI.ini")


""" useful for short tests,
    or when just a few settings customize it from here"""
#analyzer = DriverMC(analyzername="nested", model="LCDM", datasets="HD")
#analyzer.executer(nlivepoints=10)
#analyzer.postprocess()
