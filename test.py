from simplemc.DriverMC import DriverMC

"read all setting from .ini file"
inifile = "baseConfig.ini"

analyzer = DriverMC(iniFile=inifile)
analyzer.executer()
analyzer.postprocess()


"when just a few settings customize them from here"
#analyzer = DriverMC(analyzername="nested", model="LCDM", datasets="HD")
#analyzer.executer(nlivepoints=10)
#analyzer.postprocess()