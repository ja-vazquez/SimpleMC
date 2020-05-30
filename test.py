from simplemc.DriverMC import DriverMC

# inifile = "baseConfig.ini"

#analyzer = DriverMC(iniFile=inifile)
analyzer = DriverMC(analyzername="nested", model="LCDM", datasets="HD")
analyzer.executer(nlivepoints=10)