from simplemc.DriverMC import DriverMC

inifile = "tallerConfig.ini"

analyzer = DriverMC(iniFile=inifile)
analyzer.executer()
analyzer.postprocess()

fig = analyzer.plot(show=True)
txt = "Bienvenid@ al\ntaller de cosmología\ncon SimpleMC 2.0"
fig.simpleCorner(label=txt)

