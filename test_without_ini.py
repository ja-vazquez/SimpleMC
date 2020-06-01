from simplemc.DriverMC import DriverMC

sampler = "emcee"

analyzer = DriverMC(analyzername=sampler, model="LCDM", datasets="HD")
analyzer.executer(nsamp=5000, walkers=100)
analyzer.postprocess(addtxt=["Autoccorelation time failed"])