from simplemc.DriverMC import DriverMC

sampler = 'nested'
a = DriverMC(analyzername=sampler ,model="LCDM", datasets="HD")
a.executer(neuralNetwork=True)
#a.nestedRunner(neuralNetwork=True)
a.postprocess()