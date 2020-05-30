from simplemc.DriverMC import DriverMC
from simplemc.tools.Simple_Plots import Simple_plots

# sampler = 'mcmc'
# optimizer = 'genetic'
# model = "owaCDM"
# data = 'HD+BBAO'
# inifile = "baseConfig.ini"
# chainsdir = "/home/isidro/SimpleMC/SimpleMC_chains/"
# root = "LCDM_phy_HD_mcmc"
# inifile='baseConfig.ini'

#analyzer = DriverMC(iniFile=inifile)
analyzer = DriverMC(analyzername="nested", model="LCDM", datasets="HD")
analyzer.executer(nlivepoints=10)
#analyzer = DriverMC(model=model, datasets=data, analyzername=sampler, priortype='u')

#outputxt, params = analyzer.executer(nsamp=100, skip=0)
# analyzer.executer(n_individuals=100, n_generations=100, mut_prob=0.4)

#analyzer.postprocess(summary=True, addtxt=['Probando SimpleMC'])

# plt = splt.Plots1D(chainsdir, root)
