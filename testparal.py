from simplemc.run_paral import *

list_models = ['LCDM', 'owaCDM']

obs = [['BBAO'], ['HD+BBAO']]

#
# from mpi4py import MPI
#
# comm = MPI.COMM_WORLD

print("Hello! I'm rank " + str(comm.rank))

analyzer = DriverMC(model=model[comm.rank], datasets=data[comm.rank], analyzername=sampler).executer()

#outputxt, params = analyzer.executer(nsamp=100, skip=0)
#analyzer.executer()

