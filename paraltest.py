# mpirun -np 2 python3 paraltest.py

from simplemc.DriverMC import DriverMC
from mpi4py import MPI

name = MPI.Get_processor_name()
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

datasets = ["SN", "BBAO+HD"]
models = ["LCDM", "owaCDM"]

a = DriverMC(model=models[rank], datasets=datasets[rank], analyzername="nested")
a.executer(nlivepoints=30, accuracy=0.8)
a.postprocess(addtxt=["\nUsed the {} processor of {} with MPI\n".format(rank, size)])

# print (sys.argv[1], sys.argv[2])
#comd ="Run/driver_old.py phy %s %s %s 100 80000"%(sys.argv[1], sys.argv[2], rank+1)
# os.system(comd)


