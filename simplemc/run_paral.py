#mpirun -np 2 python3 Run/run_paral.py model data
from mpi4py import MPI
import os, sys

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# print (sys.argv[1], sys.argv[2])
# comd ="Run/driver_old.py phy %s %s %s 100 80000"%(sys.argv[1], sys.argv[2], rank+1)
# os.system(comd)

