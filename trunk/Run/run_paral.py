



from mpi4py import MPI
import os, sys

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

print (sys.argv[1], sys.argv[2])
comd ="Run/driver.py phy %s %s %s 50 50000"%(sys.argv[1], sys.argv[2], rank+1)
os.system(comd)

