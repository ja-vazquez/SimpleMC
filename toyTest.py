from simplemc.tools.ToyModels import ToyModels as TM
from simplemc.analyzers import NestedSampler
from simplemc.analyzers.pybambi.bambi import bambi

import numpy as np
import matplotlib.pyplot as plt
import time

import multiprocessing as mp

"""
This script calls toy distributions from the ToyModel class and make a sampling 
for these models through dynesty with and without a neural network (based on pybambi).
"""
np.random.seed(0)

# ##### SETTINGS ###########
# modelname can be {'egg', 'himmel', 'ring', 'square', 'gaussian'}
modelname = 'egg'
show_plots = True  # choose False if you are in a server
tm = TM(model=modelname)  # create a ToyModel instance
priorTransform = tm.priorTransform  # uniform prior transform
loglike = tm.loglike
dims = 2
nlive = 500

# Parallel options
nworkers = 2
pool = mp.Pool(nworkers)

# ###### FIRST SAMPLING WITH ONLY DYNESTY
sampler1 = NestedSampler(loglike, priorTransform, ndim=dims,
                        bound='multi', sample='unif', nlive=nlive,
                        pool=pool, queue_size=nworkers,
                        use_pool={'loglikelihood': False})

ti = time.time()
sampler1.run_nested(dlogz=0.01, outputname=modelname+"_dynesty")
resultnested = sampler1.results
tfnested = time.time() - ti

# ### CONFIGURE pybambi neural network
thumper = bambi(loglike, dims, learner='keras', proxy_tolerance=100.0, nlive=nlive, epochs=100,
                updInt=200, dlogz_start=10)
neural_logLike = thumper.loglikelihood
dumper = thumper.dumper

# ###### SECOND SAMPLING WITH DYNESTY + NEURAL NETWORK
sampler2 = NestedSampler(neural_logLike, priorTransform, ndim=dims,
                        bound='multi', sample='unif', nlive=nlive,
                        pool=pool, queue_size=nworkers,
                        use_pool={'loglikelihood': False})

print("\nNext sampling:")
ti = time.time()
sampler2.run_nested(dlogz=0.01, outputname=modelname+"_bambi", dumper=dumper, netError=0.1)
resultbambi = sampler2.results
tfbambi = time.time() - ti

# ###### PRINT SUMMARY OF BOTH SAMPLING PROCESSES

print("\n\nDynesty:")
resultnested.summary()

print("\n\nDynesty + ANN :")
resultbambi.summary()
print("\nTime dynesty: {:.4f} min | Time dynesty+ANN: {:.4f} min".format(tfnested/60, tfbambi/60 ))

# ### Plot results if you aren't in a server
if show_plots:
    nestdata = np.loadtxt(modelname+'_dynesty_1.txt', usecols=(2, 3))
    bambidata = np.loadtxt(modelname+'_bambi_1.txt', usecols=(2, 3))
    znest = np.zeros(len(nestdata))

    for i, row in enumerate(nestdata):
        znest[i] = tm.loglike(row)

    zbambi = np.zeros(len(bambidata))

    for i, row in enumerate(bambidata):
        zbambi[i] = tm.loglike(row)
    fig = plt.figure()

    # 3D plots
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(nestdata[:, 0], nestdata[:, 1],  znest, c='red', alpha=0.5)
    ax.scatter(bambidata[:, 0], bambidata[:, 1], zbambi, c='green', alpha=0.5)
    plt.legend(["dynesty", "dynesty + neural net"],  loc="upper right")

    plt.savefig(modelname+"_bambi3D.png")
    plt.show()

    # 2D plots
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.scatter(nestdata[:, 0], nestdata[:, 1], c='red', alpha=0.5)
    plt.scatter(bambidata[:, 0], bambidata[:, 1], c='green', alpha=0.5)
    plt.legend(["dynesty", "dynesty + neural net"],  loc="upper right")
    ax.set_aspect('equal', adjustable='box')
    plt.savefig(modelname+"_bambi2D.png")
    plt.show()

