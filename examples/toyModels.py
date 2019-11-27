#!/usr/bin/env python
from __future__ import print_function
#Add paths, we want to be able to run in either root or Run/
import sys,os
#print sys.path
sys.path=["../py","Run", "../Run"]+sys.path

import numpy as np
import time
from scipy.special import ndtri
import math

from numpy.random import RandomState
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import corner 
from PIL import Image
import nestle

from NestleConection import *

rstate = RandomState(0)

if (len(sys.argv)<4):
    print("""Usage: 

python toyModels.py sampler toymodel nlivepoints accuracy

where sampler is one of 

%s 

where

snest: singlenest (ellipsoidal), mnest: multinest, b: bambi (ANN+multinest),
bs: ANN+singlenest

and toymodel is one of 

%s

Example:

Run/toyModels.py egg mnest 50 0.7

"""%("snest, mnest, bambi, sbambi","egg, ring, square, gauss, shells, square"))
    sys.exit(1)

"""
Some toy models for test the nested samplers with and without ANN. 
"""
toymodel, samplername  = sys.argv[1:3]

nlivepoints = int(sys.argv[3])

accuracy = float(sys.argv[4])

dims = 2

path = os.getcwd() + "/toyModelsPlots"

try:
    os.mkdir(path)
except OSError:
    print("Creation of the directory %s failed or already exists" % path)

outputname = toymodel+"_"+samplername+"_"+str(nlivepoints)+'_'+str(accuracy)
#################################################################################

if toymodel =='ring':
	infx, supx, infy, supy, infz, supz = [-5., 5., -5., 5., 0., 10.]
elif toymodel == 'gauss':
	infx, supx, infy, supy, infz, supz = [-5., 5., -5., 5., 0., 1.]
elif toymodel == 'egg':
	infx, supx, infy, supy, infz, supz = [0., 1., 0., 1., -100., 300.]
elif toymodel == 'himmel':
	infx, supx, infy, supy, infz, supz = [-5., 5., -5., 5., 0., 1.]
elif toymodel == 'square':
	infx, supx, infy, supy, infz, supz = [-5., 5., -5., 5., 0., 1.]

infpr, suppr = infx, supx


#########################PRIOR TRANSFORM#########################################
def genericPriorTransform(cube):
	return cube*(suppr-	infpr)+ infpr
#################################################################################
######################eggbox prior and loglikelihood#############################

tmax = 5.0 * np.pi
constant = np.log(1.0 / tmax**2)

def eggLoglike(cube):
	t = 2.0 * tmax * cube - tmax
	a = (2.0 + np.cos(t[0]/2.0)*np.cos(t[1]/2.0))**5.0
	return (2.0 + np.cos(t[0]/2.0)*np.cos(t[1]/2.0))**5.0


################loglikes##############################################

def himmelLoglike(cube):
    return -(cube[0]**2+cube[1]-11)**2.0-(cube[0]+cube[1]**2-7)**2
   
def gaussLoglike(x):
        return -((x[0])**2+(x[1])**2/2.0- 1.0*x[0]*x[1])/2.0

def ringLoglike(x):
        r2=x[0]**2+x[1]**2
        return -(r2-4.0)**2/(2*0.5**2)

def squareLoglike(x):
	sq = 1.
	if abs(x[0])>5 or abs(x[1])>5:
		sq = 0.
	return sq
######################PLOTERS##########################################
def genericPlot2D(res, likeli, infx, supx, infy, supy, time):
	fig = plt.figure(figsize=(8., 8.))
	ax = fig.add_subplot(121, aspect=1)
	xx, yy = np.meshgrid(np.linspace(infx, supx, 50),np.linspace(infy, supy, 50))
	if toymodel != 'square':
		Z = likeli(np.array([xx, yy]))
		ax.contourf(xx, yy, Z, 12, cmap=plt.cm.Blues_r)
		ax.set_xlim(infx, supx)
		ax.set_ylim(infy, supy)
	else:
		ax.scatter(xx,yy)
		ax.set_xlim(infx-1, supx+1)
		ax.set_ylim(infy-1, supy+1)
	ax.set_title('True Log Likelihood')
	
	ax = fig.add_subplot(122, aspect=1)
	#if toymodel == 'egg' or toymodel == 'himmel':
	if toymodel == 'egg':
		ax.scatter(res.samples[:,0], res.samples[:, 1], res.logl, marker='.')
	else:
		ax.scatter(res.samples[:,0], res.samples[:, 1], np.exp(res.logl), marker='.')
	if toymodel == 'square':
		ax.set_xlim(infx-1, supx+1)
		ax.set_ylim(infy-1, supy+1)
	else:
		ax.set_xlim(infx, supx)
		ax.set_ylim(infy, supy)
	ax.set_title('Nested sampling points')

	#fig.suptitle(toymodel +' with '+ samplername + '\nNo of live points:'+\
	#			str(nlivepoints) + '\nAccuracy:'+str(accuracy) +\
	#  			'\nTime = '+str(time)+' seconds = '+str(time/60)+' minutes ')

	fig.suptitle('%s with %s \nNo of live points : %s \nAccuracy : %s \n%.3f seconds = %.3f minutes'\
				%(toymodel,samplername,str(nlivepoints),str(accuracy),time,time/60))

	fig.savefig(path+'/'+outputname+"_2D.png")
	img = Image.open(path+'/'+outputname+"_2D.png")
	img.show()


	fig2 = corner.corner(res.samples, weights=res.weights, bins=500,\
					color = 'r', show_titles=True,\
                    range=[(infx, supx), (infx, supx)])
	fig2.set_size_inches(8., 8.)
	fig2.suptitle(toymodel +' with '+ samplername + '--- No of live points:' + \
				 str(nlivepoints) + '--- Accuracy:'+str(accuracy)+'\n\n')

	fig2.savefig(path+'/'+outputname+"_CORNER.png")
	
	img = Image.open(path+'/'+outputname+"_CORNER.png")
	img.show()


def genericPlot3D(res, loglike, infx, supx, infy, supy, infz, supz, time):
	xx, yy = np.meshgrid(np.linspace(infx, supx, 200),
	                     np.linspace(infy, supy, 200))
	
	Z = np.exp(loglike(np.array([xx, yy])))
	
	fig = plt.figure(figsize=(14., 6.))
	ax = fig.add_subplot(121, projection='3d')
	if toymodel == 'egg' or toymodel == 'himmel':
		ax.scatter(xx,yy,Z,c='k',vmin=0.5,vmax=1,s=2, edgecolor='none',marker='o')
	else:
		ax.plot_surface(xx, yy, Z, rstride=1, cstride=1, linewidth=0, cmap='coolwarm')
	ax.set_xlim3d(infx, supx)
	ax.set_ylim3d(infy, supy)
	ax.set_zlim3d(infz, supz)
	ax.set_zlabel('L')
	ax.set_title('Likelihood evaluated on fine grid')
	if toymodel == 'egg':
		ax.invert_zaxis()
		ax.view_init(elev=40., azim=-35)
	

	ax = fig.add_subplot(122, projection='3d')
	#scatter
	ax.scatter(res.samples[:,0], res.samples[:, 1], np.exp(res.logl),
	           marker='.', c=np.exp(res.logl), linewidths=(0.,), cmap='coolwarm')
	ax.set_xlim(infx, supx)
	ax.set_ylim(infy, supy)
	ax.set_zlim(infz, supz)
	ax.set_zlabel('L')
	ax.set_title('Nested sampling points')
	if toymodel == 'egg':
		ax.invert_zaxis()
		ax.view_init(elev=40., azim=-35)
	
	fig.suptitle(toymodel +' with '+ samplername + '\nNo of live points : '+\
				str(nlivepoints) + '\nAccuracy : '+str(accuracy))
	
	fig.savefig(path+'/'+outputname+"_3D.png")
	img = Image.open(path+'/'+outputname+"_3D.png")
	img.show()

##############################################################################
dispatcher = {'ringLoglike' : ringLoglike, 'eggLoglike' : eggLoglike,\
			'himmelLoglike': himmelLoglike,\
			'gaussLoglike' : gaussLoglike,\
			'genericPriorTransform' : genericPriorTransform,\
			'squareLoglike' : squareLoglike}

def call_func(func):
    try:
        return dispatcher[func]
    except:
        return "Invalid function"

priorTransform = call_func('genericPriorTransform')
likeli = dispatcher[toymodel+'Loglike']
 
##############################################################################
print("The number of live points is %d "%(nlivepoints))

print('---Using NESTLE---')

if samplername in ['mnest','snest']:

	print ('-'*55)
	if samplername == 'mnest':
		print("The sampler used is MULTINEST. Feroz et al (2009)")
		method='multi'
	elif samplername == 'snest':
		print("Mukherjee, Parkinson & Liddle (2006) nested sampling")
		method='single'
		print ('-'*55)
	ti = time.time()
    #IGV: npoints = nsamp or, maybe better, npoints=dims*25
	M = nestle.sample(likeli, priorTransform, ndim=dims, method=method,\
		npoints=nlivepoints, dlogz=accuracy, callback=nestle.print_progress)
	tf = time.time()
	ttime = tf-ti

if samplername in ['bambi','sbambi']:
	from bambi import *

	print("ANN based on pybambi. Graff et al (2012).")
	if samplername =='bambi':
		method = 'mnest'
	elif samplername == 'sbambi':
		method = 'snest'
	ti = time.time()
	M = run_pyBAMBI(likeli, priorTransform, dims,\
                    sampler=method, eff=accuracy,\
                    nlive=nlivepoints, learner='keras')
	tf = time.time()
	ttime = tf-ti

################################################################################
print("\nSummary:\n")
print(M.summary())

file = open(path + '/' + outputname + "_summary.txt",'w')
file.write(M.summary())

file.write("\nElapsed time: %.3f minutes = %.3f seconds"%(ttime/60,ttime))
file.close()
    
print("\nElapsed time: %.3f minutes = %.3f seconds"%(ttime/60,ttime))       

skip = 0
saveChainNestle(M,skip,path +'/'+ outputname + ".txt")
    
################################################################################

genericPlot2D(M, likeli, infx, supx, infy, supy, ttime)
#genericPlot3D(M, likeli, infx, supx, infy, supy, infz, supz, ttime)
	
