#!/usr/bin/env python

from scipy import *
from cosmich import *
from ChainIterator import *
from RunBase import *
import pylab, sys
import string
import math as N

params1 = {'backend': 'pdf',
               'axes.labelsize': 20,
               'text.fontsize': 18,
               'xtick.labelsize': 20,
               'ytick.labelsize': 20,
               'legend.draw_frame': False,
               'legend.fontsize': 16,
               'lines.markersize': 6,
               'font.size': 20,
               'text.usetex': True}#
pylab.rcParams.update(params1)


#2D-plotting is valid for only one couple of parameters, and may use several models/datasets
#1D-plotting is valid for only one model, and may use several parameters/datasets

#-----------
dire = 'chains/DR12/' #_SimpleMC/decay/test/'
name_fig  = 'Plot'
datasetl = ['CombBAOzb1']
legend = ['BAO_all_phy']
#datasetl  = ['Planck+BBAO_GS500','*BBAO+Planck']
#legend = ['GMSampler','MCMC']

#------------------
Plot_1D   = 'True'
model_1D  = 'LCDM_phy'
params_1D = ['Om','h'] #, 'A', 'B']
NBins_1D  = 30

xrange    = 'False'
xmin_1, xmax_1 = 0.6, 0.75
xmin_2, xmax_2 = 0, 5

#-----------------------
Plot_2D   = 'False'
model_2D  = ['Quint2_model_phy']
param_2D  = ['lam','B']
NBins_2D  = 15

xrange_2D = 'False'
xmin, xmax = -1.5, -0.5
ymin, ymax = -0.025, 0.025

#---------------------------------------------------------------------------------


def colour(x):
    if x==1: return 'red'
    if x==2: return 'blue'
    if x==3: return 'black'
    if x==4: return 'magenta'
    if x==5: return 'cyan'
    if x==6: return 'orange'
    if x==7: return 'green'
    if x==8: return 'yellow'
    if x>8:  print("Increased colouring") 


def color_legend(leg):
        # """Color legend texts based on color of corresponding lines"""
        for line, txt in zip(leg.get_lines(), leg.get_texts()):
                txt.set_color(line.get_color())


def cosmodata(datasets):
        cosmodata=''
        if 'BBAO' in datasets:
                if '+BBAO' in datasets:
                   cosmodata=cosmodata+'+'
                cosmodata=cosmodata+'BAO'
        if 'GBAO' in datasets:
                if '+GBAO' in datasets:
                   cosmodata=cosmodata+'+'
                cosmodata=cosmodata+'Gal BAO'
        if 'LBAO' in datasets:
                if '+LBAO' in datasets:
                   cosmodata=cosmodata+'+'
                cosmodata=cosmodata+'Lya BAO'
        if 'SN' in datasets:
                if '+SN' in datasets:
                   cosmodata=cosmodata+'+'
                cosmodata=cosmodata+'SNe'
        if "Planck"  in datasets:
                if '+Planck' in datasets:
                   cosmodata=cosmodata+'+'
                cosmodata=cosmodata+'Planck'           
        if '6dFGS' in datasets:
                if '+6dFGS' in datasets:
                   cosmodata=cosmodata+'+'
                cosmodata=cosmodata+'6dFGS'
	if '*' in datasets:
		cosmodata='MCMC'
        return cosmodata


def get_filename(dire,model,dataset):
 	if '*' in dataset:
	    dire = 'chains_SimpleMC/'
	    dataset='BBAO+Planck'		
	return dire+ model+'_'+dataset	


if 'True' in Plot_2D:
  a=0
  for model in model_2D:    
     for datasets in datasetl:
        a+=1             
	
	if legend:
	   labels= legend[a-1]
	else:
	   labels =cosmodata(datasets)	

	fname=get_filename(dire,model,datasets)
	C=cosmochain(fname) #, nums=[1])
        C.Plot2D(param_2D[0], param_2D[1], filled=colour(a), label=labels, N=NBins_2D)
	#C.Plot2D(param_2D[0], param_2D[1],lims=[0.2,0.4,0.6,0.8], filled=colour(a), label=cosmodata(datasets), N=NBins_2D)	

  if 'True' in xrange_2D:   
     pylab.xlim(xmin,xmax)
     pylab.ylim(ymin,ymax)

  leg=pylab.legend(loc='upper left')
  leg.draw_frame(False)		# No box & colour legend
  color_legend(leg)

  pylab.tight_layout()
  pylab.xlabel(C.latexname(param_2D[0]))
  pylab.ylabel(C.latexname(param_2D[1]))
  pylab.savefig(name_fig+'_2D.pdf')
  pylab.show()



if 'True' in Plot_1D:
   a=0
   for datasets in datasetl:
      b=0
      for params in params_1D:
        b+=1
	pylab.subplot(1,len(params_1D),b)

	#fname=self.get_filename(dire,model,extra,dataset)
	C=cosmochain(dire+ model_1D+'_'+datasets,'auto')	

	if(True):                  #Compute functions= ratio of Fig 2 at z=0.57/0
	  dire2 = 'chains/DR12/'
          D=ChainIterator(dire2, 'LCDM','phy', 'CombBAOzb3')
          grlist=[]
	  wlist=[]
	  
	  for i in range(0,D.N,1000):
	    T=D.theory(i)
	    #grlist.append(T.Da_z(0.57)/T.Hinv_z(0)/(1.0*N.log(1.57)))	
	    #grlist.append(T.DaOverrd(0.51)) 
	    grlist.append(T.HIOverrd(0.61)) 
	    wlist.append(D.weight(i))

	  grlist=array(grlist)
	  wlist=array(wlist)
	  mn=(grlist*wlist).sum()/wlist.sum()
	  er=sqrt((grlist**2*wlist).sum()/wlist.sum()-mn**2)
	  #print "H(0.0)/(1.0) \over c*ln(1.57)/D_m(0.57)  = ",mn,"+/-",er
	  print "z = 0.38  = ",mn,"+/-",er	  

	a+=1
        xx,yy=C.GetHisto(params,NormPeak=True,nbins=NBins_1D)	
        pylab.plot(xx, yy, colour(a), label=cosmodata(datasets))

        if 'True' in xrange:
            xmin="xmin_"+str(b)
            xmax="xmax_"+str(b) 
            pylab.xlim(eval(xmin),eval(xmax))

        pylab.xlabel(C.latexname(str(params)))
        pylab.ylabel('prob.')

   #leg=pylab.legend(loc='upper right')
   #leg.draw_frame(False)
   #color_legend(leg)
   pylab.savefig(name_fig+'_1D.pdf')
   pylab.show()
   
else:
   print 'Nothing else to do'


