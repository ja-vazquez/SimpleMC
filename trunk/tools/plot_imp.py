#!/usr/bin/env python

from scipy import *
from cosmich import *
from ChainIterator import *
from RunBase import *
import pylab, sys
import string
import math as N

#params1 = {'backend': 'pdf',
#               'axes.labelsize': 20,
#               'text.fontsize': 18,
#               'xtick.labelsize': 20,
#               'ytick.labelsize': 20,
#               'legend.draw_frame': False,
#               'legend.fontsize': 16,
#               'lines.markersize': 6,
#               'font.size': 20,
#               'text.usetex': True}#
#pylab.rcParams.update(params1)


#2D-plotting is valid for only one couple of parameters, and may use several models/datasets
#1D-plotting is valid for only one model, and may use several parameters/datasets

#-----------
dire = 'chains_SimpleMC/'
name_fig  = 'Plot'
datasetl  = ['BBAO']

#------------------
Plot_1D   = 'True'
model_1D  = 'LCDM_phy'
params_1D = ['h']
NBins_1D  = 40

xrange    = 'False'
xmin_1, xmax_1 = 0, 3
xmin_2, xmax_2 = 0, 5

#-----------------------
Plot_2D   = 'False'
model_2D  = ['LCDM_phy']
param_2D  = ['Om','h']
NBins_2D  = 40

xrange_2D = 'False'
xmin, xmax = 0.2, 0.4
ymin, ymax = 0.6, 0.75

#---------------------------------------------------------------------------------


def colour(x):
    if x==1: return 'black'
    if x==2: return 'blue'
    if x==3: return 'red'
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
        return cosmodata


if 'True' in Plot_2D:
  a=0
  for model in model_2D:    
     for datasets in datasetl:
        a+=1             

	C=cosmochain(dire+ model+'_'+datasets)
        C.Plot2D(param_2D[0], param_2D[1], filled=colour(a), label=cosmodata(datasets), N=NBins_2D)

  if 'True' in xrange_2D:   
     pylab.xlim(xmin,xmax)
     pylab.ylim(ymin,ymax)

  leg=pylab.legend()
  leg.draw_frame(False)		# No box & colour legend
  color_legend(leg)

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
	C=cosmochain(dire+ model_1D+'_'+datasets,'auto')	

        C=ChainIterator('chains_140411_155701','LCDM','phy','BBAO+CMBP')
        grlist=[]
	wlist=[]
	for i in range(0,C.N,1000):
	    T=C.theory(i)
	    grlist.append(T.growth(10.0))
	    wlist.append(C.weight(i))

	grlist=array(grlist)
	wlist=array(wlist)
	mn=(grlist*wlist).sum()/wlist.sum()
	er=sqrt((grlist**2*wlist).sum()/wlist.sum()-mn**2)
	print "growth z=10/z=0 = ",mn,"+/-",er



	a+=1
        xx,yy=C.GetHisto(params,NormPeak=True,nbins=NBins_1D)	
        pylab.plot(xx, yy, colour(a), label=cosmodata(datasets))

        if 'True' in xrange:
            xmin="xmin_"+str(b)
            xmax="xmax_"+str(b) 
            pylab.xlim(eval(xmin),eval(xmax))

        pylab.xlabel(C.latexname(str(params)))
        pylab.ylabel('prob.')

   leg=pylab.legend(loc='upper right')
   leg.draw_frame(False)
   color_legend(leg)
   pylab.savefig(name_fig+'_1D.pdf')
   pylab.show()
   
else:
   print 'Nothing else to do'


