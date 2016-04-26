#!/usr/bin/env python

from scipy import *
import pylab,sys

## get the latest chain dir
sys.path.append("Pack")
import latest

cdir = 'chains_SimpleMC/'


nlist=['DR11LOWZ','DR11CMASS','DR11LyaAuto','DR11LyaCross','MGS','SixdFGS']
modell= ['owaCDM','waCDM','owCDM','wCDM','oLCDM','LCDM']


if len(sys.argv)>1: 
   data=sys.argv[1:]
   Paper = 'F'
else:
   Paper = 'T'
   dataset = "BBAO+SN+Planck"
   nlist.insert(0,'BetouleSN')
   modell =['WeirdCDM','LCDM_Neff','nuLCDM','DecayFrac','EarlyDE','SlowRDE'] + modell

 
#Print same for other datasets
if Paper == 'F':	
  for d in data:
    if d =="BBAO":
        dataset = "BBAO"		
    elif d == "BBAO+Planck":	
        dataset = "BBAO+Planck" 	
        modell =['StepCDM','SlowRDE','WeirdCDM','nuLCDM','nuoLCDM','nuwCDM','DecayFrac','SlowRDE'] + modell
    elif d == "BBAO+SN+Planck":
        dataset = "BBAO+SN+Planck"  
        nlist.insert(0,"BetouleSN")
        modell= ['PolyCDM','StepCDM','WeirdCDM','LCDM_Neff','nuLCDM','nuoLCDM','nuwCDM','DecayFrac','EarlyDE','SlowRDE'] + modell
    else:
        print 'Only BBAO, BBAO+Planck or BBAO+SN+Planck datasets'	

mnames = { 'LCDM':'$\\Lambda$CDM',
           'oLCDM':'$o\\Lambda$CDM',
	   'wCDM': '$w$CDM',	
           'owCDM': '$ow$CDM',
           'waCDM': '$w_0w_a$CDM', 
           'owaCDM': '$ow_0w_a$CDM', 
           'LadoCDM': 'PolyCDM',
           'LCDMmasslessnu': '($\\nu=0$)CDM',
           'nuoLCDM': '$\\nu$oCDM',
	   'nuwCDM': '$\\nu$wCDM',
	   'nuLCDM': '$\\nu$CDM',  
	   'LCDM_Neff': '$\Delta N_{\\rm eff}$',
           'EarlyDE': 'Early DE',
           'DecayFrac': 'Decay DM Fr',
           'WeirdCDM': 'Tuned Osc',
	   'StepCDM': 'Step',
	   'PolyCDM': 'Poly',
           'SlowRDE': 'Slow-roll DE'}


mdof = {   'LCDM':3,
           'oLCDM':4,
           'wCDM': 4,
           'owCDM':5,
           'waCDM': 5,
           'owaCDM': 6,
           'LadoCDM': 6,
           'LCDMmasslessnu': 4,
	   'LCDM_Neff': 5,
           'nuLCDM': 4,
	   'nuoLCDM': 5,
	   'nuwCDM': 5,	
           'EarlyDE': 4,
           'DecayFrac': 5,
           'WeirdCDM': 6,
	   'StepCDM': 7,
	   'PolyCDM': 6,
           'SlowRDE': 4}

colors= { 'DR11LOWZ':'red',
          'DR11CMASS':'blue',
          'DR11LyaAuto':'green',
          'DR11LyaCross':'yellow',
          'MGS':'orange', 
          'SixdFGS':'magenta',
          'SPlanck':'black',
          'BetouleSN':'cyan'}

nnames= { 'DR11LOWZ':'LOWZ',
          'DR11CMASS':'CMASS',
          'DR11LyaAuto':'Ly$\\alpha$ Auto',
          'DR11LyaCross':'Ly$\\alpha$ cross',
	  'MGS':'MGS',
          'SixdFGS':'6dFGS',
          'SPlanck':'Planck',
          'BetouleSN':'SN - 30'}

defdof= { 'DR11LOWZ':1,
          'DR11CMASS':2,
          'DR11LyaAuto':2,
          'DR11LyaCross':2,
	  'MGS':1,
          'SixdFGS': 1,
          'SPlanck': 2,
          'BetouleSN': 1}


patches={}
position = 18

class plot_bar:
    def __init__(self):
        self.cy=1
        self.ys=[]
        self.names=[]
	self.patches={}

        
    def add_bar(self,fname, cname, model):
	chiT, dof =0, 0
        fname=cdir+fname 
        pnames=open(fname+'.paramnames').readlines()
	if 'Neff' in fname:
	    loglsx=open(fname+'.maxlike').readlines()
	    logls2=loglsx
	else:
          for i in range(3):
            loglsx=open(fname+'_'+str(i+1)+'.maxlike').readlines()
            if (i==0):
                logls2=loglsx
            else:
                if float(loglsx[0].split()[1])<float(logls2[0].split()[1]):
                    logls2=loglsx
        logls=logls2[0].split(' ')[2:]
        chi2d={}
	print ' '
	print '++++'+ model
        for pname, logl in zip(pnames,logls):
            if "_like" in pname:
                ppname=pname.split(' ')
		if 'SPlanck' in ppname[0]:
                   chi2=0 
		else:
		   if 'Neff' in fname:
	  	       chi2= float(logl)
		       if 'Betoule' in ppname[0]:
			   chi2=chi2 -692
                   else:
                       chi2=-2*float(logl)
		       if 'Betoule' in ppname[0]:   
                           chi2=chi2 -30 	      
      		chiT+=chi2
                xname=ppname[0].replace('_like','')
		print xname, chi2
                chi2d[xname]=chi2
	print 'Min_chi2 = ', chiT+30

	
	param=mdof[model]
	for xname in nlist:
		dof +=defdof[xname]
        left=0
        for xname in nlist:
            chi2=chi2d[xname]
            color=colors[xname]
            PP=pylab.barh(self.cy-0.25,chi2,left=left,height=0.5,color=color, linewidth=0)
            self.patches[xname]=PP[0]
            left+=chi2
	   
 	    if "SN" in dataset:				
	       pylab.text(position, self.cy-0.25, r' %.2f / %s'%(chiT +30, dof-param +30)  , fontsize=15)
	    else:
	       pylab.text(position, self.cy-0.25, r' %.2f / %s'%(chiT, dof-param)  , fontsize=15)		
        self.ys.append(self.cy)
        self.cy+=1
        self.names.append(cname)



    def add_def_bar(self):
        left=0
        for xname in nlist:
            color=colors[xname]
            chi2=defdof[xname]
            PP=pylab.barh(self.cy-0.25,chi2,left=left,height=0.5,color=color, linewidth=0)
            self.patches[xname]=PP[0]
            left+=chi2
        self.ys.append(self.cy)
        self.cy+=2
        self.names.append('D.o.F.')

    def close(self):
        print self.ys
        print self.names
        pylab.yticks(self.ys,self.names)
        pylab.xlabel('$\chi^2$')
        px,py=[],[]
        for col in nlist:
            px.append(self.patches[col])
            py.append(nnames[col])
        pylab.legend(px,py,loc='lower right', bbox_to_anchor=(1.12, 0.6),fancybox=True, shadow=True)


params1 = {'backend': 'pdf',
               'axes.labelsize': 28,
               'text.fontsize': 28,
               'xtick.labelsize': 22,
               'ytick.labelsize': 22,
               'legend.draw_frame': False,
               'legend.fontsize': 16,
               'lines.markersize': 6,
               'font.size': 20,
               'text.usetex': True}#
pylab.rcParams.update(params1)

pylab.figure(figsize=(16,9))
P=plot_bar()

P.add_def_bar()
for model in modell:
    nn=mnames[model]
    if 'Neff' in model:
	P.add_bar(model+'_PL_BBAO_JLA',nn, model)	
    else: 
        P.add_bar(model+'_phy_'+ dataset,nn, model)	

P.close()
pylab.text(position, 1.5, '$\chi^2$ / d.o.f'  , fontsize=20)
pylab.xlim(0,position+2.5)
if Paper == 'F':
 pylab.savefig("Chisq_"+dataset+".pdf")
else: 
 pylab.savefig("Chisq.pdf")
pylab.show()

