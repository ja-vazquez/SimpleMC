#!/usr/bin/env python

from scipy import *
import pylab,sys


## get the latest chain dir
sys.path.append("Pack")
import latest

#cdir="/astro/u/anze/work/April/chains/"
cdir="/gpfs01/astro/workarea/jvazquez/chains/"

nlist=['BetouleSN','DR11LOWZ','DR11CMASS','DR11LyaAuto','DR11LyaCross','MGS','SixdFGS']

modell= ['LCDM', 'oLCDM', 'wCDM', 'owCDM', 'waCDM', 'owaCDM', 'nuLCDM_cosmo']

mnames = { 'LCDM':'$\\Lambda$CDM',
           'oLCDM':'$o\\Lambda$CDM',
	   'wCDM': '$w$CDM',	
           'owCDM': '$ow$CDM',
           'waCDM': '$w_0w_a$CDM', 
           'owaCDM': '$ow_0w_a$CDM', 
           'LadoCDM': 'PolyCDM',
           'LCDMmasslessnu': '($\\nu=0$)CDM',
           'nuLCDM_cosmo': '$\\nu$CDM'  }



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
          'SPlanck': 3,
          'BetouleSN': 3}


patches={}

class plot_bar:
    def __init__(self):
        self.cy=1
        self.ys=[]
        self.names=[]
        self.patches={}
        
    def add_bar(self,fname, cname):
        fname=cdir+fname 
        pnames=open(fname+'.paramnames').readlines()
        for i in range(3):
            loglsx=open(fname+'_'+str(i+1)+'.maxlike').readlines()
            if (i==0):
                logls2=loglsx
            else:
                if float(loglsx[0].split()[1])<float(logls2[0].split()[1]):
                    logls2=loglsx
        logls=logls2[0].split(' ')[2:]
        chi2d={}
        for pname, logl in zip(pnames,logls):
            if "_like" in pname:
                ppname=pname.split(' ')
                if 'Betoule' in ppname[0]:
                   chi2=-2*float(logl)-30
                   print "bchi2=",chi2
                else:
                   chi2=-2*float(logl)      
                xname=ppname[0].replace('_like','')
                chi2d[xname]=chi2
        left=0
        for xname in nlist:
            chi2=chi2d[xname] #/defdof[xname]
            color=colors[xname]
            PP=pylab.barh(self.cy-0.25,chi2,left=left,height=0.5,color=color, linewidth=0)
            self.patches[xname]=PP[0]
            left+=chi2
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
        pylab.legend(px,py,loc='lower right')


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

pylab.figure(figsize=(15,9))
P=plot_bar()

P.add_def_bar()
for model in modell:
    nn=mnames[model]
    P.add_bar(model+'_phy_BBAO+SN+Planck',nn)
P.close()
pylab.xlim(0,21)
pylab.savefig("Chisq.pdf")
pylab.show()

