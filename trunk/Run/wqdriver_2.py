#!/usr/bin/env python

## this scripts runs everything on the bnl cluster
import os, sys

datasetl=['BBAO','GBAO','GBAO_no6dF','LBAO','LBAO+CMBP','LBAO+CMBW','GBAO+CMBP','GBAO+CMBW',
         'BBAO+CMBP','BBAO+CMBW','BBAO+SN','BBAO+SN+RiessH0','SN','SN+CMBW','SN+CMBP',
         'BBAO+SN+CMBP','BBAO+SN+CMBW', 'GBAO+SN+cCMBP']
## removed JordiCDM for the time being
modell= ['LCDM','LCDMmasslessnu','nuLCDM','oLCDM','wCDM','waCDM', 'owaCDM', 'owCDM', 'PolyCDM', 'TLight']

def wqsubmit(prefact,model,datasets):
    comd = "Run/driver.py %s %s_2 %s "%(prefact,model,datasets)
    print comd
    nm=prefact+"_"+model+"_2"+"_"+datasets
    wqcom='wq sub -r job_name:%s -c "%s"'%(nm,comd)
    os.system("nohup %s >%s 2>%s &"%(wqcom,'chains/logs/'+nm+'.log','chains/logs/'+nm+'.err'))


if len(sys.argv)>1:
    if (sys.argv[1]=='hladder'):
        wqsubmit('phy','PolyCDM','GBAO+SN+cCMBP')
        wqsubmit('phy','PolyCDM','GBAO+SNx10+cCMBPx10')
        wqsubmit('phy','PolyCDM','GBAOx10+SN+cCMBPx10')
        wqsubmit('phy','PolyCDM','GBAOx10+SNx10+cCMBP')

        wqsubmit('phy','PolyCDM','GBAOx10+SN+cCMBPx10')
        wqsubmit('phy','PolyCDM','GBAO+SNx10+cCMBP')
        wqsubmit('phy','PolyCDM','GBAO+SN+cCMBPx10')
        sys.exit(0)

    if sys.argv[1]!='all':
        modell=sys.argv[1].split(',')
    if sys.argv[2]!='all':
        datasetl=sys.argv[2].split(',')

for model in modell:
    for datasets in datasetl:
    #for datasets in ['SN']:
        if datasets in ['BBAO','GBAO','LBAO']:
            prefactl=['pre','phy']
        else:
            prefactl=['phy']
        for prefact in prefactl:
            wqsubmit (prefact,model,datasets)


