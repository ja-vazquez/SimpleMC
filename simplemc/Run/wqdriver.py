#!/usr/bin/env python

# this scripts runs everything on the bnl cluster
import os
import sys
import time

# datasetl=['BBAO','GBAO','GBAO_no6dF','LBAO','LBAO+CMBP','LBAO+CMBW','GBAO+CMBP','GBAO+CMBW',
#         'BBAO+CMBP','BBAO+CMBW','BBAO+SN','BBAO+SN+RiessH0','SN','SN+CMBW','SN+CMBP',
#         'BBAO+SN+CMBP','BBAO+SN+CMBW', 'GBAO+SN+cCMBP']


# removed JordiCDM for the time being
modell = ['LCDM', 'oLCDM', 'wCDM', 'owCDM', 'waCDM',
          'owaCDM', 'PolyCDM', 'LCDMmasslessnu', 'nuLCDM']


def wqsubmit(prefact, model, datasets, nchains=10, nsamp=100000, chainsdir="chains"):
    # run three chains
    for i in range(1, nchains+1):

        comd = "Run/driver.py %s %s %s %i %i %i %s" % (
            prefact, model, datasets, i, min(100, nsamp/100), nsamp, chainsdir)
        print(comd)
        nm = model+"_"+prefact+"_"+datasets+"_"+str(i)
        wqcom = 'wq sub -r job_name:%s -c "%s"' % (nm, comd)
        #os.system("nohup %s&"%(wqcom))
        os.system("nohup %s >%s 2>%s &" % (wqcom, chainsdir +
                                           '/logs/'+nm+'.log', chainsdir+'/logs/'+nm+'.err'))
        time.sleep(0.5)


if __name__ == "__main__":
    if sys.argv[1] != 'all':
        modell = sys.argv[1].split(',')
    if sys.argv[2] != 'all':
        datasetl = sys.argv[2].split(',')

    for model in modell:
        for datasets in datasetl:
            prefactl = ['phy']
            for prefact in prefactl:
                wqsubmit(prefact, model, datasets)
