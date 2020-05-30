#!/usr/bin/env python

# this scripts runs everything on the bnl cluster
import os
import sys
from wqdriver import wqsubmit

if len(sys.argv) > 1:
    lst = sys.argv[1:]
else:
    lst = ['S3', 'S5', 'HLadder', 'smoking',
           'Early', 'poly', 'SlowRDE', 'step']

for l in lst:
    if l == 'wDM':
        wqsubmit('phy', 'wDM', 'BBAO+PlRd', 3, 30000)

    if l == 'Quint':
        for datasets in ['BBAO+Planck', 'BBAO']:
            wqsubmit('phy', 'Quint_last', datasets, 5, 60000)

    if l == 'S3':
        for datasets in ['GBAO', 'LBAO', 'BBAO']:
            wqsubmit('phy', 'LCDM', datasets, 5, 60000)
            wqsubmit('pre', 'LCDM', datasets)
            wqsubmit('pre', 'oLCDM', datasets, 5, 60000)
            wqsubmit('pre', 'oLCDM', datasets+'+PlDa', 5, 60000)

    elif l == 'S5':
        for model in ['LCDM', 'oLCDM', 'wCDM', 'owCDM', 'waCDM', 'owaCDM', 'LCDMmasslessnu', 'nuLCDM', 'nuoLCDM']:
            for datasets in ['BBAO+Planck', 'SN+Planck', 'BBAO+SN+Planck']:
                wqsubmit('phy', model, datasets, 5, 60000)

    elif l == 'onlyBAO':
        for model in ['LCDM', 'oLCDM', 'wCDM', 'owCDM', 'waCDM', 'owaCDM']:
            for datasets in ['BBAO']:
                wqsubmit('phy', model, datasets, 5, 60000)

    elif l == "mnuw":
        for model in ['nuwCDM']:
            for datasets in ['BBAO+Planck', 'SN+Planck', 'BBAO+SN+Planck']:
                wqsubmit('phy', model, datasets, 5, 60000)

    elif l == 'HLadder':
        for data in ['GBAO+SN+PlRd', 'GBAO+SNx10+PlRdx10', 'GBAOx10+SN+PlRdx10', 'GBAOx10+SNx10+PlRd', 'GBAOx10+SN+PlRd',
                     'GBAO+SNx10+PlRd', 'GBAO+SN+PlRdx10']:
            wqsubmit('phy', 'PolyCDM', data)
        wqsubmit('phy', 'owaCDM', 'GBAO+SN+PlRd')
        wqsubmit('phy', 'PolyCDM', 'GBAO+UnionSN+PlRd')
        wqsubmit('phy', 'fPolyCDM', 'GBAO+UnionSN+PlRd')
        wqsubmit('phy', 'owaCDM', 'GBAO+UnionSN+PlRd')

    elif l == 'smoking':
        for datasets in ['LBAO', 'BBAO', 'BBAO+Planck', 'BBAO+SN+Planck', 'Planck']:
            for models in ['Decay', 'DecayFrac', 'WeirdCDM']:
                wqsubmit('phy', models, datasets)

    elif l == 'decayfraconly':
        for datasets in ['BBAO', 'Planck', 'BBAO+Planck', 'BBAO+SN+Planck']:
            for models in ['Decay01', 'Decay05']:
                wqsubmit('phy', models, datasets, 8, 30000)

    elif l == 'decaybig':
        for datasets in ['BBAO+Planck', 'BBAO+SN+Planck']:
            #      for models in ['Decay','Decay01', 'Decay05']:
            for models in ['Decay01']:
                wqsubmit('phy', models, datasets, 24, 100000)
    elif l == 'decaybig2':
        for datasets in ['BBAO+Planck']:
            for models in ['Decay01']:
                wqsubmit('phy', models, datasets, 64, 10000)

    elif l == 'decaytest':
        for datasets in ['LBAO+Planck']:
            for models in ['Decay']:
                wqsubmit('phy', models, datasets, 64, 10000)

    elif l == 'poly':
        wqsubmit('phy', 'PolyCDM', 'BBAO+SN+Planck')
        wqsubmit('phy', 'PolyCDM', 'BBAO+SN+PlRd')
        wqsubmit('phy', 'TLight', 'BBAO+SN+PlRd')

    elif l == 'smokingW':  # weird alone
        for datasets in ['LBAO', 'BBAO', 'BBAO+Planck', 'BBAO+SN+Planck', 'Planck']:
            for models in ['WeirdCDM']:
                wqsubmit('phy', models, datasets)

    elif l == 'Early':
        for model in ['EarlyDE', 'EarlyDE_rd_DE']:
            for datasets in ['BBAO+SN+Planck']:
                wqsubmit('phy', model, datasets, 5, 60000)

    elif l == 'SlowRDE':
        for model in ['SlowRDE']:
            for datasets in ['BBAO+Planck', 'SN+Planck', 'BBAO+SN+Planck']:
                wqsubmit('phy', model, datasets, 5, 60000)

    elif l == 'Neff':
        wqsubmit('phy', 'NeffLCDM', 'LBAO+Planck')
        wqsubmit('phy', 'NeffLCDM', 'BBAO+Planck')
        wqsubmit('phy', 'NeffLCDM', 'BBAO+SN+Planck')

    elif l == "step":
        wqsubmit('phy', 'StepCDM', 'BBAO+SN+Planck', 8, 30000)
        wqsubmit('phy', 'StepCDM', 'BBAO+Planck', 8, 30000)

    else:
        print("WTF", l)
