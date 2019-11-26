#!/usr/bin/env python

# this scripts runs everything on the bnl cluster
import os
import sys
from wqdriver import wqsubmit

if len(sys.argv) > 1:
    lst = sys.argv[1:]
else:
    lst = ['TakadaFlat']


for l in lst:
    if l == 'TakadaFlat':
        wqsubmit('pre', 'PolyOk', 'BBAO+SN', 5, 30000)
        wqsubmit('pre', 'PolyOkc', 'BBAO+SN', 5, 30000)
        wqsubmit('pre', 'PolyOkf', 'BBAO+SN', 5, 30000)
