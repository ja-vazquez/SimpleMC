#!/usr/bin/env python
import pylab

file ='/astro/astronfs01/workarea/jvazquez/Getdist/waCDM_phy_BBAO+SN+Planck.margestats'

with open(file) as inf:
    for line in inf:
        parts = line.split() # split line into parts
        if len(parts) > 1:   # if at least 2 parts/columns
                print parts[0], parts[1], parts[2]   # print column 2


#file = open(file, 'r')
#lines = file.readlines()
#file.close()



#for line in lines:
#    parts = line.split() # split line into parts
#    if len(parts) > 1:
#        column2 = parts[1]
#        columnLength = len(column2)
#        desireValue = column2[columnLength-3]
#        print desireValue


