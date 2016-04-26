#!/usr/bin/env python
import pylab


with open('del.out') as inf:
    for line in inf:
        parts = line.split() # split line into parts
        if len(parts) > 1:   # if at least 2 parts/columns
#            if '19363' in parts:
                print 'wq rm '+ parts[0]   # print column 2
