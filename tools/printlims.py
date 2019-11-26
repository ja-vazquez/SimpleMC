#!/usr/bin/env python
import cosmich
import sys
ch = cosmich.cosmochain(sys.argv[1])
for nm in ch.paramnames:
    res = ch.GetLimits(nm)
    print(nm, res)

ch['Ol'] = 1-ch['Ok']-ch['Om']
print("Ol", ch.GetLimits("Ol"))
