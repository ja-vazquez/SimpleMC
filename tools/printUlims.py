#!/usr/bin/env python
import cosmich
import sys
ch = cosmich.cosmochain(sys.argv[1])
for nm in ch.paramnames:
    res = ch.GetLimits(
        nm, limlist=[0, 0, 0, 0.95, 0.997, 0.999], returnlims=True)
    print(nm, res)

ch['Ol'] = 1-ch['Ok']-ch['Om']
print("Ol", ch.GetLimits("Ol", lims=[0, 0, 0, 0.95, 0.997, 1.0]))
