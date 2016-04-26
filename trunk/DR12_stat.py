
import numpy as np

dir = 'data/'
files = 'combxiFBBAOlikgrid_zb'

for i in range(3):
  apar, aper = [], []
  with open(dir+files + '%d.dat'%(i+1), 'r') as f:
    for lines in f:
        vals =lines.split()[0:]
        apar.append(float(vals[0]))
        aper.append(float(vals[1]))

  print 'apar',np.mean(apar), np.std(apar)
  print 'aper',np.mean(aper), np.std(aper)
  print '--'*10
