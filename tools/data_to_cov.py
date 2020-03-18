
#Given a file with data and error bars, create its covariance matrix

import pandas as pd
import numpy  as np

file = '../data/HDiagramCompilacion-data_31.txt'

df   = pd.read_csv(file, sep ='\s+', names = ['z', 'Hz', 'error'], skiprows=[0,1,2])
sigma = np.array(df['error'].tolist())**2
cov = np.zeros((31, 31), float)
np.fill_diagonal(cov, sigma)

with open('HDiagramCompilacion-cov_31.txt', 'wb') as f:
    f.write('# 31 31 \n')
    np.savetxt(f, cov, fmt='%.2f')
