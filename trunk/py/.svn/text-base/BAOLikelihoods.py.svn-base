#
# The BAO likelihoods.
#

from TabulatedBAOLikelihood import *
from TabulatedBAODVLikelihood import *
from GaussBAODVLikelihood import *
from LCDMCosmology import *
from scipy import *


class DR11LOWZ(GaussBAODVLikelihood):
    def __init__(self):
        obh2=0.0224; Om=0.274; h=0.7; mnu=0
        fidTheory=LCDMCosmology(obh2,Om,h,mnu)
        GaussBAODVLikelihood.__init__(self,"DR11LOWZ",0.32,1264.0,25.0, fidTheory)

class DR11CMASS(TabulatedBAOLikelihood):
    def __init__(self):
        ## fiducial cosmology for LOWZ/CMASS data.
        ## see Anderson et al, page 28
        obh2=0.0224; Om=0.274; h=0.7; mnu=0 #rd=149.28
        fidTheory=LCDMCosmology(obh2,Om,h,mnu)
        ## negative col means the cols is probability and not chi2
        TabulatedBAOLikelihood.__init__(self,"DR11CMASS",'data/sdss_DR11CMASS_consensus.dat', 
                                   -2, fidTheory, 0.57 )

class DR11LyaAuto(TabulatedBAOLikelihood):
    def __init__(self):
        ## fiducial cosmology for Lya data.
        ## see e.g. Busca's email on 12/3/13

        obh2=0.0227; Om=0.27; h=0.7; mnu=0.06 #;# rd=149.77
        fidTheory=LCDMCosmology(obh2, Om, h, mnu)
        # File from 5/16 from Nicolas.
        TabulatedBAOLikelihood.__init__(self,"DR11LyaAuto",'data/chi2_surface_dr11_baseline_fit.txt', 
                                   4, fidTheory, 2.34 )

class DR11LyaCross(TabulatedBAOLikelihood):
    def __init__(self):
        obh2=0.0227; Om=0.27; h=0.7; mnu=0#;# rd=149.77
        fidTheory=LCDMCosmology(obh2, Om, h, mnu)
        TabulatedBAOLikelihood.__init__(self,"DR11LyaCross",'data/lyabaocross.scan', 
                                   2, fidTheory, 2.36 )

#Data extracted from http://arxiv.org/pdf/1106.3366.pdf 
class SixdFGS(GaussBAODVLikelihood):
    def __init__(self):
        obh2=0.02227; Om=0.27; h=0.7; mnu=0
        fidTheory=LCDMCosmology(obh2,Om,h,mnu)
        GaussBAODVLikelihood.__init__(self,"SixdFGS",0.106,456.0,27.0, fidTheory, maxchi2=4)

# SDSS Main Galaxy Sample BAO
class SDSSMGS(TabulatedBAODVLikelihood):
    def __init__(self):
        obh2=0.021547; Om=0.31; h=0.67; mnu=0
        fidTheory=LCDMCosmology(obh2,Om,h,mnu)
        TabulatedBAODVLikelihood.__init__(self,"MGS","data/chidavexi8stavePk5staverec.dat", fidTheory,0.15)

