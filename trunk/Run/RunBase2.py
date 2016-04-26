#Add paths, we want to be able to run in either root or Run/
import sys,os
sys.path.append("py")
sys.path.append("../py")


#Cosmologies
from LCDMCosmology import *
from oLCDMCosmology import *
from wLCDMCosmology import *
from PolyCDMCosmology import *

#JAV
from owa0CDMCosmology import *
from JordiCDMCosmology import *
from WeirdCDMCosmology import *
from TLightCosmology import *
from SplineLCDMCosmology import *
from DecayLCDMCosmology import *
from StepCDMCosmology import *
from EarlyDECosmology import *

#Like modules
from BOSSLikelihoods import *
from SimpleCMB import *
from CompressedSNLikelihood import *
from HubbleParameterLikelihood import *

#Composite Likelihood
from CompositeLikelihood import *

#Likelihood Multiplier
from LikelihoodMultiplier import *

#Analyzers
from MaxLikeAnalyzer import *
from MCMCAnalyzer import *


## String parser Aux routines
model_list="LCDOM, LCDMasslessnu, nuLCDM, NnuLCDM, noradLCDM, oLCDM, wCDM, waCDM, owCDM,"\
           "owaCDM, JordiCDM, WeirdCDM, TLight, Spline, PolyCDM"

def ParseModel(model):
#    if model=="LCDM_2":
    if model=="LCDM" or "LCDM_" in model:
        T=LCDMCosmology()
    elif model=="LCDMmasslessnu" or "LCDMmasslessnu_" in model:
        T=LCDMCosmology(mnu=0)
    elif model=="nuLCDM" or "nuLCDM_" in model:
        T=LCDMCosmology()
        T.setVaryMnu()
    elif model=="NnuLCDM" or "NnuLCDM_" in model:
        LCDMCosmology.rd_approx="tabulated_Nnu"
        T=LCDMCosmology()
        T.setVaryNnu()
    elif model=="noradLCDM" or "noradLCDM_" in model:
        T=LCDMCosmology(disable_radiation=True)
    elif model=="oLCDM" or "oLCDM_" in model:
        T=oLCDMCosmology()
    elif model=="wCDM" or "wCDM_" in model:
        T=wLCDMCosmology()
    elif model=="waCDM" or "waCDM_" in model:
        T=owa0CDMCosmology(varyOk=False)
    elif model=="owCDM" or "owCDM_" in model:
        T=owa0CDMCosmology(varywa=False)
    elif model=="owaCDM" or "owaCDM_" in model:
        T=owa0CDMCosmology()
    elif model=="JordiCDM" or "JordiCDM_" in model:
        T=JordiCDMCosmology()
    elif model=="WeirdCDM" or "WeirdCDM_" in model:
        T=WeirdCDMCosmology()
    elif model=="TLight" or "TLight_" in model:
        T=TLightCosmology()
    elif model=="StepCDM" or "StepCDM_" in model:
        T=StepCDMCosmology()
    elif model=="Spline" or "Spline_" in model:
        T=SplineLCDMCosmology()
    elif model=="Decay" or "Decay_" in model:
        T=DecayLCDMCosmology() 
    elif model=="PolyCDM" or "PolyCDM_" in model:
        T=PolyCDMCosmology()
    elif model=="EarlyDE" or "EarlyDE_" in model:
        T=EarlyDECosmology()
    else:
        print "Cannot recognize model", model
        sys.exit(1)

    return T

data_list="BBAO, GBAO, GBAO_no6dF, LBAO, LaBAO, LxBAO, CMBP, cCMBP, CMBW, SN, RiessH0, 6dFGS"
def ParseDataset(datasets):
    dlist=datasets.split('+')
    if 'BBAO' in dlist:
        L=CompositeLikelihood([
        DR11LOWZ(),
        DR11CMASS(),
        DR11LyaAuto(),
        DR11LyaCross(),
        SixdFGS()
        ])
        dlist.remove('BBAO')
    elif 'GBAO' in dlist:
        L=CompositeLikelihood([
        DR11LOWZ(),
        DR11CMASS(),
        SixdFGS()
        ])
        dlist.remove('GBAO')
    elif 'GBAOx10' in dlist:
        L=CompositeLikelihood([
        LikelihoodMultiplier(DR11LOWZ,DR11LOWZ(),10.0),
        LikelihoodMultiplier(DR11CMASS,DR11CMASS(),10.0),
        LikelihoodMultiplier(SixdFGS,SixdFGS(),10.0)
        ])
        dlist.remove('GBAOx10')
    elif 'GBAO_no6dF' in dlist:
        L=CompositeLikelihood([
        DR11LOWZ(),
        DR11CMASS()
        ])
        dlist.remove('GBAO_no6dF')
    elif 'LBAO' in dlist:
        L=CompositeLikelihood([
        DR11LyaAuto(),
        DR11LyaCross()
        ])
        dlist.remove('LBAO')

    elif 'LaBAO' in dlist:
        L=CompositeLikelihood([
        DR11LyaAuto(),
        ])
        dlist.remove('LaBAO')

    elif 'LxBAO' in dlist:
        L=CompositeLikelihood([
        DR11LyaCross(),
        ])
        dlist.remove('LxBAO')
    else:
        L=CompositeLikelihood([])

    if 'CMBP' in dlist:
        L.addLikelihood(PlanckLikelihood())
        dlist.remove('CMBP')

    if 'cCMBP' in dlist:
        L.addLikelihood(cPlanckLikelihood())
        dlist.remove('cCMBP')
 
    if 'cCMBPx10' in dlist:
        L.addLikelihood(LikelihoodMultiplier(cPlanckLikelihood, cPlanckLikelihood(),10.0))
        dlist.remove('cCMBPx10')

    if 'CMBW' in dlist:
        L.addLikelihood(WMAP9Likelihood())
        dlist.remove('CMBW')

    if 'SN' in dlist:
        L.addLikelihood(BetouleSN())
        dlist.remove('SN')

    if 'SNx10' in dlist:
        L.addLikelihood(LikelihoodMultiplier(BetouleSN,BetouleSN(),10.0))
        dlist.remove('SNx10')

    if 'RiessH0' in dlist:
        L.addLikelihood(RiessH0())
        dlist.remove('RiessH0')
    
    if '6dFGS' in dlist:
        L.addLikelihood(SixdFGS())
        dlist.remove('6dFGS')    
    
    if len(dlist)>0:
        print "Cannot parse data, unrecognizable part:", dlist
        sys.exit(1)

    return L

