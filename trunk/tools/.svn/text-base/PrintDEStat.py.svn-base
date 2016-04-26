#!/usr/bin/env python
##! This prints various DE statistics for the paper
from RunBase import *

if True:
    L=CompositeLikelihood([
            DR11LOWZ(),
            DR11CMASS(),
            DR11LyaAuto(),
            DR11LyaCross()])
#,
 #           SixdFGS()
  #          ])

    #### OLCDM
    T=oLCDMCosmology(kwargs_LCDM={'disable_radiation':True})
    T.setVaryPrefactor()
    L.setTheory(T)
    A=MaxLikeAnalyzer(L, noErrors=True)
    T.printParameters(A.params)
    chi2=L.loglike()*(-2)
    print 'loglike in oLCDM=',L.loglike(), L.theory().Ocb


    #### OLCDM
    T=oLCDMCosmology(zeroDE=True,kwargs_LCDM={'disable_radiation':True})
    T.setVaryPrefactor()
    T.updateParams([Parameter("Om",1.0)])
    L.setTheory(T)
    A=MaxLikeAnalyzer(L, noErrors=True)
    T.printParameters(A.params)
    chi2null=L.loglike()*(-2)
    print 'loglike in oCDM=',L.loglike(), L.theory().Ocb



    print '******************************'
    dchi2=chi2-chi2null
    print 'CHI2 difference:',dchi2, 'sqrt=',sqrt(-dchi2)
    print '******************************'

if False:

    deltaOl=0.003
    T=oLCDMCosmology({'disable_radiation':True})
    T.setVaryPrefactor()
    for n in T.freeParameters():
        print n.name

    ## first we need to get prior ratio
    okmin, okmax=Ok_par.bounds
    ommin, ommax=Om_par.bounds
    Np=10000000
    c=0
    for i in range(Np):
        ok=random.uniform(okmin,okmax)
        om=random.uniform(ommin,ommax)
        ol=1-ok-om
        if (abs(ol)<deltaOl):
            c+=1
        
    priorrat=float(c)/float(Np)
    print 'prior rat=',priorrat
    from cosmich import *
    sys.path.append("Pack")
    import latest
    C=cosmochain(latest.latest_chain+'/oLCDM_pre_BBAO.txt')
    Ol=1.0-C['Om']-C['Ok']
    w=C.chain[:,0]
    
    postrat=w[where(abs(Ol)<deltaOl)].sum()/w.sum()
    print (abs(Ol)<deltaOl).sum(),'aa'
    print w[where(abs(Ol)<deltaOl)].sum(),w.sum(), 'post ratio=',postrat
    dEv=postrat/priorrat
    print 'Ev, log Ev, 1/Ev =',dEv, log(dEv),1/dEv

