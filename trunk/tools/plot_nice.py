#!/usr/bin/env python

from scipy import *
from cosmich import *
import pylab, sys

choice=sys.argv[1]

params1 = {'backend': 'pdf',
               'axes.labelsize': 20,
               'text.fontsize': 18,
               'xtick.labelsize': 20,
               'ytick.labelsize': 20,
               'legend.draw_frame': False,
               'legend.fontsize': 16,
               'lines.markersize': 6,
               'font.size': 20,
               'text.usetex': True}#
pylab.rcParams.update(params1)

#dire='/astro/u/anze/work/April/chains/'
dire='chains/'

if choice=='1':

    C=cosmochain('/home/anze/work/Planck/PLA/base/planck_lowl_lowLike/base_planck_lowl_lowLike','auto')
    C['h']=C['H0*']/100.0
    C.Plot2D('h','omegam*',filled='green',lw=2)

    cosmochain(dire+'LCDM_phy_GBAO.txt').Plot2D('h','Om', filled='red', lims=(0.5,1.0,0.0,1.0),N=60)
    cosmochain(dire+'/LCDM_phy_LBAO.txt').Plot2D('h','Om', filled='blue', lims=(0.5,1.0,0.0,1.0),N=60)
    C=cosmochain(dire+'/LCDM_phy_BBAO.txt')
    C.Plot2D('h','Om', filled='black',N=60)
    print C.GetLimits('h'),'XX1'
    print C.GetLimits('Om'),'XX2'

    pylab.plot([],[],'r-',label='Galaxy BAO')
    pylab.plot([],[],'b-',label='Lyman-$\\alpha$ BAO')
    pylab.plot([],[],'k-',label='Combined BAO')
    pylab.plot([],[],'g-',label='Planck CMB (PLA)')
    pylab.legend(loc='upper left')

    pylab.xlabel('$h$')
    pylab.ylabel('$\Omega_m$')
    pylab.savefig('lcdm.pdf')

if choice=='1.1':

    cosmochain('chains/LCDM_phy_CMBP.txt').Plot2D('h','Om', filled=1)

    C=cosmochain('/home/anze/work/Planck/PLA/base/planck_lowl_lowLike/base_planck_lowl_lowLike','auto')
    i=C.parcol['H0*']
    C.chain[:,i]/=100
    j=C.parcol['omegabh2']
    k=C.parcol['omegach2']
    C.chain[:,k]=(C.chain[:,j]+C.chain[:,k])/C.chain[:,i]**2
    C.Plot2D(i,k,filled='green',lw=2)

    C=cosmochain('/home/anze/work/Planck/PLA/base_Alens/planck_lowl_lowLike/base_Alens_planck_lowl_lowLike','auto')
    i=C.parcol['H0*']
    C.chain[:,i]/=100
    j=C.parcol['omegabh2']
    k=C.parcol['omegach2']
    C.chain[:,k]=(C.chain[:,j]+C.chain[:,k])/C.chain[:,i]**2
    C.Plot2D(i,k,filled='blue',lw=2)

    C=cosmochain('/home/anze/work/Planck/PLA/base/WMAP/base_WMAP','auto')
    i=C.parcol['H0*']
    C.chain[:,i]/=100
    j=C.parcol['omegabh2']
    k=C.parcol['omegach2']
    C.chain[:,k]=(C.chain[:,j]+C.chain[:,k])/C.chain[:,i]**2
    C.Plot2D(i,k,filled='orange',lw=2)


    pylab.xlabel('$h$')
    pylab.ylabel('$\Omega_m$')
    pylab.savefig('lcdm11.pdf')

if choice=='2.0':
    #cosmochain(dire+'/LCDM_phy_BBAO.txt').Plot2D('h','Om', filled='magenta',label="BAO")
    cosmochain(dire+'/LCDM_phy_BBAO+CMBP.txt').Plot2D('h','Om', filled='cyan',label="BAO+Planck")    
    cosmochain(dire+'/LCDM_phy_SN+CMBP.txt').Plot2D('h','Om', filled='orange',label="SN+Planck")    
    cosmochain(dire+'/LCDM_phy_BBAO+SN+CMBP.txt').Plot2D('h','Om', filled='black',label="BAO+SN+Planck")    
    pylab.xlim(0.64,0.74)
    pylab.ylim(0.23,0.37)
    pylab.legend()
    pylab.xlabel('$h$')
    pylab.ylabel('$\\Omega_m$')
    pylab.savefig('deplot20.pdf')

if choice=='2.1':
    #cosmochain(dire+'/oLCDM_phy_BBAO.txt').Plot2D('Om','Ok', filled='magenta',label="BAO")
    cosmochain(dire+'/oLCDM_phy_BBAO+CMBP.txt').Plot2D('Om','Ok', filled='cyan',label="BAO+Planck")    
    cosmochain(dire+'/oLCDM_phy_SN+CMBP.txt').Plot2D('Om','Ok', filled='orange',label="SN+Planck", N=40)    
    cosmochain(dire+'/oLCDM_phy_BBAO+SN+CMBP.txt').Plot2D('Om','Ok', filled='black',label="BAO+SN+Planck")    
    pylab.ylim(-0.05,0.05)
    pylab.legend()
    pylab.xlabel('$\\Omega_m$')
    pylab.ylabel('$\\Omega_k$')
    pylab.savefig('deplot21.pdf')

if choice=='2.2':
    #cosmochain(dire+'/wCDM_phy_BBAO.txt').Plot2D('Om','w', filled='magenta',label="BAO")
    cosmochain(dire+'/wCDM_phy_BBAO+CMBP.txt').Plot2D('Om','w', filled='cyan',label="BAO+Planck",N=30)    
    cosmochain(dire+'/wCDM_phy_SN+CMBP.txt').Plot2D('Om','w', filled='orange',label="SN+Planck",N=30)    
    cosmochain(dire+'/wCDM_phy_BBAO+SN+CMBP.txt').Plot2D('Om','w', filled='black',label="BAO+SN+Planck",N=30)    
    pylab.legend()
    pylab.xlabel('$\\Omega_m$')
    pylab.ylabel('$w$')
    xll,xlh=pylab.xlim()
    pylab.plot ([xll,xlh],[-1,-1],'k:')
    pylab.ylim(-1.4,-0.7)
    pylab.savefig('deplot22.pdf')

if choice=='2.3':
    #cosmochain(dire+'/owCDM_phy_BBAO.txt').Plot2D('w','Ok', filled='magenta',label="BAO")
    cosmochain(dire+'/owCDM_phy_BBAO+CMBP.txt').Plot2D('w','Ok', filled='cyan',label="BAO+Planck")    
    cosmochain(dire+'/owCDM_phy_SN+CMBP.txt').Plot2D('w','Ok', filled='orange',label="SN+Planck")    
    cosmochain(dire+'/owCDM_phy_BBAO+SN+CMBP.txt').Plot2D('w','Ok', filled='black',label="BAO+SN+Planck")    
    pylab.legend()
    pylab.xlabel('$w$')
    pylab.ylabel('$\\Omega_k$')
    pylab.ylim(-0.05,0.05)
    pylab.xlim(-1.6,-0.5)
    pylab.plot ([-1,-1],[-0.05,0.05],'k:')
    pylab.plot ([-1.6,-0.5],[0,0],'k:')

    pylab.savefig('deplot23.pdf')

if choice=='2.4':
    #C=cosmochain(dire+'/owaCDM_phy_BBAO.txt')
    #C['w_z']=C['w']+0.33*C['wa']
    #C.Plot2D('w_z','wa', filled='magenta',label="BAO",N=20)
    D=cosmochain(dire+'/owaCDM_phy_BBAO+CMBP.txt')
    D['w_z']=D['w']+0.33*D['wa']
    D.Plot2D('w_z','wa', filled='cyan',label="BAO+Planck",N=20)    
    E=cosmochain(dire+'/owaCDM_phy_SN+CMBP.txt')
    E['w_z']=E['w']+0.33*E['wa']
    E.Plot2D('w_z','wa', filled='orange',label="SN+Planck",N=20)    
    F=cosmochain(dire+'/owaCDM_phy_BBAO+SN+CMBP.txt')
    F['w_z']=F['w']+0.33*F['wa']
    F.Plot2D('w_z','wa', filled='black',label="BAO+SN+Planck",N=40)
    pylab.legend()
    pylab.xlabel('$w_{(z=0.5)}$')
    pylab.ylabel('$w_a$')
    pylab.xlim(-2.3,-0.5)

    pylab.savefig('deplot24.pdf')



if choice=='3':

    #for da,col,col2,N,lab in [('GBAO','r-','red',40,'Galaxy BAO'),('LBAO','b-','blue',40,'Lyman-$\\alpha$ BAO'),('BBAO','k-',2,40,'Combined BAO')]:
    for da,col,col2,N,lab in [('BBAO','k-',2,40,'Combined BAO')]:
        C=cosmochain(dire+'oLCDM_pre_%s'%da)
        C['Ol']=1-C['Ok']-C['Om']
        pylab.figure(0)
        xx,yy=C.GetHisto('Ol',NormPeak=False)
        print lab
        pylab.plot (xx,yy,col, label=lab,lw=4)

        pylab.figure(1)
        xx,yy=C.GetHisto('Pr',NormPeak=False)
        pylab.plot (xx,yy,col, label=lab,lw=4)

        pylab.figure(2)
        if (da=='BBAO'):
            pylab.figure(2)
            C.Plot2D('Om','Ol', filled='red', lims=(0.0,0.7,0.0,1.3),N=N)
        print da,'XX1',C.GetLimits('Ol')
        C=cosmochain(dire+'oLCDM_pre_%s+Pl_Da'%da)
        C['Ol']=1-C['Om']
        pylab.figure(0)
        xx,yy=C.GetHisto('Ol',NormPeak=False, nbins=100)
        pylab.plot (xx,yy,col.replace('-','--'),lw=2)
        pylab.figure(1)
        xx,yy=C.GetHisto('Pr',NormPeak=False, nbins=100)
        pylab.plot (xx,yy,col.replace('-','--'),lw=2)

        print da,'XX2',C.GetLimits('Ol')
        if False:
            C=cosmochain(dire+'wLCDM_pre_%s.txt'%da)
            C['Ol']=1-C['Om']
            pylab.figure(0)
            xx,yy=C.GetHisto('Ol',NormPeak=False)
            print lab
            pylab.plot (xx,yy,col+'-', lw=4)

            pylab.figure(1)
            xx,yy=C.GetHisto('Pr',NormPeak=False)
            pylab.plot (xx,yy,col+'-', lw=4)


    pylab.figure(0)
    pylab.xlabel('$\\Omega_\\Lambda$')
    pylab.ylabel('${\\rm prob.}$')
    pylab.xlim(-0.2,1.5)
    pylab.ylim(0,22.)
    pylab.legend(loc="upper left")
    pylab.savefig('pre0.pdf')

    pylab.figure(1)
    pylab.xlabel('$c/(H_0 r_d)$')
    pylab.ylabel('${\\rm prob.}$')
    pylab.xlim(18,45)
    pylab.legend()
    pylab.savefig('pre1.pdf')

    pylab.figure(2)
    pylab.plot([0,1],[1,0],'k:')
    pylab.xlabel('$\Omega_m$')
    pylab.ylabel('$\Omega_\Lambda$')
    pylab.legend()
    pylab.savefig('pre2.pdf')

if choice=='4':
    #cosmochain(dire+'/oLCDM_phy_BBAO.txt').Plot2D('Om','Ok', filled='red')
    cosmochain(dire+'/oLCDM_phy_BBAO+CMBW.txt').Plot2D('Om','Ok', filled='blue')    
    cosmochain(dire+'/oLCDM_phy_BBAO+CMBP.txt').Plot2D('Om','Ok', filled='black')    
    cosmochain(dire+'/oLCDM_phy_BBAO+SN+RiessH0.txt').Plot2D('Om','Ok', filled='green')    
    cosmochain(dire+'/oLCDM_phy_SN+CMBP.txt').Plot2D('Om','Ok', filled='orange')    
    pylab.xlabel('$\Omega_m$')
    pylab.ylabel('$\Omega_k$')
    pylab.savefig('olcdm2.pdf')

if choice=="tl":
    C=cosmochain(dire+'/TLight_phy_BBAO+SN.txt')
    C.Plot1D("beta",N=70)
    print C.GetLimits("beta")
    pylab.xlim(0.8,1.2)
    pylab.xlabel("$\\beta$")
    pylab.ylabel("$p(\\beta)$")
    pylab.savefig('tlight.pdf')

    pylab.show()
    


if choice=='5':
    cosmochain(dire+'wCDM_phy_GBAO.txt').Plot2D('Om','w', filled='red', lims=(0.0,1.0,-2.0,0.0),N=30)
    cosmochain(dire+'wCDM_phy_LBAO.txt').Plot2D('Om','w', filled='blue', lims=(0.0,1.0,-2.0,0.0),N=30)
    cosmochain(dire+'wCDM_phy_BBAO.txt').Plot2D('Om','w', filled='black')
    pylab.xlabel('$\Omega_m$')
    pylab.ylabel('$w$')
    pylab.savefig('wcdm.pdf')

if choice=='6':
    cosmochain(dire+'/wCDM_phy_BBAO.txt').Plot2D('Om','w', filled='red')
    cosmochain(dire+'/wCDM_phy_BBAO+CMBW.txt').Plot2D('Om','w', filled='blue')    
    cosmochain(dire+'/wCDM_phy_BBAO+CMBP.txt').Plot2D('Om','w', filled='black')    
    cosmochain(dire+'/wCDM_phy_BBAO+SN.txt').Plot2D('Om','w', filled='green')    
    cosmochain(dire+'/wCDM_phy_SN+CMBP.txt').Plot2D('Om','w', filled='orange')    
    pylab.xlabel('$\Omega_m$')
    pylab.ylabel('$w$')
    pylab.savefig('wcdm2.pdf')

if choice=='7':
    cosmochain(dire+'/nuLCDM_phy_GBAO+CMBP.txt').GetHisto('mnu',nbins=40, plot='r-')
    cosmochain(dire+'/nuLCDM_phy_BBAO+CMBP.txt').GetHisto('mnu',nbins=40, plot='b-')
    cosmochain(dire+'/nuLCDM_phy_BBAO+SN+CMBP.txt').GetHisto('mnu',nbins=40, plot='k-')


    pylab.xlabel('$\\sum m_\\nu$')
    pylab.ylabel('prob.')

    pylab.savefig('nu.pdf')


def plotH(y,chainname,color, txt):
    
    if "PLA" in chainname:
        C=cosmochain(chainname,'auto')
        C['h']=C['H0*']/100
        m,p1,p2,p3,m1,m2,m3=C.GetLimits('h')
    elif "N:" in chainname:
        m,e=map(float,chainname.split()[1:])
        p1,p2,p3=e,2*e,3*e
        m1,m2,m3=p1,p2,p3
    else:
        C=cosmochain(dire+chainname)
        m,p1,p2,p3,m1,m2,m3=C.GetLimits('h')
    pylab.errorbar(m,y,xerr=[[m1],[p1]],fmt='-o'+color,lw=4,capthick=4, capsize=10,ms=10)
    #pylab.errorbar(m,y,xerr=[[m2],[p2]],fmt='-'+color,lw=1,capthick=2,capsize=5)
    #pylab.errorbar(m,y,xerr=[[m3],[p3]],fmt='-'+color,lw=1,capthick=2,capsize=5)
    print txt,(p1+m1)/(2*m), m,p1,p2,p3, m1,m2,m3
    pylab.text(m+p1+0.005,y-0.05,txt,verticalalignment='baseline')

if choice=='h':

    plotH(0.8,'LadoCDM_phy_GBAO+SN+cCMBP.txt','b', '$\\mbox{PolyCDM GBAO+SN}+\mathbf{r_d}$')
    plotH(1,'LadoCDM_phy_BBAO+SN+CMBP.txt','b', 'PolyCDM BAO+SN+Planck')
    plotH(1.2,'owaCDM_phy_BBAO+SN+CMBP.txt','b','o$w_0w_a$-CDM BAO+SN+Planck')
#    plotH(1.4,'LadoCDM_phy_BBAO+SN+CMBW.txt','r','PolyCDM BAO+SN+WMAP')
#    plotH(1.6,'owaCDM_phy_BBAO+SN+CMBW.txt','r','o$w_0w_a$CDM BAO+SN+WMAP')

    plotH(1.6,'/home/anze/work/Planck/PLA/base/planck_lowl_lowLike/base_planck_lowl_lowLike','g', '$\Lambda$CDM PLA Planck')
    plotH(1.8,'/home/anze/work/Planck/PLA/base/WMAP/base_WMAP','g', '$\Lambda$CDM PLA WMAP')

    plotH(2.2,'N: 0.738 0.024','m','Riess++')
    plotH(2.4,'N: 0.743 0.021','m','Freedmann++')
    plotH(2.6,'N: 0.725 0.025','m','Efstathiou')
    #plotH(3.0,'N: 0.72 0.03','m','Efstathiou NGC4258')

    pylab.xlim(0.60,0.9)
    pylab.ylim(0.6,2.8)
    pylab.grid(axis='both')
    pylab.yticks([])
    pylab.xlabel("$h$")
    pylab.savefig('Ayaaaa-boom.pdf')
    pylab.show()

