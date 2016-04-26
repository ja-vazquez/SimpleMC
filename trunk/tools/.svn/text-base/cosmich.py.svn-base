#!/usr/bin/env python
#
# Finally a getdist improvement
#
#

from glob import *
from scipy import *
import pylab
import numpy.fft as fft
from  scipy.interpolate import UnivariateSpline

def myloadtxt(fname):
    ## fixed the annoying "setting an array element with a sequence."
    ## when chain still runnig.
    da=open(fname).readlines()[:-1]
    da=array([map(float, line.split()) for line in da])
    return da

class cosmochain:
    
    def __init__ (self,root, nums='auto', skip_=None,temp=1.0, balance=True, weightfunc=None, kirkby=False):

        if (balance and ("PLA" in root or "neffch" in root)):
            balance=False

        ## get file list
        if type(nums)==type('auto'):
            if nums=='auto':
                flist = glob (root+'_?.txt')+glob (root+'_??.txt')
        else:
            if (nums==None):
                flist=[root]
            elif (type(nums[0])==type(1)):
                flist = map (lambda x: root+'_'+str(x)+'.txt',nums)

        if len(flist)==0:
            print "Bad chain spec."
            print root+'_?.txt'
        
        ## get parameter names
        try:
            self.paramnames=[n.split()[0] for n in open(flist[0].replace('.txt','.paramnames')).readlines()]
        except:
            try:
                lines=open(root+'.paramnames').readlines()
                self.paramnames=[n.split()[0] for n in lines]
                self.latexnames=[' '.join(n.split()[1:]) for n in lines]
                self.lname={}
                for n,ln in zip(self.paramnames,self.latexnames):
                    self.lname[n]=ln
            except:
                 try:
                     self.paramnames=[n.split()[0] for n in open(root+'.params').readlines()]
                 except:
                     print "neither params nor paramnames"
                     self.paramnames=[]
        self.parcol={}
        for i,n in enumerate(self.paramnames):
            print i,n
            self.parcol[n]=i+2
        print "Got ",len(self.paramnames), "parameters."
        data=[]

        for fname in flist:
            print "Reading", fname,"...",
            if (kirkby):
                skip=3
                print "kirkby style ",
                cdata=open(name).readlines()[skip:-1]
                cdata=array([ [1,0]+map(float,x.split()) for x in cdata])

            else:
                da=myloadtxt(fname)
                if (skip_==None):
                    finlike = da[-1,1]
                    ii=0
                    while (da[ii,1]>finlike):
                        ii+=1
                    #plus 10
                    skip=ii+20
                    cdata=da[skip:-1]
                else:
                    skip=skip_
                print "skipping:",skip
                if (balance):
                    print "balancing..."
                    ## make weights the same in average to avoid issues with 
                    ## temperature
                    cdata[:,0]/=cdata[:,0].mean()
                cdata=da[skip:-1]

            data=data+list(cdata)
            
        self.chain=array(data)

        if weightfunc!=None:
            print "Reweighting"
            for i,line in enumerate(self.chain):
                self.chain[i,0]=weightfunc(line)
                
        print len(self.chain), len(self.chain[0])
        del data
        self.N=len(self.chain)

        try:
            self.bestarg = self.chain[:,1].argmin()
            self.best = self.chain[self.bestarg]
        except:
            print "WTF?"

        if (temp<>1):
            like = self.chain[:,1]
            like = like-like.min()
            self.chain[:,0] *= exp(-(temp-1.0)*like)
            
    def latexname(self,name):
        print self.lname[name]
        return '$'+self.lname[name]+'$'

    def slist(self,name):
        return self.chain[:,self.parcol[name]]

    def __getitem__ (self,key):
        if (type(key)==type("p")):
            key=self.parcol[key]
        return self.chain[:,key]

    def __setitem__ (self,key,res):
        if (type(key)==type("p")):
            if self.parcol.has_key(key):
                key=self.parcol[key]
            else:
                ## need to add a column
                N=len(self.chain)
                self.chain=concatenate((self.chain,zeros((N,1))),1)
                Nc=len(self.chain[0])-1
                self.parcol[key]=Nc
                self.paramnames.append(key)
                key=Nc
                
        self.chain[:,key]=res



    def BestSample(self):
        ii=self.chain[:,1].argmin()
        return self.chain[ii,:]

    def Plot1D(self, p1, sty='r-', label="",N=50):
        xx,yy=self.GetHisto(p1, nbins=N)
        pylab.plot(xx,yy,sty, label="",lw=2)


    def Plot2D(self, p1, p2, N=60, lims=None, conts=[0.68,0.95,0.997], filled=True, lw=2, bp=False, blur=None, nch1=None,nch2=None, label=""):
        pl =zeros((N,N))


        if (nch1==None):
            if (type(p1)==type("p")):
                p1=self.parcol[p1]
            xx=self.chain[:,p1]
        else:
            xx=nch1
            
        if (nch2==None):
            if (type(p2)==type("p")):
                p2=self.parcol[p2]
            yy=self.chain[:,p2]
        else:
            yy=nch2

        we = self.chain[:,0]

        if (lims==None):
            xmin=xx.min()
            xmax=xx.max()
            de = (xmax-xmin)/100
            if (de==0): 
                de=xmax/100
            xmin -= de
            xmax += de


            ymin=yy.min()
            ymax=yy.max()
            de = (ymax-ymin)/100
            if (de==0):
                de=ymax/100
            ymin -= de
            ymax += de
        else:
            xmin,xmax, ymin,ymax=lims

        out=0
        for x, y,w in zip(xx,yy,we):
            i1=int((x-xmin)/(xmax-xmin)*N)
            i2=int((y-ymin)/(ymax-ymin)*N)
            try:
                pl[i2,i1]+=w
            except:
                out+=w

        if (out>0):
            print "warn: out =", out/we.sum()
        

        b=pl.flatten()

        b=b.tolist()
        b.sort(reverse=True)
        b=array(b)

        c=b*1.0

        c=c.cumsum()
        c/=c[-1]

        l1=0
        l2=0
        l3=0

        for val,cum in zip(b,c):
            if (cum>conts[0]) and (l1==0):
                l1=val
            if (cum>conts[1]) and (l2==0):
                l2=val
            if (cum>conts[2]) and (l3==0):
                l3=val

        print l1,l2,l3

        #pl=smline2(pl)

        #pylab.imshow(pl,extent=(xmin,ymin,xmax,ymax),origin='lower',aspect='auto')

        if type(filled)==type('string'):
            print 'lw=',lw
            pylab.contour(pl,[l3,l2,l1],extent=(xmin,xmax,ymin,ymax),origin='lower',aspect='auto', colors=filled, linewidths=lw)
            if (label!=""):
                pylab.plot ([],[],color=filled,linewidth=lw,label=label)
	#JAV
            if filled == "green":
                pylab.contourf(pl,[l3,l2,l1],extent=(xmin,xmax,ymin,ymax),origin='lower',aspect='auto', cmap=pylab.get_cmap('Greens'))
                pylab.contour(pl,[l3,l2,l1],extent=(xmin,xmax,ymin,ymax),origin='lower',aspect='auto', colors='green', lwidth=lw)

        else:
            if filled==1:
                pylab.contourf(pl,[l3,l2,l1],extent=(xmin,xmax,ymin,ymax),origin='lower',aspect='auto', cmap=pylab.get_cmap('cool'))
                pylab.contour(pl,[l3,l2,l1],extent=(xmin,xmax,ymin,ymax),origin='lower',aspect='auto', colors='red', lwidth=lw)
            elif filled==2:
                pylab.contourf(pl,[l3,l2,l1],extent=(xmin,xmax,ymin,ymax),origin='lower',aspect='auto', color='red', cmap=pylab.get_cmap('autumn'))
                #pylab.imshow(pl,extent=(xmin,xmax,ymin,ymax),origin='lower', aspect='auto', cmap=pylab.get_cmap('autumn'))
            elif filled==0:
                pylab.contour(pl,[l3,l2,l1],extent=(xmin,xmax,ymin,ymax),origin='lower',aspect='auto', colors='black', lwidth=lw)
            elif filled==-1:
                pylab.contour(pl,[l3,l2,l1],extent=(xmin,xmax,ymin,ymax),origin='lower',aspect='auto', colors='red', lwidth=lw)
            elif filled==-2:
                pylab.contour(pl,[l3,l2,l1],extent=(xmin,xmax,ymin,ymax),origin='lower',aspect='auto', colors='blue', lwidth=lw)
            elif filled==-3:
                pylab.contour(pl,[l3,l2,l1],extent=(xmin,xmax,ymin,ymax),origin='lower',aspect='auto', colors='green', lwidth=lw)

        if (bp):
            pylab.plot(self.best[p1], self.best[p2],'ro')

        return (xmin,xmax, ymin,ymax)


    def GetLimits (self, param, ML=False, nch=None, limlist=[0.5-0.997/2,0.5-0.95/2, 0.5-0.68/2,0.5,0.5+0.68/2, 0.5+0.95/2,0.5+0.997/2],
                    returnlims=False):
        if nch!=None:
            lis = zip(nch, self.chain[:,0])
        else:
            if (type(param)==type("st")):
                param=self.parcol[param]
            lis = zip(self.chain[:,param], self.chain[:,0])

        lis.sort()
        lis=array(lis)
        pars= array(lis[:,0])
        wei = array(lis[:,1])
        wei = wei.cumsum()
        wei = wei/wei[-1]
        lims =[]
        for lim in limlist:
            a=where(wei>lim)[0][0]
            lims.append(pars[a])
        if (ML):
            print "Central moved:", lims[3],
            if nch!=None:
                lims[3]=nch[self.bestarg]
            else:
                lims[3]=self.best[param]
            print "->", lims[3]

        if returnlims:
            return lims
        m,ep1,ep2,ep3, em1,em2,em3=lims[3],lims[4]-lims[3],lims[5]-lims[3],lims[6]-lims[3],lims[3]-lims[2],lims[3]-lims[1],lims[3]-lims[0]
        return (m,ep1,ep2, ep3, em1,em2,em3)

    def GetMarginal (self, param, val):
        return ((self.chain[:,param]>val)*self.chain[:,0]).sum()/self.chain[:,0].sum()
    
            

    
    def GetHisto (self,param, nbins=50, nch=None, mnval=None, mxval=None, smooth=None, NormPeak=False, plot=None, lw=4):
        "Returns histogram for plotting."
        if nch!=None:
            line = nch
        else:
            if (type(param)==type("st")):
                param=self.parcol[param]

            line = self.chain[:,param]

        #print self.chain[0,:]

        if mnval==None:
            mnval = line.min()
        if mxval==None:    
            mxval = line.max() #(to add the last one)
        if (mnval==mxval):
            return None,None
        mxval*=1.001
        #print mnval, mxval
        stp = (1.0*mxval-1.0*mnval)/nbins
        tmp = map(int,(line-mnval)/stp)

        histo = zeros((nbins,))


        for ii in range(len(tmp)):
            if (tmp[ii]>=0) and (tmp[ii]<nbins):
                histo[tmp[ii]]+=self.chain[ii,0]
        xval = array(map(lambda x:mnval+(x+0.5)*stp,range(nbins)))
        yval = histo/stp
        print xval, mnval,mxval,stp

        if smooth:
            yvalpad = array([0,0,0,0]+list(yval)+[0,0,0,0])
            if smooth==1:
                yval = (yvalpad[4:nbins+4]+yvalpad[3:nbins+3]+yvalpad[5:nbins+5])/3.0
            if smooth==2:
                yval = (yvalpad[4:nbins+4]+yvalpad[3:nbins+3]+yvalpad[5:nbins+5]+yvalpad[2:nbins+2]+yvalpad[6:nbins+6])/5.0

        if (NormPeak):
            yval/=yval.max()
        else:
            area = yval.sum()*stp
            yval/=area
            

        if (plot!=None):
            pylab.plot(xval,yval, plot, lw=lw)
        return xval,yval
    
    def GetLimitsOld (self,param,lims=[0.6826894920,0.9544997360,0.9973002039]):
        "returns median and pairs corresponding to lims"
        line = map (lambda x:(x[param+1],x[0]),self.chain)
        #print line
        line.sort()
        sw=self.chain[:,0].sum()
        lowl = [None]*len(lims)
        highl = [None]*len(lims)
        lowtrig = 0.5-array(lims)/2
        hightrig = 0.5+array(lims)/2

        #print sw, hightrig

        sum = 0
        mean = None
        for qq in line:
            sum+=qq[1]
            if (not mean) and (sum>=0.5*sw):
                mean=qq[0]
            for ii in range(len(lims)):
                if (not lowl[ii]) and (sum>=lowtrig[ii]*sw):
                    lowl[ii]=qq[0]
                if (not highl[ii]) and (sum>=hightrig[ii]*sw):
                    highl[ii]=qq[0]
        lowl-=mean
        highl-=mean
        str = "%2.0f^{+%2.0f+%2.0f+%2.0f}_{%2.0f%2.0f%2.0f}"%(mean,highl[0],highl[1],highl[2],lowl[0],lowl[1],lowl[2])
        print str


        return mean, lowl, highl

    def GetCovariance(self,parlist):
        N=len(parlist)
        cov=zeros((N,N))
        mn=zeros(N)
        sw=0.0
        for el in self.chain:#[:100]:
            sw+=el[0]
            vec=array([el[i] for i in parlist])
            mn+=vec*el[0]
            cov+=el[0]*array([[v1*v2 for v1 in vec] for v2 in vec])
            #print cov[0,0], vec[0],sw

        mn/=sw
        cov/=sw

        #print mn,cov[0,0]
        cov-=array([[v1*v2 for v1 in mn] for v2 in mn])
        #print mn,cov[0,0]
        #stop
        return mn, cov

    def plotAll (self,color, justlines=False, parlist=None):
            cc=0
            if not parlist:
                N=len(self.paramnames)
                parlist=range(N)
            else:
                N=len(parlist)
                for i,el in enumerate(parlist):
                    if type(el)==type('string'):
                        parlist[i]=self.parcol[el]-2
                    
            for ic,i in enumerate(parlist):
                for jc,j in enumerate(parlist):
                    cc+=1
                    if (ic<jc):
                        continue
                    if (i<0) or (j<0):
                        continue

                    pylab.subplot(N,N,cc)
                    if (ic==jc):
                        xv,yv=self.GetHisto(i+2)
                        pylab.plot(xv,yv,'-',color=color)
                    elif (ic>jc):
                        print i,j,'aaa',N
                        self.Plot2D(j+2,i+2,filled=color)
                        if (jc==0):
                            if not justlines:
                                pylab.ylabel(self.paramnames[i])
                            
                    if (ic==N-1) and (not justlines):
                        pylab.xlabel(self.paramnames[j])


def smline (x,y,mnmx):

    ##first lets.pad
    y=log(y+1e-30)
    N=len(y)
    y=array([y[0]]*N+list(y)+[y[-1]]*N)
    rft = fft.rfft(y)
    Nx=len(rft)
    k=linspace(0,1,Nx)
    rft*=exp(-k*k/(2*0.2**2))
    y = fft.irfft(rft)
    y=y[N:2*N]
    y=exp(y)
    return x,y
                       

def smline2(pic):
    Nx=len(pic[0])
    Ny=len(pic)
    picp=zeros((3*Nx,3*Ny))
    picp[Nx:2*Nx, Ny:2*Ny]=log(pic+1e-10)
    rft=fft.rfft2(picp)
    Nxf=len(rft[0])
    Nyf=len(rft)
    print Nx, Nxf
    print Ny, Nyf
    print rft.sum()
    for i in range(Nxf):
        for j in range(Nyf):
            kx=i*1.0/(Nxf)
            ky=j*1.0/(Nyf)
            k=sqrt(kx*kx+ky*ky)
            rft[j,i]*=exp(-k*k/(0.1*3.0**2))
    print rft.sum()
    picp=fft.irfft2(rft)
    pic=exp(picp[Nx:2*Nx, Ny:2*Ny])
    return pic
