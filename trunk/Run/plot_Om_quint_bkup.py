
#!/usr/bin/env python
import pylab
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os, sys, time

dir = '/gpfs01/astro/workarea/jvazquez/cosmomc_Quint/camb/'

params1 = {'backend': 'pdf',
               'axes.labelsize': 16,
               'text.fontsize': 18,
               'xtick.labelsize': 20,
               'ytick.labelsize': 20,
               'legend.draw_frame': False,
               'legend.fontsize': 12,
               'lines.markersize': 6,
               'font.size': 20,
               'text.usetex': True}#
pylab.rcParams.update(params1)


def colour(x):
    if x==1: return 'red'
    if x==2: return 'blue'
    if x==3: return 'black'
    if x==4: return 'magenta'
    if x==5: return 'cyan'
    if x==6: return 'orange'
    if x==7: return 'green'
    if x==8: return 'yellow'
    if x==9: return 'purple'
    if x>9: return 'black' # print("Increased colouring")

root = 'test_'

names = ['test_ede_Om_100.dat'] #'Om_Quint_B1','Om_Quint_B2','Om_Quint_B3','Om_Quint_B4','Om_Quint_B5']
rs = ['144.910']#'135.47','140.47','143.12','144.72','146.90']
rs_lcdm = 147.78 #149.01
s=0
lnat, Omt, Oedet, Oedet0, inter_Oede = [], [], [], [], []
for name in names:
  pnames=open(dir+name+'.dat').readlines()
  lpnames = len(pnames)
  lna, Om, Oede, Oede0 = [], [], [], []  
  for l in range(lpnames):
     if ((l< lpnames/3.) or (l> lpnames/3. and l%50 == 0)):      
       vals =pnames[l].split()[0:]
       lna.append(float(vals[0])) 
       Oede0.append(float(vals[1]))	
       Oede.append((1. - float(vals[1]))**(-0.5))
       Om.append(float(rs[s])*(1. - float(vals[1]))**(-0.5)/rs_lcdm)     
  lnat.append(lna)
  Omt.append(Om)
  Oedet.append(Oede)
  Oedet0.append(Oede0)
  s+=1

inter_Oede0_0 = interp1d(lnat[0], Oedet0[0])
#inter_Oede0_1 = interp1d(lnat[1], Oedet0[1])
#inter_Oede0_2 = interp1d(lnat[2], Oedet0[2])
#inter_Oede0_3 = interp1d(lnat[3], Oedet0[3])
#inter_Oede0_4 = interp1d(lnat[4], Oedet0[4])

inter_Oede_0 = interp1d(lnat[0], Oedet[0])
#inter_Oede_1 = interp1d(lnat[1], Oedet[1])
#inter_Oede_2 = interp1d(lnat[2], Oedet[2])
#inter_Oede_3 = interp1d(lnat[3], Oedet[3])
#inter_Oede_4 = interp1d(lnat[4], Oedet[4])


namesq = ['test_ede_Da_0100.dat']#'Quint_B1','Quint_B2','Quint_B3','Quint_B4','Quint_B5']
lnaat, Dat, Dht = [], [], []
pnames=open(dir+'test_Da_ede_lcdm'+'.dat').readlines()
lpnames = len(pnames)
Da_lcdm, Dh_lcdm = [], []

for l in range(lpnames): 
   vals =pnames[l].split()[0:]
   Da_lcdm.append(float(vals[1]))
   Dh_lcdm.append(float(vals[2]))



s=0
for name in namesq:
    pnames=open(dir+name+'.dat').readlines()
    lpnames = len(pnames)
    lnaa, Da, Dh = [], [], []
    print lpnames, len(Da_lcdm)
    for l in range(lpnames):
       vals =pnames[l].split()[0:]
       lnaa.append(float(vals[0]))  
       Da.append((float(vals[1]))/(Da_lcdm[l])) 
       Dh.append((float(vals[2]))/(Dh_lcdm[l]))
    s+=1       

    lnaat.append(lnaa)
    Dat.append(Da) 
    Dht.append(Dh)

inter_Da_0 = interp1d(lnaat[0], Dat[0])
#inter_Da_1 = interp1d(lnaat[1], Dat[1])
#inter_Da_2 = interp1d(lnaat[2], Dat[2])
#inter_Da_3 = interp1d(lnaat[3], Dat[3])
#inter_Da_4 = interp1d(lnaat[4], Dat[4])

inter_Dh_0 = interp1d(lnaat[0], Dht[0])
#inter_Dh_1 = interp1d(lnaat[1], Dht[1])
#inter_Dh_2 = interp1d(lnaat[2], Dht[2])
#inter_Dh_3 = interp1d(lnaat[3], Dht[3])
#inter_Dh_4 = interp1d(lnaat[4], Dht[4])

ar = np.log(1./(1+1060.0))
nn0, nn1, nn2, nn3, nn4, nn5 = [], [], [], [], [], []
for i in range(0,len(names)):
   n0 = 'inter_Oede0_%i(ar)'%(i)
   n1 = 'inter_Oede_%i(ar)'%(i) 
   n2 = 'inter_Da_%i(ar)'%(i)
   n3 = 'inter_Dh_%i(ar)'%(i)
   n4 = 'inter_Da_%i(ar)'%(i)
   n5 = '(1+0.41*inter_Oede0_%i(ar))/(1-0.1*inter_Oede0_%i(ar))'%(i, i)

   nn0.append(eval(n0))
   nn1.append(eval(n1)*float(rs[i])/rs_lcdm)
   nn2.append(eval(n1)*eval(n2))
   nn3.append(eval(n1)*eval(n3))
   nn4.append(eval(n4)*rs_lcdm/float(rs[i]))
   nn5.append(eval(n5))
  
#-------------------------------------

#if False:
# fig =pylab.figure(figsize=(16,15))

# ax = fig.add_subplot(3,3,1)
# for i in range(0,len(names)):
#     ax.plot(lnat[i], Oedet0[i], color=colour(i+1), label = "$r_s$ =%3.2f,"%(float(rs[i])))
# ax.plot([ar,ar],[0,0.8])
# ax.grid(True)
# plt.legend(loc="upper left")
# plt.ylim([0,0.8])
# plt.xlim([-14,0])
# plt.xlabel("$\ln a$")
# plt.ylabel("$\Omega_{ede}$")

# ax2 = fig.add_subplot(3,3,4)
# for i in range(0,len(names)):
#    ax2.plot(lnat[i],Omt[i],color=colour(i+1))
# ax2.grid(True)
# plt.xlim([-14,0])
# #plt.ylim(ymax=1.2)
# plt.ylim([0.9,1.2])
# plt.xlabel("$\ln a$")
# plt.ylabel("$[r_{s,ede}~/~\sqrt{1-\Omega_{ede}}]~/~r_{s,LCDM}$")

# ax3 = fig.add_subplot(3,3,5)
# for i in range(0,len(names)):
#     nameO = 'inter_Oede_%i(lnaat[%i])'%(i,i)
#     ax3.plot(lnaat[i], np.array(Dat[i])*np.array(eval(nameO)),color=colour(i+1))   
# ax3.grid(True)
# plt.ylim([0.95,1.2])
# plt.xlabel("$\ln a$")
# plt.ylabel("$[D_{A,ede}~/~\sqrt{1-\Omega_{ede}}]~/~D_{A,LCDM}$") 

# ax4 = fig.add_subplot(3,3,7)
# ax4.plot(nn0, nn1)
# ax4.grid(True)
# plt.xlabel("$\Omega_{ede}(z_{drag})$")
# plt.ylabel("$[r_{s,ede}~/~\sqrt{1-\Omega_{ede}}]~/~r_{s,LCDM} (z_{drag})$")

# ax5 = fig.add_subplot(3,3,8)
# ax5.plot(nn0, nn2)
# ax5.grid(True)
# plt.xlabel("$\Omega_{ede}(z_{drag})$")
# plt.ylabel("$[D_{A,ede}~/~\sqrt{1-\Omega_{ede}}]~/~D_{A,LCDM} (z_{drag})$")

# ax6 = fig.add_subplot(3,3,9)
# ax6.plot(nn0, nn3)
# ax6.grid(True)
# plt.xlabel("$\Omega_{ede}(z_{drag})$")
# plt.ylabel("$[D_{H,ede}~/~\sqrt{1-\Omega_{ede}}]~/~D_{H,LCDM} (z_{drag})$")

# ax7 = fig.add_subplot(3,3,6)
# for i in range(0,len(namesq)):
#    ax7.plot(lnaat[i],np.array(Dht[i])*(rs_lcdm/float(rs[i])),color=colour(i+1))
##        label = "$r_s$ =%3.2f,"%(float(rs[i])))
## plt.legend(loc="lower right")
# ax7.grid(True)
 #plt.xlabel("$\ln a$")
 #plt.ylabel("$[D_H/r_s]_{ede}~/~[D_H/r_s]_{LCDM}$")


 #plt.tight_layout()
 #plt.savefig('Dh_Da'+".pdf")
 #plt.show()


#fig =pylab.figure(figsize=(6,5))
#ax = fig.add_subplot(1,1,1)
#ax.plot(nn0, nn4, label = '$[D_A/r_s]~/~[D_A/r_s]_{lcdm}$')
#ax.plot(nn0, nn5, label = '$\sim \sqrt{1+\Omega_{ede}}$')
#plt.xlabel("$\Omega_{ede}(z_{drag})$")
#plt.legend(loc="upper left")
##plt.ylabel("$[D_A/r_s]~/~[D_A/r_s]_{lcdm}$") #/~\sqrt{1+\Omega_{ede}}$")
#plt.tight_layout()
#plt.savefig('Final'+".pdf")
#plt.show()


#if False: 
# ax2 = fig.add_subplot(1,3,2)
# for i in range(0,len(namesq)):
#    ax2.plot(lnaat[i],np.array(Dat[i])*(rs_lcdm/float(rs[i])),color=colour(i+1))
# ax2.grid(True)
# plt.xlabel("$\ln a$")
# plt.ylabel("$[D_A/r_s]_{ede}~/~[D_A/r_s]_{LCDM}$")


# ax3 = fig.add_subplot(1,3,3)
# for i in range(0,len(namesq)):
#    ax3.plot(lnaat[i],np.array(Dht[i])*(rs_lcdm/float(rs[i])),color=colour(i+1),
#	label = "$r_s$ =%3.2f,"%(float(rs[i])))
# plt.legend(loc="lower right")
# ax3.grid(True)
# plt.xlabel("$\ln a$")
# plt.ylabel("$[D_H/r_s]_{ede}~/~[D_H/r_s]_{LCDM}$")


# plt.tight_layout()
# plt.savefig('Dh_Da'+".pdf")
# plt.show()

