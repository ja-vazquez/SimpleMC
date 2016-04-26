
#!/usr/bin/env python
import pylab
import matplotlib.pyplot as plt

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

lnat, Omt = [], []
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
#names = ['Om_Quint','Om_Quint2','Om_Quint3','Om_Quint4','Om_Quint5','Om_Quint6']
#rs = ['143.12','143.15','143.16','143.26','143.22','143.30']

names = ['Om_Quint_B','Om_Quint_B2','Om_Quint_B3','Om_Quint_B4','Om_Quint_B5']
rs = ['135.47','135.38','135.43','135.48','135.43']

fig =pylab.figure(figsize=(8,6))
ax = fig.add_subplot(1,1,1)
s=0
for name in names:
  pnames=open(dir+root+name+'.dat').readlines()
  lpnames = len(pnames)
  lna, Om = [], []  
  for l in range(lpnames):
     if ((l< lpnames/3.) or (l> lpnames/3. and l%50 == 0)):      
       vals =pnames[l].split()[0:]
       lna.append(float(vals[0])) 
       Om.append(float(rs[s])*(1. - float(vals[1]))**(-0.5)/149.01)      
  lnat.append(lna)
  Omt.append(Om)
  s+=1

for i in range(0,len(names)):
    ax.plot(lnat[i],Omt[i],color=colour(i+1))
#ax.plot([-9.,-9.],[0.9,1.2],'k-')
ax.grid(True)
plt.xlim([-14,0])
#plt.ylim(ymax=1.2)
plt.ylim([0.9,1.2])
plt.xlabel("$\ln a$")
plt.ylabel("$r_{s,EDE}/\sqrt{1-\Omega_{EDE}}/r_{s,LCDM}$")

plt.show()
