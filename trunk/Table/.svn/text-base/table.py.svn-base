import numpy as np
import string

def value(file,name):

    data=np.loadtxt(file+'.margestats', skiprows=3, comments='#', dtype='S') #usecols=[0,1,2],dtype='S')#{'names': ('param', 'avg', 'std'), 'formats': ('S10','f13','f13')})
    par=data[:,0]
    avg=data[:,1].astype(np.float)
    std=data[:,2].astype(np.float)

    for i in range(np.size(par)):
        if(par[i]==name):
            n=order(std[i])
            return "%.*f (%i) " %(n,avg[i],np.rint(10**n*std[i]))
    return "--"


#WE SHOULD RETURN A STRING, NOT A NUMBER!!
def order(std):
    for i in range(-8,0):
        x=np.int(std/10.**i)
#        print i,std,x
        if (x<25):
            return -i
    return 0
    print 'Number too large'


#data=np.loadtxt(file+'.params')
#if (string.find(name,data) == -1):
#        ltx='---'

#omh2_avg,omh2_std=value('lcdm_planck','omegamh2')
#om_avg,  om_std  =value('lcdm_planck','omegam')
#h0_avg,  h0_std  =value('lcdm_planck','H0')
#ok_avg,  ok_std  =value('lcdm_planck','omegak')
#w_avg,   w_std   =value('lcdm_planck','w')
#wa_avg,  wa_std  =value('lcdm_planck','wa')


outputfile = file('table.dat','w')

outputfile.write(r"\documentclass{article}"+'\n')
outputfile.write(r"\usepackage{graphicx}"+'\n')
outputfile.write(r"\usepackage{geometry}"+'\n')
outputfile.write(r"\usepackage{bm}"+'\n')
outputfile.write(r"\geometry{verbose,tmargin=90pt,bmargin=90pt,lmargin=90pt,rmargin=90pt}"+'\n')
outputfile.write(r"\begin{document}"+'\n')
outputfile.write(r"\clearpage"+'\n')
outputfile.write(r"\thispagestyle{empty}"+'\n')

#    outputfile.write(r"\begin{center}"+'\n')
#    outputfile.write(r"\begin{tabular}[t]{lr}"+'\n')
#    outputfile.write(r"\textbf{Model: %s} & \textbf{Data: %s} \\" %(cosmomodel,cosmodata) +'\n')
#    outputfile.write(r"\end{tabular}"+'\n')
#    outputfile.write(r"\end{center}"+'\n')
outputfile.write(r"\begin{table}[!h]"+'\n')
outputfile.write(r"\tiny"+'\n')
    #outputfile.write(r"\centering"+'\n')
    #outputfile.write(r"\resizebox{0.4\columnwidth}{!}{%"+'\n')
#    outputfile.write(r"\begin{minipage}[t]{0.5\linewidth}"+'\n')
outputfile.write(r"\begin{tabular}[t]{llcccccccccccc}"+'\n')
outputfile.write(r"\hline"+'\n')
outputfile.write(r"Model & Data & $\Omega_m$ & $\Omega_b h^2$ & $h$ & $\Omega_k$ & $w$ & $w_a$ & $m_{\nu}$ & $q$ & $z_a$ & $z_b$ & $\Omega_1$ & $\Omega_2$ \\"+'\n')

tablelist=np.loadtxt('input.txt',dtype='S')
for i in range(len(tablelist)):

    model=tablelist[i]
    if (string.find(model,'LadoCDM') != -1):
        cosmomodel='LadoCDM'
    elif (string.find(model,'JordiCDM') != -1):
        cosmomodel='JordiCDM'
    elif (string.find(model,'LCDMmasslessnu') != -1):
        cosmomodel='$\\bm{(\\Sigma\\nu=0)}$CDM'
    elif (string.find(model,'nuLCDM') != -1):
        cosmomodel='$\\bm{\\nu}$CDM'
    elif (string.find(model,'owaCDM') != -1):
        cosmomodel='$\\bm{ow_0w_a}$CDM'
    elif (string.find(model,'waCDM') != -1):
        cosmomodel='$\\bm{w_0w_a}$CDM'
    elif (string.find(model,'wCDM') != -1):
        cosmomodel='$\\bm{w}$CDM'
    elif (string.find(model,'oLCDM') != -1):
        cosmomodel='$\\bm{o}$CDM'
    elif (string.find(model,'LCDM') != -1):
        cosmomodel='$\\bm{\\Lambda}$CDM'
    else:
        print 'model not found!',model

    cosmodata=''
#    if (string.find(model,'CMBP') != -1):
#        cosmodata=' Planck'
#    elif (string.find(model,'CMBW') != -1):
#        cosmodata=' WMAP'

    if (string.find(model,'phy') != -1):
        cosmodata=cosmodata+' phy'
    elif (string.find(model,'pre') != -1):
        cosmodata=cosmodata+' pre'


    if (string.find(model,'CMBP') != -1):
        cosmodata=cosmodata+' Planck'
    elif (string.find(model,'CMBW') != -1):
        cosmodata=cosmodata+' WMAP'


    if (string.find(model,'BBAO') != -1):
        cosmodata=cosmodata+' AllBAO'
    if (string.find(model,'GBAO') != -1):
        cosmodata=cosmodata+' GalBAO'
    if (string.find(model,'LBAO') != -1):
        cosmodata=cosmodata+' LyaBAO'
    elif (string.find(model,'SN') != -1):
        cosmodata=cosmodata+' SN'
    elif (string.find(model,'RiessH0') != -1):
        cosmodata=cosmodata+' H0'


    str_om=value(tablelist[i],'Om')
    str_ob=value(tablelist[i],'Obh2')
    str_h =value(tablelist[i],'h')
    str_ok=value(tablelist[i],'Ok')
    str_wa=value(tablelist[i],'wa')
    str_w =value(tablelist[i],'w')
    str_mnu=value(tablelist[i],'mnu')
    str_q=value(tablelist[i],'q')
    str_za=value(tablelist[i],'za')
    str_zb=value(tablelist[i],'zb')
    str_om1=value(tablelist[i],'Om1')
    str_om2=value(tablelist[i],'Om2')
    outputfile.write(r"%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\" %(cosmomodel,cosmodata,str_om,str_ob,str_h,str_ok,str_w,str_wa,str_mnu,str_q,str_za,str_zb,str_om1,str_om2) +'\n')

#    outputfile.write(r"\multicolumn{2}{c}{Model: %s, Data: %s} \\" %(cosmomodel,cosmodata) +'\n')
#    outputfile.write(r"\hline"+'\n')
outputfile.write(r"\end{tabular}"+'\n')
#    outputfile.write(r"\end{minipage}"+'\n')
#    outputfile.write(r"}"+'\n')
outputfile.write(r"\end{table}"+'\n')
outputfile.write(r"\end{document}"+'\n')

#tablelist=np.loadtxt('input.txt',dtype='S')
#for i in range(len(tablelist)):
#    do_table(tablelist[i])


