#!/usr/bin/env python


from Simple_Plots import Simple_Plots


dir_name   = 'chains/'
roots      = ['wCDM_phy_BBAO+JLA', 'wCDM_phy_BBAO+Pantheon+Planck_15']
params     = ['h', 'w', 'Ol', 'Age']
param_pair = ['h', 'w']


S = Simple_Plots(dir_name, roots, params)
S.label = ['BBAO+JLA', 'BBAO+Pantheon+PLK15']
S.Plots1D()





"""
            #C = ChainIterator('chains_140411_155701',
            #                  'LCDM', 'phy', 'BBAO+CMBP')
            #grlist = []
            #wlist = []
            #for i in range(0, C.N, 1000):
            #    T = C.theory(i)
            #    grlist.append(T.growth(10.0))
            #    wlist.append(C.weight(i))

            #grlist = array(grlist)
            #wlist = array(wlist)
            #mn = (grlist*wlist).sum()/wlist.sum()
            #er = sqrt((grlist**2*wlist).sum()/wlist.sum()-mn**2)
            #print("growth z=10/z=0 = ", mn, "+/-", er)




# -----------------------
Plot_2D = 'False'
model_2D = ['LCDM_phy']
param_2D = ['Om', 'h']
NBins_2D = 40

xrange_2D = 'False'
xmin, xmax = 0.2, 0.4
ymin, ymax = 0.6, 0.75



if 'FALSE':
    a = 0
    for model in model_2D:
        for datasets in datasetl:
            a += 1

            C = cosmochain(dire + model+'_'+datasets)
            C.Plot2D(param_2D[0], param_2D[1], filled=colour(
                a), label=cosmodata(datasets), N=NBins_2D)

    if 'True' in xrange_2D:
        pylab.xlim(xmin, xmax)
        pylab.ylim(ymin, ymax)

    leg = pylab.legend()
    leg.draw_frame(False)		# No box & colour legend
    color_legend(leg)

    pylab.xlabel(C.latexname(param_2D[0]))
    pylab.ylabel(C.latexname(param_2D[1]))
    pylab.savefig(name_fig+'_2D.pdf')
    pylab.show()

else:
    print('Nothing else to do')

"""
