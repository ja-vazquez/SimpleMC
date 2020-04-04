#!/usr/bin/env python


from Simple_Plots import Simple_Plots


dir_name   = 'chains/'
roots      = ['wCDM_phy_BBAO+JLA'] #, 'wCDM_phy_BBAO+Pantheon+Planck_15']
params_1D  = ['h', 'w', 'Ol', 'Age']
params_2D  = [['h', 'w'], ['Om', 'h']]


S = Simple_Plots(dir_name, roots)
S.label = ['BBAO+JLA', 'BBAO+Pantheon+PLK15']

S.Show_limits(params_1D)
S.Covariance(params_1D)
S.Plots1D(params_1D)
S.Plots2D(params_2D)
S.plotAlls(params_1D)


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

"""

