from simplemc.CosmoCalc import CosmoCalc

#Computes and plot (optional)
#['HIOverrd', 'DVOverrd', 'DaOverrd', 'Hubble', 'HubInvOverz', 'SNIa', 'fs8']


#C = CosmoCalc('LCDM', 'Hubble')
#C.run_plot(lw='2')

#C = CosmoCalc('LCDM', 'Age')
#C.Age()

#C = CosmoCalc('LCDM', 'DaOverrd', 'h', 0.4, 0.9)
#C.run_plot(lw='2')


#C = CosmoCalc('owaCDM', 'fs8', 'wa', -0.5, 0.5, 5, zmax=0.8)
#C.run_plot(lw='1')


C = CosmoCalc('owaCDM', 'fs8', 'wa', -0.5, 0.5, 5, zmax=3.1, plot_data=True)
C.run_plot(lw='1')

