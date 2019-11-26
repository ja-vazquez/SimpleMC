import numpy as np

#outputfile = file('table.dat','w')
# for l in range(221):
#   outputfile.write("%f %f %f %f %f\n" %(xint[l],zint_low2[l],zint_low1[l],zint_upp1[l],zint_upp2[l]))
# outputfile.close()


def value(file, name):

    # {'names': ('param', 'avg', 'std'), 'formats': ('S10','f13','f13')})
    data = np.loadtxt(file+'.margestats', skiprows=3,
                      comments='#', usecols=[0, 1, 2], dtype='S')
    par = data[:, 0]
    avg = data[:, 1].astype(np.float)
    std = data[:, 2].astype(np.float)

    for i in range(np.size(par)):
        if(par[i] == name):
            return avg[i], std[i]


outputfile = file('table15.dat', 'w')

outputfile.write(r"\begin{table*}"+'\n')
outputfile.write(r"\centering"+'\n')
outputfile.write(r"\begin{tabular}{llllllll}"+'\n')
outputfile.write(r"\hline"+'\n')
outputfile.write(
    r"Cosmological & Data Sets & $\Omega_{\rm m} h^{2}$ & $\Omega_{\rm m}$ & $H_{0}$ & $\Omega_{\rm K}$ & $w_{0}$ & $w_{a}$ \\"+'\n')
outputfile.write(r"Model & & & & km s$^{-1}$ Mpc$^{-1}$ & & & \\"+'\n')
#outputfile.write(r"Model & & & & & & & \\"+'\n')
outputfile.write(r"\hline"+'\n')


# LCDM
omh2_avg, omh2_std = value('lcdm_planck_bossiso', 'omegamh2')
om_avg,  om_std = value('lcdm_planck_bossiso', 'omegam')
h0_avg,  h0_std = value('lcdm_planck_bossiso', 'H0')
outputfile.write(r"$\Lambda$CDM & Planck + CMASS-iso + LOWZ & %.4f (%i) & %.3f (%i) & %.1f (%i) & \nodata & \nodata & \nodata \\" %
                 (omh2_avg, np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std)) + '\n')
omh2_avg, omh2_std = value('lcdm_planck_boss', 'omegamh2')
om_avg,  om_std = value('lcdm_planck_boss', 'omegam')
h0_avg,  h0_std = value('lcdm_planck_boss', 'H0')
outputfile.write(r"$\Lambda$CDM & Planck + CMASS + LOWZ & %.4f (%i) & %.3f (%i) & %.1f (%i) & \nodata & \nodata & \nodata \\" %
                 (omh2_avg, np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std)) + '\n')
omh2_avg, omh2_std = value('lcdm_planck_bao', 'omegamh2')
om_avg,  om_std = value('lcdm_planck_bao', 'omegam')
h0_avg,  h0_std = value('lcdm_planck_bao', 'H0')
outputfile.write(r"$\Lambda$CDM & Planck + BAO & %.4f (%i) & %.3f (%i) & %.1f (%i) & \nodata & \nodata & \nodata \\" %
                 (omh2_avg, np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std)) + '\n')
omh2_avg, omh2_std = value('lcdm_planck_boss_snu', 'omegamh2')
om_avg,  om_std = value('lcdm_planck_boss_snu', 'omegam')
h0_avg,  h0_std = value('lcdm_planck_boss_snu', 'H0')
outputfile.write(r"$\Lambda$CDM & Planck + CMASS + LOWZ + SN & %.4f (%i) & %.3f (%i) & %.1f (%i) & \nodata & \nodata & \nodata \\" %
                 (omh2_avg, np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std)) + '\n')
omh2_avg, omh2_std = value('lcdm_planck_bao_snu', 'omegamh2')
om_avg,  om_std = value('lcdm_planck_bao_snu', 'omegam')
h0_avg,  h0_std = value('lcdm_planck_bao_snu', 'H0')
outputfile.write(r"$\Lambda$CDM & Planck + BAO + SN & %.4f (%i) & %.3f (%i) & %.1f (%i) & \nodata & \nodata & \nodata \\" %
                 (omh2_avg, np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std)) + '\n')
omh2_avg, omh2_std = value('lcdm_wmap_bao_snu', 'omegamh2')
om_avg,  om_std = value('lcdm_wmap_bao_snu', 'omegam')
h0_avg,  h0_std = value('lcdm_wmap_bao_snu', 'H0')
outputfile.write(r"$\Lambda$CDM & WMAP + BAO + SN & %.4f (%i) & %.3f (%i) & %.1f (%i) & \nodata & \nodata & \nodata \\" %
                 (omh2_avg, np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std)) + '\n')
omh2_avg, omh2_std = value('lcdm_ewmap_bao_snu', 'omegamh2')
om_avg,  om_std = value('lcdm_ewmap_bao_snu', 'omegam')
h0_avg,  h0_std = value('lcdm_ewmap_bao_snu', 'H0')
outputfile.write(r"$\Lambda$CDM & \textit{e}WMAP + BAO + SN & %.4f (%i) & %.3f (%i) & %.1f (%i) & \nodata & \nodata & \nodata \\" % (
    omh2_avg, np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std)) + '\n')
outputfile.write(r"\hline"+'\n')

# oCDM
omh2_avg, omh2_std = value('ocdm_planck_bossiso', 'omegamh2')
om_avg,  om_std = value('ocdm_planck_bossiso', 'omegam')
h0_avg,  h0_std = value('ocdm_planck_bossiso', 'H0')
ok_avg,  ok_std = value('ocdm_planck_bossiso', 'omegak')
outputfile.write(r"oCDM & Planck + CMASS-iso + LOWZ & %.4f (%i) & %.3f (%i) & %.1f (%i) & %+.4f (%i) & \nodata & \nodata \\" %
                 (omh2_avg, np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), ok_avg, np.int(10000*ok_std)) + '\n')
omh2_avg, omh2_std = value('ocdm_planck_boss', 'omegamh2')
om_avg,  om_std = value('ocdm_planck_boss', 'omegam')
h0_avg,  h0_std = value('ocdm_planck_boss', 'H0')
ok_avg,  ok_std = value('ocdm_planck_boss', 'omegak')
outputfile.write(r"oCDM & Planck + CMASS + LOWZ & %.4f (%i) & %.3f (%i) & %.1f (%i) & %+.4f (%i) & \nodata & \nodata \\" %
                 (omh2_avg, np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), ok_avg, np.int(10000*ok_std)) + '\n')
omh2_avg, omh2_std = value('ocdm_planck_bao', 'omegamh2')
om_avg,  om_std = value('ocdm_planck_bao', 'omegam')
h0_avg,  h0_std = value('ocdm_planck_bao', 'H0')
ok_avg,  ok_std = value('ocdm_planck_bao', 'omegak')
outputfile.write(r"oCDM & Planck + BAO & %.4f (%i) & %.3f (%i) & %.1f (%i) & %+.4f (%i) & \nodata & \nodata \\" %
                 (omh2_avg, np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), ok_avg, np.int(10000*ok_std)) + '\n')
omh2_avg, omh2_std = value('ocdm_planck_boss_snu', 'omegamh2')
om_avg,  om_std = value('ocdm_planck_boss_snu', 'omegam')
h0_avg,  h0_std = value('ocdm_planck_boss_snu', 'H0')
ok_avg,  ok_std = value('ocdm_planck_boss_snu', 'omegak')
outputfile.write(r"oCDM & Planck + CMASS + LOWZ + SN & %.4f (%i) & %.3f (%i) & %.1f (%i) & %+.4f (%i) & \nodata & \nodata \\" %
                 (omh2_avg, np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), ok_avg, np.int(10000*ok_std)) + '\n')
omh2_avg, omh2_std = value('ocdm_planck_bao_snu', 'omegamh2')
om_avg,  om_std = value('ocdm_planck_bao_snu', 'omegam')
h0_avg,  h0_std = value('ocdm_planck_bao_snu', 'H0')
ok_avg,  ok_std = value('ocdm_planck_bao_snu', 'omegak')
outputfile.write(r"oCDM & Planck + BAO + SN & %.4f (%i) & %.3f (%i) & %.1f (%i) & %+.4f (%i) & \nodata & \nodata \\" %
                 (omh2_avg, np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), ok_avg, np.int(10000*ok_std)) + '\n')
omh2_avg, omh2_std = value('ocdm_wmap_bao_snu', 'omegamh2')
om_avg,  om_std = value('ocdm_wmap_bao_snu', 'omegam')
h0_avg,  h0_std = value('ocdm_wmap_bao_snu', 'H0')
ok_avg,  ok_std = value('ocdm_wmap_bao_snu', 'omegak')
outputfile.write(r"oCDM & WMAP + BAO + SN & %.4f (%i) & %.3f (%i) & %.1f (%i) & %+.4f (%i) & \nodata & \nodata \\" %
                 (omh2_avg, np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), ok_avg, np.int(10000*ok_std)) + '\n')
omh2_avg, omh2_std = value('ocdm_ewmap_bao_snu', 'omegamh2')
om_avg,  om_std = value('ocdm_ewmap_bao_snu', 'omegam')
h0_avg,  h0_std = value('ocdm_ewmap_bao_snu', 'H0')
ok_avg,  ok_std = value('ocdm_ewmap_bao_snu', 'omegak')
outputfile.write(r"oCDM & \textit{e}WMAP + BAO + SN & %.4f (%i) & %.3f (%i) & %.1f (%i) & %+.4f (%i) & \nodata & \nodata \\" % (
    omh2_avg, np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), ok_avg, np.int(10000*ok_std)) + '\n')
outputfile.write(r"\hline"+'\n')

# wCDM
omh2_avg, omh2_std = value('wcdm_planck_bossiso', 'omegamh2')
om_avg,  om_std = value('wcdm_planck_bossiso', 'omegam')
h0_avg,  h0_std = value('wcdm_planck_bossiso', 'H0')
w_avg,   w_std = value('wcdm_planck_bossiso', 'w')
outputfile.write(r"$w$CDM & Planck + CMASS-iso + LOWZ & %.4f (%i) & %.3f (%i) & %.1f (%i) & \nodata & %.2f (%i) & \nodata \\" %
                 (omh2_avg, np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), w_avg, np.int(100*w_std)) + '\n')
omh2_avg, omh2_std = value('wcdm_planck_boss', 'omegamh2')
om_avg,  om_std = value('wcdm_planck_boss', 'omegam')
h0_avg,  h0_std = value('wcdm_planck_boss', 'H0')
w_avg,   w_std = value('wcdm_planck_boss', 'w')
outputfile.write(r"$w$CDM & Planck + CMASS + LOWZ & %.4f (%i) & %.3f (%i) & %.1f (%i) & \nodata & %.2f (%i) & \nodata \\" %
                 (omh2_avg, np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), w_avg, np.int(100*w_std)) + '\n')
omh2_avg, omh2_std = value('wcdm_planck_bao', 'omegamh2')
om_avg,  om_std = value('wcdm_planck_bao', 'omegam')
h0_avg,  h0_std = value('wcdm_planck_bao', 'H0')
w_avg,   w_std = value('wcdm_planck_bao', 'w')
outputfile.write(r"$w$CDM & Planck + BAO & %.4f (%i) & %.3f (%i) & %.1f (%i) & \nodata & %.2f (%i) & \nodata \\" %
                 (omh2_avg, np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), w_avg, np.int(100*w_std)) + '\n')
omh2_avg, omh2_std = value('wcdm_planck_boss_snu', 'omegamh2')
om_avg,  om_std = value('wcdm_planck_boss_snu', 'omegam')
h0_avg,  h0_std = value('wcdm_planck_boss_snu', 'H0')
w_avg,   w_std = value('wcdm_planck_boss_snu', 'w')
outputfile.write(r"$w$CDM & Planck + CMASS + LOWZ + SN & %.4f (%i) & %.3f (%i) & %.1f (%i) & \nodata & %.2f (%i) & \nodata \\" %
                 (omh2_avg, np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), w_avg, np.int(100*w_std)) + '\n')
omh2_avg, omh2_std = value('wcdm_planck_bao_snu', 'omegamh2')
om_avg,  om_std = value('wcdm_planck_bao_snu', 'omegam')
h0_avg,  h0_std = value('wcdm_planck_bao_snu', 'H0')
w_avg,   w_std = value('wcdm_planck_bao_snu', 'w')
outputfile.write(r"$w$CDM & Planck + BAO + SN & %.4f (%i) & %.3f (%i) & %.1f (%i) & \nodata & %.2f (%i) & \nodata \\" %
                 (omh2_avg, np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), w_avg, np.int(100*w_std)) + '\n')
omh2_avg, omh2_std = value('wcdm_wmap_bao_snu', 'omegamh2')
om_avg,  om_std = value('wcdm_wmap_bao_snu', 'omegam')
h0_avg,  h0_std = value('wcdm_wmap_bao_snu', 'H0')
w_avg,   w_std = value('wcdm_wmap_bao_snu', 'w')
outputfile.write(r"$w$CDM & WMAP + BAO + SN & %.4f (%i) & %.3f (%i) & %.1f (%i) & \nodata & %.2f (%i) & \nodata \\" %
                 (omh2_avg, np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), w_avg, np.int(100*w_std)) + '\n')
omh2_avg, omh2_std = value('wcdm_ewmap_bao_snu', 'omegamh2')
om_avg,  om_std = value('wcdm_ewmap_bao_snu', 'omegam')
h0_avg,  h0_std = value('wcdm_ewmap_bao_snu', 'H0')
w_avg,   w_std = value('wcdm_ewmap_bao_snu', 'w')
outputfile.write(r"$w$CDM & \textit{e}WMAP + BAO + SN & %.4f (%i) & %.3f (%i) & %.1f (%i) & \nodata & %.2f (%i) & \nodata \\" % (
    omh2_avg, np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), w_avg, np.int(100*w_std)) + '\n')
outputfile.write(r"\hline"+'\n')

# owCDM
omh2_avg, omh2_std = value('owcdm_planck_bossiso', 'omegamh2')
om_avg,  om_std = value('owcdm_planck_bossiso', 'omegam')
h0_avg,  h0_std = value('owcdm_planck_bossiso', 'H0')
ok_avg,  ok_std = value('owcdm_planck_bossiso', 'omegak')
w_avg,   w_std = value('owcdm_planck_bossiso', 'w')
outputfile.write(r"o$w$CDM & Planck + CMASS-iso + LOWZ & %.4f (%i) & %.3f (%i) & %.1f (%i) & %+.4f (%i) & %.2f (%i) & \nodata \\" % (omh2_avg,
                                                                                                                                     np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), ok_avg, np.int(10000*ok_std), w_avg, np.int(100*w_std)) + '\n')
omh2_avg, omh2_std = value('owcdm_planck_boss', 'omegamh2')
om_avg,  om_std = value('owcdm_planck_boss', 'omegam')
h0_avg,  h0_std = value('owcdm_planck_boss', 'H0')
ok_avg,  ok_std = value('owcdm_planck_boss', 'omegak')
w_avg,   w_std = value('owcdm_planck_boss', 'w')
outputfile.write(r"o$w$CDM & Planck + CMASS + LOWZ & %.4f (%i) & %.3f (%i) & %.1f (%i) & %+.4f (%i) & %.2f (%i) & \nodata \\" % (omh2_avg,
                                                                                                                                 np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), ok_avg, np.int(10000*ok_std), w_avg, np.int(100*w_std)) + '\n')
omh2_avg, omh2_std = value('owcdm_planck_bao', 'omegamh2')
om_avg,  om_std = value('owcdm_planck_bao', 'omegam')
h0_avg,  h0_std = value('owcdm_planck_bao', 'H0')
ok_avg,  ok_std = value('owcdm_planck_bao', 'omegak')
w_avg,   w_std = value('owcdm_planck_bao', 'w')
outputfile.write(r"o$w$CDM & Planck + BAO & %.4f (%i) & %.3f (%i) & %.1f (%i) & %+.4f (%i) & %.2f (%i) & \nodata \\" % (omh2_avg, np.int(
    10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), ok_avg, np.int(10000*ok_std), w_avg, np.int(100*w_std)) + '\n')
omh2_avg, omh2_std = value('owcdm_planck_boss_snu', 'omegamh2')
om_avg,  om_std = value('owcdm_planck_boss_snu', 'omegam')
h0_avg,  h0_std = value('owcdm_planck_boss_snu', 'H0')
ok_avg,  ok_std = value('owcdm_planck_boss_snu', 'omegak')
w_avg,   w_std = value('owcdm_planck_boss_snu', 'w')
outputfile.write(r"o$w$CDM & Planck + CMASS + LOWZ + SN & %.4f (%i) & %.3f (%i) & %.1f (%i) & %+.4f (%i) & %.2f (%i) & \nodata \\" % (omh2_avg,
                                                                                                                                      np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), ok_avg, np.int(10000*ok_std), w_avg, np.int(100*w_std)) + '\n')
omh2_avg, omh2_std = value('owcdm_planck_bao_snu', 'omegamh2')
om_avg,  om_std = value('owcdm_planck_bao_snu', 'omegam')
h0_avg,  h0_std = value('owcdm_planck_bao_snu', 'H0')
ok_avg,  ok_std = value('owcdm_planck_bao_snu', 'omegak')
w_avg,   w_std = value('owcdm_planck_bao_snu', 'w')
outputfile.write(r"o$w$CDM & Planck + BAO + SN & %.4f (%i) & %.3f (%i) & %.1f (%i) & %+.4f (%i) & %.2f (%i) & \nodata \\" % (omh2_avg,
                                                                                                                             np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), ok_avg, np.int(10000*ok_std), w_avg, np.int(100*w_std)) + '\n')
omh2_avg, omh2_std = value('owcdm_wmap_bao_snu', 'omegamh2')
om_avg,  om_std = value('owcdm_wmap_bao_snu', 'omegam')
h0_avg,  h0_std = value('owcdm_wmap_bao_snu', 'H0')
ok_avg,  ok_std = value('owcdm_wmap_bao_snu', 'omegak')
w_avg,   w_std = value('owcdm_wmap_bao_snu', 'w')
outputfile.write(r"o$w$CDM & WMAP + BAO + SN & %.4f (%i) & %.3f (%i) & %.1f (%i) & %+.4f (%i) & %.2f (%i) & \nodata \\" % (omh2_avg, np.int(
    10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), ok_avg, np.int(10000*ok_std), w_avg, np.int(100*w_std)) + '\n')
omh2_avg, omh2_std = value('owcdm_ewmap_bao_snu', 'omegamh2')
om_avg,  om_std = value('owcdm_ewmap_bao_snu', 'omegam')
h0_avg,  h0_std = value('owcdm_ewmap_bao_snu', 'H0')
ok_avg,  ok_std = value('owcdm_ewmap_bao_snu', 'omegak')
w_avg,   w_std = value('owcdm_ewmap_bao_snu', 'w')
outputfile.write(r"o$w$CDM & \textit{e}WMAP + BAO + SN & %.4f (%i) & %.3f (%i) & %.1f (%i) & %+.4f (%i) & %.2f (%i) & \nodata \\" % (omh2_avg, np.int(
    10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), ok_avg, np.int(10000*ok_std), w_avg, np.int(100*w_std)) + '\n')
outputfile.write(r"\hline"+'\n')

# w0waCDM
omh2_avg, omh2_std = value('w0wacdm_planck_bossiso', 'omegamh2')
om_avg,  om_std = value('w0wacdm_planck_bossiso', 'omegam')
h0_avg,  h0_std = value('w0wacdm_planck_bossiso', 'H0')
w_avg,   w_std = value('w0wacdm_planck_bossiso', 'w')
wa_avg,  wa_std = value('w0wacdm_planck_bossiso', 'wa')
outputfile.write(r"$w_0w_a$CDM & Planck + CMASS-iso + LOWZ & %.4f (%i) & %.3f (%i) & %.1f (%i) & \nodata & %.2f (%i) & %.2f (%i) \\" % (omh2_avg,
                                                                                                                                        np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), w_avg, np.int(100*w_std), wa_avg, np.int(100*wa_std)) + '\n')
omh2_avg, omh2_std = value('w0wacdm_planck_boss', 'omegamh2')
om_avg,  om_std = value('w0wacdm_planck_boss', 'omegam')
h0_avg,  h0_std = value('w0wacdm_planck_boss', 'H0')
w_avg,   w_std = value('w0wacdm_planck_boss', 'w')
wa_avg,  wa_std = value('w0wacdm_planck_boss', 'wa')
outputfile.write(r"$w_0w_a$CDM & Planck + CMASS + LOWZ & %.4f (%i) & %.3f (%i) & %.1f (%i) & \nodata & %.2f (%i) & %.2f (%i) \\" % (omh2_avg,
                                                                                                                                    np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), w_avg, np.int(100*w_std), wa_avg, np.int(100*wa_std)) + '\n')
omh2_avg, omh2_std = value('w0wacdm_planck_bao', 'omegamh2')
om_avg,  om_std = value('w0wacdm_planck_bao', 'omegam')
h0_avg,  h0_std = value('w0wacdm_planck_bao', 'H0')
w_avg,   w_std = value('w0wacdm_planck_bao', 'w')
wa_avg,  wa_std = value('w0wacdm_planck_bao', 'wa')
outputfile.write(r"$w_0w_a$CDM & Planck + BAO & %.4f (%i) & %.3f (%i) & %.1f (%i) & \nodata & %.2f (%i) & %.2f (%i) \\" % (omh2_avg,
                                                                                                                           np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), w_avg, np.int(100*w_std), wa_avg, np.int(100*wa_std)) + '\n')
omh2_avg, omh2_std = value('w0wacdm_planck_boss_snu', 'omegamh2')
om_avg,  om_std = value('w0wacdm_planck_boss_snu', 'omegam')
h0_avg,  h0_std = value('w0wacdm_planck_boss_snu', 'H0')
w_avg,   w_std = value('w0wacdm_planck_boss_snu', 'w')
wa_avg,  wa_std = value('w0wacdm_planck_boss_snu', 'wa')
outputfile.write(r"$w_0w_a$CDM & Planck + CMASS + LOWZ + SN & %.4f (%i) & %.3f (%i) & %.1f (%i) & \nodata & %.2f (%i) & %.2f (%i) \\" % (omh2_avg,
                                                                                                                                         np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), w_avg, np.int(100*w_std), wa_avg, np.int(100*wa_std)) + '\n')
omh2_avg, omh2_std = value('w0wacdm_planck_bao_snu', 'omegamh2')
om_avg,  om_std = value('w0wacdm_planck_bao_snu', 'omegam')
h0_avg,  h0_std = value('w0wacdm_planck_bao_snu', 'H0')
w_avg,   w_std = value('w0wacdm_planck_bao_snu', 'w')
wa_avg,  wa_std = value('w0wacdm_planck_bao_snu', 'wa')
outputfile.write(r"$w_0w_a$CDM & Planck + BAO + SN & %.4f (%i) & %.3f (%i) & %.1f (%i) & \nodata & %.2f (%i) & %.2f (%i) \\" % (omh2_avg,
                                                                                                                                np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), w_avg, np.int(100*w_std), wa_avg, np.int(100*wa_std)) + '\n')
omh2_avg, omh2_std = value('w0wacdm_wmap_bao_snu', 'omegamh2')
om_avg,  om_std = value('w0wacdm_wmap_bao_snu', 'omegam')
h0_avg,  h0_std = value('w0wacdm_wmap_bao_snu', 'H0')
w_avg,   w_std = value('w0wacdm_wmap_bao_snu', 'w')
wa_avg,  wa_std = value('w0wacdm_wmap_bao_snu', 'wa')
outputfile.write(r"$w_0w_a$CDM & WMAP + BAO + SN & %.4f (%i) & %.3f (%i) & %.1f (%i) & \nodata & %.2f (%i) & %.2f (%i) \\" % (omh2_avg,
                                                                                                                              np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), w_avg, np.int(100*w_std), wa_avg, np.int(100*wa_std)) + '\n')
omh2_avg, omh2_std = value('w0wacdm_ewmap_bao_snu', 'omegamh2')
om_avg,  om_std = value('w0wacdm_ewmap_bao_snu', 'omegam')
h0_avg,  h0_std = value('w0wacdm_ewmap_bao_snu', 'H0')
w_avg,   w_std = value('w0wacdm_ewmap_bao_snu', 'w')
wa_avg,  wa_std = value('w0wacdm_ewmap_bao_snu', 'wa')
outputfile.write(r"$w_0w_a$CDM & \textit{e}WMAP + BAO + SN & %.4f (%i) & %.3f (%i) & %.1f (%i) & \nodata & %.2f (%i) & %.2f (%i) \\" % (
    omh2_avg, np.int(10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), w_avg, np.int(100*w_std), wa_avg, np.int(100*wa_std)) + '\n')
outputfile.write(r"\hline"+'\n')

# ow0waCDM
omh2_avg, omh2_std = value('ow0wacdm_planck_bossiso', 'omegamh2')
om_avg,  om_std = value('ow0wacdm_planck_bossiso', 'omegam')
h0_avg,  h0_std = value('ow0wacdm_planck_bossiso', 'H0')
ok_avg,  ok_std = value('ow0wacdm_planck_bossiso', 'omegak')
w_avg,   w_std = value('ow0wacdm_planck_bossiso', 'w')
wa_avg,  wa_std = value('ow0wacdm_planck_bossiso', 'wa')
outputfile.write(r"o$w_0w_a$CDM & Planck + CMASS-iso + LOWZ & %.4f (%i) & %.3f (%i) & %.1f (%i) & %+.4f (%i) & %.2f (%i) & %.2f (%i) \\" % (omh2_avg, np.int(10000 *
                                                                                                                                                             omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), ok_avg, np.int(10000*ok_std), w_avg, np.int(100*w_std), wa_avg, np.int(100*wa_std)) + '\n')
omh2_avg, omh2_std = value('ow0wacdm_planck_boss', 'omegamh2')
om_avg,  om_std = value('ow0wacdm_planck_boss', 'omegam')
h0_avg,  h0_std = value('ow0wacdm_planck_boss', 'H0')
ok_avg,  ok_std = value('ow0wacdm_planck_boss', 'omegak')
w_avg,   w_std = value('ow0wacdm_planck_boss', 'w')
wa_avg,  wa_std = value('ow0wacdm_planck_boss', 'wa')
outputfile.write(r"o$w_0w_a$CDM & Planck + CMASS + LOWZ & %.4f (%i) & %.3f (%i) & %.1f (%i) & %+.4f (%i) & %.2f (%i) & %.2f (%i) \\" % (omh2_avg, np.int(10000 *
                                                                                                                                                         omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), ok_avg, np.int(10000*ok_std), w_avg, np.int(100*w_std), wa_avg, np.int(100*wa_std)) + '\n')
omh2_avg, omh2_std = value('ow0wacdm_planck_bao', 'omegamh2')
om_avg,  om_std = value('ow0wacdm_planck_bao', 'omegam')
h0_avg,  h0_std = value('ow0wacdm_planck_bao', 'H0')
ok_avg,  ok_std = value('ow0wacdm_planck_bao', 'omegak')
w_avg,   w_std = value('ow0wacdm_planck_bao', 'w')
wa_avg,  wa_std = value('ow0wacdm_planck_bao', 'wa')
outputfile.write(r"o$w_0w_a$CDM & Planck + BAO & %.4f (%i) & %.3f (%i) & %.1f (%i) & %+.4f (%i) & %.2f (%i) & %.2f (%i) \\" % (omh2_avg, np.int(10000*omh2_std),
                                                                                                                               om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), ok_avg, np.int(10000*ok_std), w_avg, np.int(100*w_std), wa_avg, np.int(100*wa_std)) + '\n')
omh2_avg, omh2_std = value('ow0wacdm_planck_boss_snu', 'omegamh2')
om_avg,  om_std = value('ow0wacdm_planck_boss_snu', 'omegam')
h0_avg,  h0_std = value('ow0wacdm_planck_boss_snu', 'H0')
ok_avg,  ok_std = value('ow0wacdm_planck_boss_snu', 'omegak')
w_avg,   w_std = value('ow0wacdm_planck_boss_snu', 'w')
wa_avg,  wa_std = value('ow0wacdm_planck_boss_snu', 'wa')
outputfile.write(r"o$w_0w_a$CDM & Planck + CMASS + LOWZ + SN & %.4f (%i) & %.3f (%i) & %.1f (%i) & %+.4f (%i) & %.2f (%i) & %.2f (%i) \\" % (omh2_avg, np.int(10000 *
                                                                                                                                                              omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), ok_avg, np.int(10000*ok_std), w_avg, np.int(100*w_std), wa_avg, np.int(100*wa_std)) + '\n')
omh2_avg, omh2_std = value('ow0wacdm_planck_bao_snu', 'omegamh2')
om_avg,  om_std = value('ow0wacdm_planck_bao_snu', 'omegam')
h0_avg,  h0_std = value('ow0wacdm_planck_bao_snu', 'H0')
ok_avg,  ok_std = value('ow0wacdm_planck_bao_snu', 'omegak')
w_avg,   w_std = value('ow0wacdm_planck_bao_snu', 'w')
wa_avg,  wa_std = value('ow0wacdm_planck_bao_snu', 'wa')
outputfile.write(r"o$w_0w_a$CDM & Planck + BAO + SN & %.4f (%i) & %.3f (%i) & %.1f (%i) & %+.4f (%i) & %.2f (%i) & %.2f (%i) \\" % (omh2_avg, np.int(10000*omh2_std),
                                                                                                                                    om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), ok_avg, np.int(10000*ok_std), w_avg, np.int(100*w_std), wa_avg, np.int(100*wa_std)) + '\n')
omh2_avg, omh2_std = value('ow0wacdm_wmap_bao_snu', 'omegamh2')
om_avg,  om_std = value('ow0wacdm_wmap_bao_snu', 'omegam')
h0_avg,  h0_std = value('ow0wacdm_wmap_bao_snu', 'H0')
ok_avg,  ok_std = value('ow0wacdm_wmap_bao_snu', 'omegak')
w_avg,   w_std = value('ow0wacdm_wmap_bao_snu', 'w')
wa_avg,  wa_std = value('ow0wacdm_wmap_bao_snu', 'wa')
outputfile.write(r"o$w_0w_a$CDM & WMAP + BAO + SN & %.4f (%i) & %.3f (%i) & %.1f (%i) & %+.4f (%i) & %.2f (%i) & %.2f (%i) \\" % (omh2_avg, np.int(10000*omh2_std),
                                                                                                                                  om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), ok_avg, np.int(10000*ok_std), w_avg, np.int(100*w_std), wa_avg, np.int(100*wa_std)) + '\n')
omh2_avg, omh2_std = value('ow0wacdm_ewmap_bao_snu', 'omegamh2')
om_avg,  om_std = value('ow0wacdm_ewmap_bao_snu', 'omegam')
h0_avg,  h0_std = value('ow0wacdm_ewmap_bao_snu', 'H0')
ok_avg,  ok_std = value('ow0wacdm_ewmap_bao_snu', 'omegak')
w_avg,   w_std = value('ow0wacdm_ewmap_bao_snu', 'w')
wa_avg,  wa_std = value('ow0wacdm_ewmap_bao_snu', 'wa')
outputfile.write(r"o$w_0w_a$CDM & \textit{e}WMAP + BAO + SN & %.4f (%i) & %.3f (%i) & %.1f (%i) & %+.4f (%i) & %.2f (%i) & %.2f (%i) \\" % (omh2_avg, np.int(
    10000*omh2_std), om_avg, np.int(1000*om_std), h0_avg, np.int(10*h0_std), ok_avg, np.int(10000*ok_std), w_avg, np.int(100*w_std), wa_avg, np.int(100*wa_std)) + '\n')
outputfile.write(r"\hline"+'\n')

outputfile.write(r"\end{tabular}"+'\n')
outputfile.write(r"\caption{Cosmological constraints by different datasets in the cosmological models $\Lambda$CDM, oCDM, $w$CDM, o$w$CDM, $w_0w_a$CDM, and o$w_0w_a$CDM. We compare the cosmological constraints from combining Planck with acoustic scale from BOSS galaxies as well as lower and higher redshift BAO measurements from the 6-degree field galaxy redshift survey (6DF) and the BOSS-Lyman alpha forest (Ly$\alpha$F), respectively. We also compare how these combinations benefit from the constraining power of type-Ia Supernovae from the Union 2 compilation by the Supernovae Cosmology Project (SN). The WMAP and \textit{e}WMAP cases have been added for comparison. As in Table~\ref{tab:bigcos}, 'CMASS-iso' means the isotropic measurement from the CMASS sample, whereas the anisotropic one is referred to simply as 'CMASS'. 'LOWZ' is the isotropic measurement from the LOWZ sample. 'BAO' stands for the combination CMASS + LOWZ + 6DF + Ly$\alpha$F.}"+'\n')
# $H_0$ is in units of km s$^{-1}$ Mpc$^{-1}$.
outputfile.write(r"\label{tab:bigcos2}"+'\n')
outputfile.write(r"\end{table*}"+'\n')

outputfile.close()
