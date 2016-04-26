# To produce some plots of the convergence values

set term post eps enhanced color "Times Roman" 25
set output "Converg.eps"
set size 1.1, 1.1

set style line 1 lt 1 lw 3 pt 3 lc rgb "#0000FF"
set style line 2 lt 1 lw 3 pt 3 lc rgb "#DC143C"
set style line 3 lt 1 lw 3 pt 3 lc rgb "#008000"
set style line 4 lt 1 lw 3 pt 3 lc rgb "#CCFF00"
set style line 6 lt 2 lw 1 pt 3 lc rgb "#000000"

set logscale y
set format y "10^{%L}"
set format x "%2.0t{/Symbol \264}10^{%L}"
#set format x "10^{%L}"

set xlabel "{# Like evaluations}" font "Times-Roman,30"
set ylabel "{Conver}" font "Times-Italic,25" 
set xtics 0,20000,10e5

plot 'MCMC_owCDM.converge' w l title "MCMC" ls 1 ,  'chains/owCDM_phy_BBAO+Planck_GS2000.converge' u 2:3  w l title "GSampler" ls 2, 'chains/owCDM_phy_BBAO+Planck_GS2000_2.converge' u 2:3 w l title "{See caption}" ls 3
