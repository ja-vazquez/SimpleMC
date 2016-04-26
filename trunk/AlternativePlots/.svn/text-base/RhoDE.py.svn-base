#!/usr/bin/env python

from contour import *

def compute_rho(T,z) :
    # dark energy density / rhoc
    return T.Rho_de(1/(1+z))



if len(sys.argv)<4 :
    print "usage:"
    print sys.argv[0]," chains_dir cosmo data"
    print "example:"
    print sys.argv[0]," chains_140411_155701 owaCDM BBAO+CMBP"
    sys.exit(12)


chaindir=sys.argv[1]
cosmo=sys.argv[2]
data=sys.argv[3]

z,rho_1sigma_minus,rho_best,rho_1sigma_plus = contour(chaindir,cosmo,data,compute_rho)

if False :
    pylab.plot(z,rho_best)
    pylab.plot(z,rho_1sigma_minus)
    pylab.plot(z,rho_1sigma_plus)
    pylab.show()


# write contours
filename="AlternativePlots/rhode-%s-%s-%s.list"%(cosmo,data,chaindir)
file=open(filename,"w")
file.write("# rho=rho_DE(z)/rho_crit\n")
file.write("# z  rho(-1sigma) rho(best) rho(+1sigma)\n")
for i in range(z.shape[0]) :
    file.write("%f %f %f %f\n"%(z[i],rho_1sigma_minus[i],rho_best[i],rho_1sigma_plus[i]))
file.close()
print "wrote",filename

sys.exit(0)

