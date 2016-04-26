#!/usr/bin/env python

from contour import *

def compute_hz(T,z) :
    # H(z)/(1+z)
    return 100*T.h*sqrt(T.RHSquared_a(1/(1+z)))/(1+z)

if len(sys.argv)<4 :
    print "usage:"
    print sys.argv[0]," chains_dir cosmo data"
    print "example:"
    print sys.argv[0]," chains_140411_155701 owaCDM BBAO+CMBP"
    sys.exit(12)

chaindir=sys.argv[1]
cosmo=sys.argv[2]
data=sys.argv[3]

z,hz_1sigma_minus,hz_best,hz_1sigma_plus = contour(chaindir,cosmo,data,compute_hz)

if False :
    pylab.plot(z,hz_best)
    pylab.plot(z,hz_1sigma_minus)
    pylab.plot(z,hz_1sigma_plus)
    pylab.show()


# write contours
filename="AlternativePlots/hz-%s-%s-%s.list"%(cosmo,data,chaindir)
file=open(filename,"w")
file.write("# hz=H(z)/(1+z)\n")
file.write("# z  hz(-1sigma) hz(best) hz(+1sigma)\n")
for i in range(z.shape[0]) :
    file.write("%f %f %f %f\n"%(z[i],hz_1sigma_minus[i],hz_best[i],hz_1sigma_plus[i]))
file.close()
print "wrote",filename

sys.exit(0)

