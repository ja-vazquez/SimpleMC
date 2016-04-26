#!/usr/bin/env python

# first, let's commit
import os, time, subprocess


def mysystem(com):
    print "Running ", com
    os.system(com)
#get svn version
svnversion=subprocess.Popen('svnversion', stdout=subprocess.PIPE).stdout.read()[:-1]
if ':' in svnversion:
    svnversion=svnversion[:svnversion.find(':')]+"M"

date=time.strftime("%y%m%d")
nm="chains_%s_%s"%(date,svnversion)

open('Pack/latest.py','a').write('latest_chain="%s"\n'%nm)
mysystem("cp -rv chains "+nm)

#mysystem("tar zcvf /astro/astronfs01/www/anze/achains/%s.tgz %s/"%(nm,nm))
#mysystem ("svn commit Pack/latest.py -m 'new chains available'")

print "tar zcvf /astro/astronfs01/www/anze/achains/%s.tgz %s/"%(nm,nm)
print "svn commit Pack/latest.py -m 'new chains available'"




