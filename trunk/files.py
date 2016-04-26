#!/usr/bin/env python
import pylab

models =['Decay','EarlyDE','LCDM','nuLCDM','nuoLCDM','nuwCDM','oLCDM','owaCDM','owCDM','PolyCDM','SlowRDE','StepCDM','TLight','waCDM','wCDM','WeirdCDM']

for model in models:
	print 'tar -cvf %s.tar %s*'%(model, model)
  	print 'gzip %s.tar'%(model)
	print 'mv %s.tar.gz /astro/astronfs01/www/jvazquez/chains_SimpleMC'%(model) 
