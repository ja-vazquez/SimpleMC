from __future__ import print_function
import corner
import numpy as np

class Plotter:
	"""This class make plots with getdist or corner."""

	def __init__(self, outputname, chainsdir, samplername, ndim, truths, weights, plottername):

		self.outputname = outputname
		self.chainsdir = chainsdir
		self.samplername = samplername
		self.ndim = ndim
		self.truths = truths
		self.weights = weights
		self.plottername = plottername

		self.parsfile = chainsdir + '/' + outputname.rstrip(".txt")+".paramnames"

		self.readFile()

		if plottername == "getdist":
			self.getdistPlotter()
		elif plottername == "corner":
			self.cornerPlotter()
		elif plottername == "cosmich":
			self.cosmichPlotter()

		self.showFig()


	def cornerPlotter(self):
		"""
		Corner method. 
		"""
		print("Plotting with Corner!")
		figure = corner.corner(self.samples, labels=self.latexnames[0:self.ndim],\
	                       bins = 20,\
	                       weights = self.weights,\
	                       color='g',\
	                       quantiles=[0.5],\
	                       show_titles=True,\
	                       title_fmt = '.4f',\
	                       smooth1d=True,\
	                       smooth=True,\
	                       fill_contours=True,\
	                       plot_contours =True,\
	                       plot_density=True,\
	                       #truths=truths, truth_color='#4682b4',\
	                       title_kwargs={"fontsize": 12})
		figure.savefig(self.chainsdir+'/'+self.outputname+".png")
		figure.savefig(self.chainsdir+'/'+self.outputname+".pdf")

		

	def getdistPlotter(self):
		"""
	    Getdist method.
		"""
		from getdist import plots, MCSamples, mcsamples
		import matplotlib.pyplot as plt
		print("Plotting with Getdist!")

		root = self.rootDirWithFile()
		#print(self.chainsdir, self.outputname)
		mcsamplefile = mcsamples.loadMCSamples(root, ini=None, jobItem=None, no_cache=False)
		g = plots.getSubplotPlotter()
		g.settings.figure_legend_frame = False
		g.settings.alpha_filled_add=0.4
		g.settings.title_limit_fontsize = 13
		g.triangle_plot(mcsamplefile, self.paramnames[0:self.ndim+2], filled_compare=True, 
			contour_colors=['green'],title_limit=1)

		plt.savefig(self.chainsdir+'/'+self.outputname+".png")
		plt.savefig(self.chainsdir+'/'+self.outputname+".pdf")

	def cosmichPlotter(self):
		"""
        cosmich plotter method
		"""
		import cosmich
		root = self.rootDirWithFile()+".txt"
		cosmich.cosmochain(root, nums=None)

	def readFile(self):
		"""
 		This method reads the samples and the .param file. 
		"""
		labelsfile = open(self.parsfile,'r')
		self.latexnames = []
		self.paramnames = []

		for item in labelsfile: 
		    self.latexnames.append('$'+item.split('			')[1].strip('\n')+'$')
		    self.paramnames.append(item.split('			')[0].strip('\n'))

		labelsfile.close()

		if self.samplername == 'mh':
			npchain = np.loadtxt(self.chainsdir + '/'+self.outputname+'_1.txt')
		else:
			npchain = np.loadtxt(self.chainsdir + '/'+self.outputname)
			self.outputname = self.outputname.rstrip(".txt")

		self.samples = npchain[:,2:self.ndim+2]

		

	def showFig(self):
		"""
		Displays the plots after the samples are generated.
		"""
		from PIL import Image
		img = Image.open(self.chainsdir+'/'+self.outputname+".png")
		img.show()

	def rootDirWithFile(self):
		"""
		Finds the path.
		"""
		import os
		dirpath = os.getcwd()
		root = dirpath+'/'+self.chainsdir +'/'+ self.outputname
		return root