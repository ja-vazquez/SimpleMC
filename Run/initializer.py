from RunBase import ParseModel, ParseDataset
from numpy import array

def iniReader(inifile):
    import sys
    settings = []
    if sys.version_info > (3, 0):
        import configparser
        config = configparser.ConfigParser()
        config.read(inifile)
        chainsdir = config['custom']['chainsdir']
        model = config['custom']['model']
        prefact = config['custom']['prefact']
        datasets = config['custom']['datasets']
        samplername = config['custom']['sampler']
        plottername = config['custom']['plotter']
        settings.extend([chainsdir,model,prefact,datasets,samplername,plottername])
        
        if samplername == 'mh':
            nsamp = int(config['metropolis']['nsamp'])
            skip = int(config['metropolis']['skip'])
            settings.extend([nsamp,skip])
        
        if samplername in ['snest','mnest','sbambi','bambi']:
            print('\nUsing NESTLE')
            nlivepoints = int(config['nested']['nlivepoints'])
            accuracy = float(config['nested']['accuracy'])
            priortype = config['nested']['priortype'] 
            settings.extend([nlivepoints,accuracy,priortype])
            
            print("The number of live points is: %d"%(nlivepoints))
            if samplername == 'mnest':
                print('\nThe sampler used is MULTINEST. Feroz et al (2009)\n')
            elif samplername == 'snest': 
                print('\nMukherjee, Parkinson & Liddle (2006) nested sampling\n')
            if samplername in ['sbambi','bambi']: 
                print("\nANN based on pybambi. Graff et al (2012).\n")
            if priortype == 'g': 
                print("Using %d-sigmas for the gaussian prior.\n"%(nsigma))
            else: 
                print("Using flat priors...")
           
    else:
        #python2
        print("We highly recommend you to use Python3, however SuperMC can run with Python2.")
        import ConfigParser
        f = open(inifile,'r')
        config = ConfigParser.ConfigParser()
        config.readfp(f)
        chainsdir = config.get('custom','chainsdir')
        model = config.get('custom','model')
        prefact = config.get('custom','prefact')
        datasets = config.get('custom','datasets')
        samplername = config.get('custom','sampler')
        plottername = config.get('custom','plotter')

        settings.extend([chainsdir,model,prefact,datasets,samplername,plottername])
        
        if samplername == 'mh':
            nsamp = int(config.get('metropolis','nsamp'))
            skip = int(config.get('metropolis','skip'))
            settings.extend([nsamp,skip])
            print("The sampler used is Metropolis-Hastings")
            
        elif samplername in ['snest','mnest','sbambi','bambi']:
            print('\nUsing NESTLE')
            accuracy = float(config.get('nested','accuracy'))
            nlivepoints = int(config.get('nested','nlivepoints'))
            priortype = config.get('nested','priortype')
            settings.extend([nlivepoints,accuracy,priortype])
            if samplername == 'mnest':
                print('\nThe sampler used is MULTINEST. Feroz et al (2009)\n')
            elif samplername == 'snest': 
                print('\nMukherjee, Parkinson & Liddle (2006) nested sampling\n')
            if samplername in ['sbambi','bambi']: 
                print("\nANN based on pybambi. Graff et al (2012).\n")
        f.close()

    print("Running:\n", str(settings).lstrip('[').rstrip(']').strip(','))

    return settings

    
def TLinit(model,datasets):
    """
    Returns T evaluated at the model, and L in at the datasets.
    
    Parameters:

    model:      str chosen model
    datasets:   str chosen datasets
    """
    T = ParseModel(model)
    L = ParseDataset(datasets)        
    L.setTheory(T)
    return T, L

def priorValues(T,L):
    """Returns some values for the priorTransform """
    bounds = SetPriors(T)[3]
    means = array(SetPriors(T)[1])
    nsigma = 4. 
    return bounds, means, nsigma


def priorsTransform(theta, bounds, priortype):
    """prior Transform for gaussian and flat priors"""
    priors = []
    if priortype == 'g':
        for c, bound in enumerate(bounds):
            mu = means[c]
            sigma = (bound[1]-bound[0])/n
            priors.append(mu+sigma*(ndtri(theta[c])))
    elif priortype == 'g2':
        for c, bound in enumerate(bounds):
            mu = means[c]
            sigma = (bound[1]-bound[0])/n
            priors.append( -0.5*((theta[c]-mu)/sigma)**2)
    else:
        for c, bound in enumerate(bounds):
        #When theta 0-> append bound[0], if theta 1-> append bound[1]
            priors.append(theta[c]*(bound[1]-bound[0])+bound[0])
            #At this moment, np.array(priors) has shape (dims,) 
    return array(priors)


def getDims(T):
    """Returns the numbers of dimensions and a parameters list"""
    #IGV: We need the names of parameters in a list and, on the other hand,
    #    the dimensions. I don't found a fancy way, probably it exists.  
    freeP = T.freeParameters()
    listP = []
    for dims,item in enumerate(freeP):
        listP.append(item.name)
    dims+=1
    return dims, listP

def SetPriors(T):
    """After setTheory, you can set the priors"""
    parameters = T.freeParameters()
    names=[]
    values=[]
    errorlist=[]
    boundlist=[]
    latexnames=[]
    for parameter in parameters:
        names.append(parameter.name)
        values.append(parameter.value)
        errorlist.append(parameter.error)
        boundlist.append(parameter.bounds)
        latexnames.append(parameter.Ltxname)
    return [names, values, errorlist, boundlist, latexnames]

#IGV:
def instantiatePars(T,value):
    """This method returns an instance of 
    Parameter objects with the sampler values"""
    aux=[]
    for item in value:
        aux.append(item)
    names, values, errorlist, boundlist, latexnames = SetPriors(T)
    instances = []
    for i in range(len(names)):
        instances.append(Parameter(names[i], 
            aux[i], errorlist[i], boundlist[i], latexnames[i]))
    return instances


def writteSummary(chainsdir, outputname, time,*args):
    file = open(chainsdir + '/' + outputname + "_Summary" + ".txt",'w')
    file.write('Summary:\n')
    for item in args:
        if type(item) is list:
            for element in item:
                file.write(str(element)+'\n')
        else:
            file.write(str(item)+'\n')
    print("\nElapsed time: %.3f minutes = %.3f seconds"%(time/60,time))  
    file.write('\nElapsed time: %.3f minutes = %.3f seconds'%(time/60,time))  
    file.close()


#next logprior def is for emcee
#def logposterior(theta, loglikelihood):
    """
    The natural logarithm of the joint posterior.
    
    Args:
        theta (tuple): a sample containing individual parameter values
        data (list): the set of data/observations
        sigma (float): the standard deviation of the data points
        x (list): the abscissa values at which the data/model is defined
    """
#    lp = priorsTransform(theta) # get the prior
    
    # if the prior is not finite return a probability of zero (log probability of -inf)
#    if not np.isfinite(lp):
#        return -np.inf
    
    # return the likeihood times the prior (log likelihood plus the log prior)
#    return lp + loglikelihood(theta)

