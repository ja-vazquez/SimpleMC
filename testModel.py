from simplemc.cosmo.Parameter import Parameter
from simplemc.DriverMC import DriverMC

### This scripts create a generic model with SimpleModel and SimpleCosmoModel classes at backend

# 1) Define your parameters objects
# name string,  value intermediate, step size,
# (bound_inf, bound_sup),  LaTeX name
# The name of the variables must to be the same at the name of the Parameter
m = Parameter("m", 0, 0.05, (0, 0.1), "m_0")
# m = Parameter("m", 3, 0.5, (0,5), "m_0")
b = Parameter("b", 3, 0.05, (0, 5), "b_0")

# 2) Create a list with your parameters objects
parameterlist = [m, b]

# 3) Define a method that reads a list of parameters,
# unzip them and return the a function of x with the
# parameters.
def model(parameterlist, x):
    m, b = parameterlist
    return m*x+b

cosmo_model = 'Ocb/a**3+Omrad/a**4+NuContrib+(1.0-Om-m)'

# 4) Use SimpleMC as usually, but with model = custom_model

analyzer = DriverMC(model='simple', datasets='dline', analyzername='mcmc',
                    custom_parameters=parameterlist, custom_function=model)

# analyzer = DriverMC(model='simple_cosmo', datasets='SN', analyzername='mcmc',
#                     custom_parameters=parameterlist, custom_function=cosmo_model)


analyzer.executer(nsamp=1000)