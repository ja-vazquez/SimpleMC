from simplemc.models.SimpleModel import SimpleModel
from simplemc.cosmo.Parameter import Parameter
from simplemc.DriverMC import DriverMC

### This scripts create a generic model (without cosmology)
### with the SimpleModel class.

# 1) Define your parameters objects
# name string,  value intermediate, step size,
# (bound_inf, bound_sup),  LaTeX name
m = Parameter("m", 3, 0.5, (0,5), "m_0")
b = Parameter("b", 3, 0.5, (0,5), "b_0")

# 2) Create a list with your parameters objects
parameterlist = [m, b]

# 3) Define a method that reads a list of parameters,
# unzip them and return the a function of x with the
# parameters.
def model(parameterlist, x):
    m, b = parameterlist
    return m*x+b

# 4) Use SimpleMC as usually, but with model = custom_model
analyzer = DriverMC(model='custom', datasets='dline', analyzername='mcmc',
                    custom_parameters=parameterlist, custom_function=model)
analyzer.executer()
fig = analyzer.plot(show=True)
fig.simpleCorner()


