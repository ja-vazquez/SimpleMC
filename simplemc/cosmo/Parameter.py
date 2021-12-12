
class Parameter:
    """
    A simple class for dealing with Parameter.
    Parameter has a name, a value, an error and some bounds
    Names have also latex names.

    Parameters
    ----------
    name : str
        Name to identify the parameter.

    value : float
        Default value. In MCMC will be the starting value.

    err: float
        Estimated error, in MCMC will be the step of the chain.
        Default is 0, but better to write any positive number.

    bounds: list
        Priors. (minimum value, maximum value).
        List of pairs of bounds.
        Default is None, but computed as '(value-5*err, value+5*err)'.

    Ltxname: str
        Latex name, use mainly for plotting.
        Default is None, and in this case uses the 'name' string.


    Example
    -------
    The hubble parameter
    h_par = Parameter('h', 0.6821,  0.05,   (0.4, 1.0),    'h')
    """

    def __init__(self, name, value, err=0.0, bounds=None, Ltxname=None):

        # Initialize name and Latex name.
        self.name = name
        if Ltxname:
            self.Ltxname = Ltxname
        else:
            self.Ltxname = name
        self.value = value

        # Initialize the estimated of error.
        self.error = err

        # Initialize the priors.
        if bounds == None:
            self.bounds = (value - 5*err, value + 5*err)
        else:
            self.bounds = bounds


    def sameParam(self, param2):
        return self.name == param2.name


    def setLatexName(self, Ltx):
        self.Ltxname = Ltx


    def setValue(self, val):
        self.value = val


    def setError(self, err):
        self.error = err


    def setBounds(self, low, high):
        self.bounds = [low, high]
