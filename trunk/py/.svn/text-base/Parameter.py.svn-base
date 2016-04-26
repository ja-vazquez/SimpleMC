#
# Simple class for dealing with Parameter.
# Parameter has a name, a value, an error and some bounds
# Names are also latex names.

class Parameter:
    def __init__(self,name,value,err=0, bounds=None,Ltxname=None):
        self.name=name
        self.Ltxname = Ltxname
        self.value=value
        # this is the estimate of error
        self.error=err
        if bounds==None:
            self.bounds=(value-5*err, value+5*err)
        else:
            self.bounds=bounds
    
    def sameParam(self,param2):
        return self.name==param2.name
  
    def setLatexName(self,Ltx):
        self.Ltxname=Ltx

    def setValue(self,val):
        self.value=val

    def setError(self,err):
        self.error=err

    def setBounds (self, low,high):
        self.bounds=[low,high]


