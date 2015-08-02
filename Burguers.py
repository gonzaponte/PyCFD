'''
    Test Newton method by solving Burguers equation.
'''
from Newton import *

class Burguers(Newton):

    def Residuals( self, T ):
        return None

    def Jacobian( self, T ):
        return None
    
    def TrueT( self ):
        return None


