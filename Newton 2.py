'''
    Abstract implementation of the Newton method for solving non-linear
    systems of equations.
'''
from __future__ import division
import numpy as np
from math import *

class Newton:
    '''
        Abstract class for solving a system of non-linear equations.
        The actual class must inheritate from this one and implement the
        Residuals and Jacobian methods in order to produce the result. An
        extra method (TrueT) may be implemented in order to compare the
        solution with the real values (through the RMS).
        The computation is performed by initializating the class. No extra
        steps are needed. The results are stores as attributes:
        - data: the vector of results
        - rmsR: RMS of the residuals
        - rmsT: RMS of the solution wrt the real value (if provided)
        - iter: number of iterations performed
        Also, these data may be printed with the PrintResult method.
    '''
    def __init__( self, T0, RMSmax = 1e-3, ITmax = 1e4 ):
        self.rmsmax = RMSmax
        self.itmax  = int(ITmax)
        self._Compute(np.array(T0))

    def Residuals( self, T ):
        raise NotImplementedError('Residuals not implemented.')

    def Jacobian( self, T ):
        raise NotImplementedError('Jacobian not implemented.')

    def TrueT( self ):
        return None
#        raise NotImplementedError('True value not implemented.')

    def _R( self, T = None ):
        return np.array( self.Residuals(T) )

    def _J( self, T = None ):
        return np.array( self.Jacobian(T) )

    def _T( self ):
        return np.array( self.TrueT() )

    def _Compute( self, T0 ):
        n     = len(T0)
        rmsT  = rmsR = float('inf')
        it    = 0
        Ttrue = self._T()
        compare = not (Ttrue is None)
        while True :
            R    = self._R(T0)
            J    = self._J(T0)
            T    = T0 - np.linalg.inv(J).dot( R )
            rmsR = sqrt( sum( R*R ) / n )
            rmsT = sqrt( sum( (T-Ttrue)*(T-Ttrue) ) / n ) if compare else None
            T0   = T
            it  += 1
            if it == self.itmax or (compare and rmsT < self.rmsmax):
                break

        # save data
        self.data = T
        self.rmsR = rmsR
        self.rmsT = rmsT
        self.iter = it

    def PrintResult( self ):
        print 'Solution of {0} after {1} iterations with an RMS residual of {2}:\nT = {3}\n'.format( self.__class__.__name__, self.iter, self.rmsT, self.data )
