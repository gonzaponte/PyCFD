'''
    Test Newton method with a solar cell (or something like that).
'''
from Newton import *

class FlatPlateCollector(Newton):
    cs = (0.06823, 0.05848, 0.01509, 2.0, 0.11696, 2.05, 0.2534, 0.06698)
    def Residuals( self, T ):
        T0, T1, T2 = T
        c0, c1, c2, c3, c4, c5, c6, c7 = self.cs
        
        R1 = T0**4 + c0 * T0 -      T1**4 - c1 * T1 - c2
        R2 = T0**4 + c1 * T0 - c3 * T1**4 - c4 * T1 + T2**4 + c1 * T2
        R3 = T1**4 + c1 * T1 - c5 * T2**4 - c6 * T2 + c7  
        return R1, R2, R3

    def Jacobian( self, T ):
        T0, T1, T2 = T
        c0, c1, c2, c3, c4, c5, c6, c7 = self.cs
        
        J0 = [ 4 * T0**3 + c0, - 4 *      T1**3 - c1,   0.                  ]
        J1 = [ 4 * T0**3 + c1, - 4 * c3 * T1**3 - c4,   4 *      T2**3 + c1 ]
        J2 = [ 0.            ,   4 *      T1**3 + c1, - 4 * c5 * T2**3 - c6 ]
        return J0, J1, J2
    
    def TrueT( self ):
        return [.415,.379,.334]

trial = FlatPlateCollector( [.3]*3, 1e-4, 10 )
trial.PrintResult()

