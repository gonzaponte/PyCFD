'''
    Test Newton method by solving 2D- Burguers equation.
'''
from __future__ import division
from Newton import *

class Burguers(Newton):
    Nx  = 5
    Ny  = 5
    N   = (Nx-2) * (Ny-2) * 2
    dx  = 0.1
    dy  = 0.1
    s1x = 0.5/dx
    s2x = 1/dx**2
    s1y = 0.5/dy
    s2y = 1/dy**2
    Re  = 10.
    ReI = 1/Re
    def Residuals( self, T ):
        Tu, Tv = T
        Ru = [ Tu[j,k] * s1x * (Tu[j+1,k]-Tu[j-1,k]) + Tv[j,k] * s1y * (Tu[j+1,k]-Tu[j-1,k]) - ReI * ( s2x * (Tu[j-1,k] - 2 * Tu[j,k] + Tu[j+1,k]) + s2y * (Tu[j-1,k] - 2 * Tu[j,k] + Tu[j+1,k]) ) for j in range(1,Nx-1) for k in range(1,Ny-1) ]
        Rv = [ Tu[j,k] * s1x * (Tv[j+1,k]-Tv[j-1,k]) + Tv[j,k] * s1y * (Tv[j+1,k]-Tv[j-1,k]) - ReI * ( s2x * (Tv[j-1,k] - 2 * Tv[j,k] + Tv[j+1,k]) + s2y * (Tv[j-1,k] - 2 * Tv[j,k] + Tv[j+1,k]) ) for j in range(1,Nx-1) for k in range(1,Ny-1) ]

    def Jacobian( self, T ):
        return None
    
    def TrueT( self ):
        return None


