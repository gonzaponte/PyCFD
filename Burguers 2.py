'''
    Test Newton method by solving 2D- Burguers equation.
'''
from __future__ import division
from Newton import *
from math import *

class Burguers(Newton):
    Nx  = 5
    Ny  = 5
    Nxm = Nx - 2
    Nym = Ny - 2
    N   = Nxm * Nym * 2
    Nm  = N - 2
    dx  = 0.1
    dy  = 0.1
    Re  = 10.
    ReI = 1/Re
    s1x = 0.5/dx
    s2x = ReI/dx**2
    s1y = 0.5/dy
    s2y = ReI/dy**2

    xiter = range(1,Nx-1)
    yiter = range(1,Ny-1)

    def matrix_shaped( self, T ):
        return np.array( [T[ i0 : i0 + self.Nym ] for i0 in range(0,self.N//2,self.Nym)] )

    def Residuals( self, T ):
        Tu = self.matrix_shaped(T[::2])
        Tv = self.matrix_shaped(T[1::2])

        R = []
        for j in self.xiter:
            jp = j+1
            jm = j-1
            for k in self.yiter:
                Ru = Tu[j,k] * self.s1x * (Tu[jp,k]-Tu[jm,k]) + Tv[j,k] * self.s1y * (Tu[jp,k]-Tu[jm,k]) - self.s2x * (Tu[jm,k] - 2 * Tu[j,k] + Tu[jp,k]) - self.s2y * (Tu[jm,k] - 2 * Tu[j,k] + Tu[jp,k])
                Rv = Tu[j,k] * self.s1x * (Tv[jp,k]-Tv[jm,k]) + Tv[j,k] * self.s1y * (Tv[jp,k]-Tv[jm,k]) - self.s2x * (Tv[jm,k] - 2 * Tv[j,k] + Tv[jp,k]) - self.s2y * (Tv[jm,k] - 2 * Tv[j,k] + Tv[jp,k])
                R.extend( [Ru,Rv] )
        return R

    def Jacobian( self, T ):
        Tu = self.matrix_shaped(T[::2])
        Tv = self.matrix_shaped(T[1::2])

        J = [ [0. for i in range(self.N)] for j in range(self.N) ]
        for j in self.xiter:
            jp = j+1
            jm = j-1
            for k in self.yiter:
                kp = k+1
                km = k-1
                uindex  = 2 * (jm * self.Nym + k)
                vindex  = uindex + 1
                umindex = uindex - 2 * self.Nym
                vmindex = umindex + 1
                upindex = uindex + 2 * self.Nym
                vpindex = upindex + 1

                J[uindex][uindex] = 2 * ( self.s2x + self.s2y ) + self.s1x * ( Tu[jp,k] - Tu[jm,k] )
                J[uindex][vindex] = self.s1y * ( Tu[j,kp] - Tu[j,km] )
                J[vindex][vindex] = 2 * ( self.s2x + self.s2y ) + self.s1y * ( Tv[j,kp] - Tv[j,km] )
                J[vindex][uindex] = self.s1x * ( Tv[jp,k] - Tv[jm,k] )

                if umindex >= 0:
                    J[uindex][umindex] = J[vindex][vmindex] = - self.s1x * Tu[j,k] - self.s2x
                if vpindex < N:
                    J[uindex][upindex] = J[vindex][vpindex] = + self.s1x * Tu[j,k] - self.s2x
                if uindex > 1:
                    J[uindex][uindex-2] = J[vindex][vindex-2] = - self.s1y *Tv[j,k] - self.s2y
                if vindex < self.Nm:
                    J[uindex][uindex+2] = J[vindex][vindex+2] = + self.s1y *Tv[j,k] - self.s2y

        return J

    def TrueT( self ):
        a0 = a1 = 110.13
        a2, a3 = 0.
        a4 = 1.0 * 2
        lb = 5.
        x0 = 1.0
        Re = - 2.0 * self.ReI
        ymax = pi/6/
        phi   = lambda x,y: a0 + a1*x + a2 * y + a3 * x * y + a4 * cosh(lb*(x-x0)) * cos(lb * y)
        phi_x = lambda x,y: a1 + a3 * y + a4 * lb * sinh(lb*(x-x0)) * cos(lb*y)
        phi_y = lambda x,y: a2 + a3 * x -    a4 * lb * cosh(lb*(x-x0)) * sin(lb*y)

        dx = 2./(self.Nx-1)
        dy = 2./(self.Ny-1)
        true = []
        for j in range(self.Nx):
            x = -1. + j*dx
            for k in self.yiter:
                u = Re * phi_x(x,y)/phi(x,y)
                v = Re * phi_y(x,y)/phi(x,y)
                true.extend( [u,v] )
        return true


        return None
