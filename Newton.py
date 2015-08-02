'''
    Solving non-linear system of equations by using Newtons method
'''
from __future__ import division
import ROOT as rt
import numpy as np
import time as tm
from math import *

class Newton:
    def __init__( self, T0, RMSmax = 1e-3, ITmax = 1e4 ):
        self.rmsmax = RMSmax
        self.itmax  = int(ITmax)
        self._Compute(np.array(T0))

    def Residuals( self, T ):
        raise NotImplementedError('Residuals not implemented.')
    
    def Jacobian( self, T ):
        raise NotImplementedError('Jacobian not implemented.')

    def TrueT( self ):
        raise NotImplementedError('True value not implemented.')

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
        while rmsT > self.rmsmax and it < self.itmax :
            R    = self._R(T0)
            J    = self._J(T0)
            T    = T0 - np.linalg.inv(J).dot( R )
            rmsR = sqrt( sum( R*R ) / n )
            rmsT = sqrt( sum( (T-Ttrue)*(T-Ttrue) ) / n )
            print rmsR, rmsT
            T0   = T
            it  += 1
        
        # save data
        self.data = T
        self.rmsR = rmsR
        self.rmsT = rmsT
        self.iter = it

    def PrintResult( self ):
        print 'Solution of {0} after {1} iterations with an RMS residual of {2}:\nT = {3}\n'.format( self.__class__.__name__, self.iter, self.rmsT, self.data )

class FlatPlateCollector(Newton):
    cs = [0.06823, 0.05848, 0.01509, 2.0, 0.11696, 2.05, 0.2534, 0.06698]
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


class Burguers:
    def Residuals( self, T ):
        pass


trial = FlatPlateCollector( [.3]*3, 1e-4, 10 )
trial.PrintResult()




















'''
### SIMULATION VARIABLES

dx = 1e-1                  # spatial step
dt = 500                   # time step

t0 = 0.                    # initial time
t1 = 3e3                   # final time

x0 = 0.                    # initial point
x1 = 1.1                   # final point

nx = int( (x1-x0) / dx )   # number of points in space
nt = int( (t1-t0) / dt )   # number of points in time

ntrue = 10                 # number of terms used in true expression

alpha = 1e-5               # diffusion coefficient
s1    = alpha * dt / dx**2 # constant 1
s2    = 1 - 2*s1           # constant 2

data = np.zeros( (nt,nx) ) # simulated data matrix
true = np.zeros( (nt,nx) ) # true data matrix

assert s1 <= 0.5,'Not stable. S1 must be less or equal than 0.5 and it is '+ str(s1)

###### BOUNDARY CONDITIONS
T0 = 100.                       # input temperature

data[0,:]  = true[0,:]  = 0.    # at time 0, T = 0

data[1:,0]  = true[1:,0]  = 1.0 # T = T0 for the boundary at all times
data[1:,-1] = true[1:,-1] = 1.0 # i.e. perfect temperature source.

data[:1,0]  = true[:1,0]  = 0.5 # start with half temperature
data[:1,-1] = true[:1,-1] = 0.5 # start with half temperature


### SIMULATION
for t in range(nt-1):
    tr = (t+1) * dt 
    for x in range(1,nx-1):
        xn = x / ( nx - 1 ) 
        data[t+1,x] = s2 * data[t,x] + s1 * ( data[t,x-1] + data[t,x+1] ) # FTCS scheme
        true[t+1,x] = 1 - sum( 4 / pi / ( 2*m+1 ) * sin( (2*m+1)*pi*xn ) * exp(-alpha*pi**2*tr*(2*m+1)**2) for m in range(ntrue) )

data *= T0
true *= T0

### PLOTTING
cv = rt.TCanvas()
hdata = [ rt.TH1F( 'data_' + str(i),';x;T', nx, x0, x1 ) for i in range(nt) ]
htrue = [ rt.TH1F( 'true_' + str(i),';x;T', nx, x0, x1 ) for i in range(nt) ]

# visual format
for h in hdata: h.SetMaximum(1.1*T0); h.SetMinimum(0.)
for h in htrue: h.SetMinimum(0.); h.SetMaximum(1.1*T0); h.SetLineColor(rt.kRed)

for t in range(nt):
    for x in range(nx):
        hdata[t].SetBinContent(x+1,data[t,x])
        htrue[t].SetBinContent(x+1,true[t,x])
    hdata[t].Draw()
    htrue[t].Draw('same')
    cv.Update()
#    raw_input()
    tm.sleep(0.05)


    

'''
