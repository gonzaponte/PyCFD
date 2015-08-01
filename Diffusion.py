'''
    One-dimensional diffusion problem solved by using a FTCS scheme.
'''
from __future__ import division
import ROOT as rt
import numpy as np
import time as tm
from math import *

### SIMULATION VARIABLES

dx = 1e-1
dt = 5e-1

t0 = 0.
t1 = 500.

x0 = 0.
x1 = 4.

nx = int( (x1-x0) / dx )
nt = int( (t1-t0) / dt )

alpha = 1e-2
s1  = alpha * dt / dx**2
s2 = 1 - 2*s1 

assert s1 < 1.,'Nope. S1 = '+ str(s1)

data = np.zeros( (nt,nx) )

###### BOUNDARY CONDITIONS
T0 = 100.

data[0,:] = 0.

data[:,0] = T0
data[:,-1] = T0

### SIMULATION
for t in range(nt-1):
    for x in range(1,nx-1):
        data[t+1,x] = s2 * data[t,x] + s1 * ( data[t,x-1] + data[t,x+1] )

### PLOTTING
cv = rt.TCanvas()
hs = [ rt.TH1F( str(i),';x;T', nx, x0, x1 ) for i in range(nt) ]
for h in hs: h.SetMinimum(0.)
for h in hs: h.SetMaximum(1.1*T0)
for t in range(nt):
    for x in range(nx):
        hs[t].SetBinContent(x+1,data[t,x])
    hs[t].Draw()
    cv.Update()
#    raw_input()
    tm.sleep(0.05)


    

