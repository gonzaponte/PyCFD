'''
    One-dimensional diffusion problem solved by using a FTCS scheme.
'''
from __future__ import division
import ROOT as rt
import numpy as np
import time as tm
from math import *

### SIMULATION VARIABLES

dx = 2e-2
dt = 1e-1

t0 = 0.
t1 = 50.

x0 = 0.
x1 = 1.

nx = int( (x1-x0) / dx )
nt = int( (t1-t0) / dt )

alpha = 1e-4
s1  = alpha * dt / dx**2
s2 = 1 - 2*s1 

assert s1 < 1.,'Nope. S1 = '+ str(s1)

data = np.zeros( (nt,nx) )
true = np.zeros( (nt,nx) )

###### BOUNDARY CONDITIONS
T0 = 100.

data[0,:]  = true[0,:] = 0.

data[:,0]  = true[:,0] = T0
data[:,-1] = true[:,-1] = T0

### SIMULATION
for t in range(nt-1):
    tr = (t+1) * dt 
    for x in range(1,nx-1):
        xn = x / ( nx - 1 ) 
        data[t+1,x] = s2 * data[t,x] + s1 * ( data[t,x-1] + data[t,x+1] )
        true[t+1,x] = 100 - sum( 400 / pi / ( 2*m+1 ) * sin( (2*m+1)*pi*xn ) * exp(-alpha*pi**2*tr*(2*m+1)**2) for m in range(20) )

### PLOTTING
cv = rt.TCanvas()
hdata = [ rt.TH1F( 'data_' + str(i),';x;T', nx, x0, x1 ) for i in range(nt) ]
htrue = [ rt.TH1F( 'true_' + str(i),';x;T', nx, x0, x1 ) for i in range(nt) ]

for h in hdata: h.SetMinimum(0.)
for h in hdata: h.SetMaximum(1.1*T0)
for h in htrue: h.SetMinimum(0.)
for h in htrue: h.SetMaximum(1.1*T0)
for h in htrue: h.SetLineColor(rt.kRed)

for t in range(nt):
    for x in range(nx):
        hdata[t].SetBinContent(x+1,data[t,x])
        htrue[t].SetBinContent(x+1,true[t,x])
    hdata[t].Draw()
    htrue[t].Draw('same')
    cv.Update()
#    raw_input()
    tm.sleep(0.05)


    

