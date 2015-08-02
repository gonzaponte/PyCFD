'''
    One-dimensional diffusion problem solved by using a FTCS scheme.
'''
from __future__ import division
import ROOT as rt
import numpy as np
import time as tm
from math import *

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


    

