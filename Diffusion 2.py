'''
    One-dimensional diffusion problem.
'''
from __future__ import division
import ROOT as rt
import numpy as np
import time as tm
from math import *

def modifiable():
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
#        raw_input()
        tm.sleep(0.05)

def FCTS( alpha, s, dx, tmax, T0 = 100, draw = False ):
    '''
        Solves the diffusion equation with a FCTS scheme.
        This scheme is stable for s<= 0.5 and optimal
        for s = 1/6. The simulation variables are:
        - dx => spatial step
        - dt => time step
        - tmax => maximum time
        - alpha => diffusion coefficient
        - s => stability variable
        - T0 => temperature in the extremes
        - nx => number of points in space
        - data => data at time level n
    '''
    assert s <= 0.5,'Not stable. s must be less or equal than 0.5 and it is '+ str(s)

    ###### DEPENDENT VARIABLES
    dt = s * dx**2 / alpha
    ss = 1 - 2*s
    nx = int( 1 / dx ) + 1
    data = np.array([0.] * nx)

    pstr  = ['Solving diffusion problem using a FCTS scheme. The input parameters are:']
    pstr += ['dx = ' + str(dx)]
    pstr += ['dt = ' + str(dt)]
    pstr += ['tmax = ' + str(tmax)]
    pstr += ['s = ' + str(s)]
    pstr += ['alpha = ' + str(alpha)]

    print '\n'.join(pstr)

    ###### BOUNDARY CONDITIONS
    data[1:-1] = 0. # T = 0 for t = 0
    data[0] = data[-1] = T0 # T = T0 for the boundary

    ###### PLOTTING
    if draw:
        cv = rt.TCanvas()
        hdata = rt.TH1F( 'data',';x;T', nx, 0, 1.1 )
        hdata.SetMaximum(1.1*T0)
        hdata.SetMinimum(0)
        hdata.SetLineColor(rt.kBlack)

    t = 0.
    ###### SIMULATION
    while t<tmax:
        data[1:-1] = [ ss * data[x] + s * ( data[x-1] + data[x+1] ) for x in range(1,nx-1) ]
        if draw:
            map( hdata.SetBinContent,range(1,nx+1),data )
            hdata.Draw()
            cv.Update()
            tm.sleep(0.05)
        t += dt

def DuFort( alpha, s, dx, tmax, T0 = 100, draw = False ):
    '''
        Solves the diffusion equation a DuFort-Frankel scheme.
        This scheme is unconditionally stable and optimal
        for s = sqrt(1/12). The simulation variables are:
        - dx => spatial step
        - dt => time step
        - tmax => maximum time
        - alpha => diffusion coefficient
        - s => stability variable
        - T0 => temperature in the extremes
        - nx => number of points in space
        - data => data at time level n
    '''

    ###### DEPENDENT VARIABLES
    dt = s * dx**2 / alpha
    sp = 2*s/(1+2*s)
    sm = (1-2*s)/(1+2*s)
    nx = int( 1 / dx ) + 1
    data0 = np.array([0.] * nx)
    data1 = np.array([0.] * nx)

    pstr  = ['Solving diffusion problem using a DuFort-Frankel scheme. The input parameters are:']
    pstr += ['dx = ' + str(dx)]
    pstr += ['dt = ' + str(dt)]
    pstr += ['tmax = ' + str(tmax)]
    pstr += ['s = ' + str(s)]
    pstr += ['alpha = ' + str(alpha)]

    print '\n'.join(pstr)

    ###### BOUNDARY CONDITIONS
    data0[1:-1] = data1[1:-1] = 0. # T = 0 for t = 0
    data0[0] = data0[-1] = T0 # T = T0 for the boundary
    data1[0] = data1[-1] = T0 # T = T0 for the boundary

    ###### PLOTTING
    if draw:
        cv = rt.TCanvas()
        hdata = rt.TH1F( 'data',';x;T', nx, 0, 1.1 )
        hdata.SetMaximum(1.1*T0)
        hdata.SetMinimum(0)
        hdata.SetLineColor(rt.kBlack)

    t = 0.
    ###### SIMULATION
    while t<tmax:
        data1[1:-1], data0 = [ sp * (data1[x-1] + data1[x+1]) + sm * data0[x] for x in range(1,nx-1) ], data1[:]
        if draw:
            map( hdata.SetBinContent,range(1,nx+1),data1 )
            hdata.Draw()
            cv.Update()
            tm.sleep(0.05)
        t += dt

def CrankNicholson( alpha, s, dx, tmax, T0 = 100, draw = False ):
    '''
        Solves the diffusion equation a DuFort-Frankel scheme.
        This scheme is unconditionally stable and optimal
        for s = sqrt(1/12). The simulation variables are:
        - dx => spatial step
        - dt => time step
        - tmax => maximum time
        - alpha => diffusion coefficient
        - s => stability variable
        - T0 => temperature in the extremes
        - nx => number of points in space
        - data => data at time level n
    '''

    ###### DEPENDENT VARIABLES
    dt = s * dx**2 / alpha
    sp = 2*s/(1+2*s)
    sm = (1-2*s)/(1+2*s)
    nx = int( 1 / dx ) + 1
    data0 = np.array([0.] * nx)
    data1 = np.array([0.] * nx)

    pstr  = ['Solving diffusion problem using a DuFort-Frankel scheme. The input parameters are:']
    pstr += ['dx = ' + str(dx)]
    pstr += ['dt = ' + str(dt)]
    pstr += ['tmax = ' + str(tmax)]
    pstr += ['s = ' + str(s)]
    pstr += ['alpha = ' + str(alpha)]

    print '\n'.join(pstr)

    ###### BOUNDARY CONDITIONS
    data0[1:-1] = data1[1:-1] = 0. # T = 0 for t = 0
    data0[0] = data0[-1] = T0 # T = T0 for the boundary
    data1[0] = data1[-1] = T0 # T = T0 for the boundary

    ###### PLOTTING
    if draw:
        cv = rt.TCanvas()
        hdata = rt.TH1F( 'data',';x;T', nx, 0, 1.1 )
        hdata.SetMaximum(1.1*T0)
        hdata.SetMinimum(0)
        hdata.SetLineColor(rt.kBlack)

    t = 0.
    ###### SIMULATION
    while t<tmax:
        data1[1:-1], data0 = [ sp * (data1[x-1] + data1[x+1]) + sm * data0[x] for x in range(1,nx-1) ], data1[:]
        if draw:
            map( hdata.SetBinContent,range(1,nx+1),data1 )
            hdata.Draw()
            cv.Update()
            tm.sleep(0.05)
        t += dt



#FCTS(1e-2,1/6, 8e-2,1e4,1000,True)
#DuFort(1e-2,1/12**0.5, 8e-2,1e4,1000,True)
