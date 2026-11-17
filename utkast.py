import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

#Uppgift 2
def robotcircle(b,wl,wr,x0,y0,o0):# Där o0 = theta0
    Dw = (wr-wl)/b
    w = (wr+wl)/2
    xbar = x0 - (w/Dw)*np.sin(o0)
    ybar = y0 + (w/Dw)*np.cos(o0)
    R = w/Dw
    Tlap = 2*np.pi*b/(wr-wl)
    return xbar, ybar, R, Tlap

#Uppgift 3
#b

def odefun(t,U,b,vL,vR):
    v1 = (vR+vL)*np.cos(U[0])/2
    v2 = (vR+vL)*np.sin(U[0])/2
    v3 = (vL-vR)/b
    return np.array([v1,v2,v3])

def forward_euler(f,t0,tf,U0,N):
    tlist = []
    Ulist = []
    h = (tf-t0)/N
    g = U0
    t = t0
    i = 0
    while i < N:
        gnew = g + h*f(t,U)
        Ulist.append(U)
        g = gnew
        tlist.append(t)
        t += h
        i += 1
    return tlist, Ulist
N = 1000
wl =4
wr = 2
b = 1
U0 = [0,0,0]
t0=0
tf=robotcircle(b,wl,wr,U0[0],U0[1],U0[2])[3]
vL = lambda t: wl
vR = lambda t: wr
f = lambda t, U: odefun(t,U,b,vL,vR)
t, U = forward_euler(f,t0,tf,U0,N)





