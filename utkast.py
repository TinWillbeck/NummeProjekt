from collections import UserList

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
    return xbar, ybar, np.abs(R), np.abs(Tlap)

#Uppgift 3
#b
#
# U = x,y,theta
#
#


def odefun(t,U,b,vL,vR):
    v1 = (vR(t) + vL(t))*np.cos(U[2])/2
    v2 = (vR(t) + vL(t))*np.sin(U[2])/2
    v3 = (vL(t) - vR(t))/b
    return np.array([v1,v2,v3])

def forward_euler(f,t0,tf,U0,N):
    tlist = []
    Ulist = []
    h = (tf-t0)/N
    g = U0
    t = t0
    i = 0
    while i < N:
        gnew = g + h*f(t, g)
        Ulist.append(g)
        g = gnew
        tlist.append(t)
        t += h
        i += 1
    return tlist, Ulist
def x():
    return x
def y():
    return y
def o():
    return o
N = 1000
wl =4
wr = 2
b = 1
U = [x,y,o]
U0 = [0,0,0]
#U =[]
t0=0
xbar,ybar,R,tf=robotcircle(b,wl,wr,U0[0],U0[1],U0[2])
vl = lambda t: wl
vr = lambda t: wr
f = lambda t, U: odefun(t,U,b,vl,vr)
t, U = forward_euler(f,t0,tf,U0,N)
U = np.array(U)
plt.plot(U[:,0],U[:,1])
plt.axis('equal')
#plt.show()

def error(function,f,t0,tf,U0,R,xbar,ybar):
    errlist = []
    N = 1
    i = 0
    while i < 20:
        t, U = function(f, t0, tf, U0, N)
        U = np.array(U)
        X = U[:,0]
        Y = U[:,1]
        err = np.abs(np.sqrt((X-xbar)**2 + (Y-ybar)**2) - R)
        N = N*2
        errmax = np.max(err)
        errlist.append(errmax)
        i += 1
    errvec = np.array(errlist)
    kvot = (errvec[:-1]/errvec[1:])
    print(kvot)
error(forward_euler,f,t0,tf,U0,R,xbar,ybar)

#Uppgift 3e Runge kutta
def rk4(f,t0,tf,U0,N):
    tlist = []
    Ulist = []
    h = (tf - t0) / N
    g = U0
    t = t0
    i = 0
    gnew = 0
    while i < N:
        k1 = f(t, g)
        k2 = f(t+h/2, g + h*k1/2)
        k3 = f(t+h/2, g + h*k2/2)
        k4 = f(t + h, g + h*k3)
        #print(k1,k2,k3,k4)
        temp = (k1+2*k2+2*k3+k4)
        gnew = g + h*temp/6


        #print('här')
        #print(gnew)
        Ulist.append(g)
        g = gnew
        tlist.append(t)
        t += h
        i += 1
    return tlist, Ulist

error(rk4,f,t0,tf,U0,R,xbar,ybar)