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
plt.figure()
#plt.show()

def error(function,f,t0,tf,U0,R,xbar,ybar):
    errlist = []
    N = 2
    i = 0
    plotlist = []
    while i < 20:
        t, U = function(f, t0, tf, U0, N)
        U = np.array(U)
        X = U[:,0]
        Y = U[:,1]
        err = np.abs(np.sqrt((X-xbar)**2 + (Y-ybar)**2) - R)
        # print(err)
        plotlist.append(N)
        errmax = np.max(err)
        errlist.append(errmax)
        # print(errmax)
        N = N*2
        i += 1
    errvec = np.array(errlist)
    kvot = (errvec[:-1]/errvec[1:])
    plt.loglog(plotlist,errlist)
    plt.axis('equal')
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
    while i < N:
        k1 = f(t, g)
        k2 = f(t+h/2, g + h*k1/2)
        k3 = f(t+h/2, g + h*k2/2)
        k4 = f(t + h, g + h*k3)
        #print(k1,k2,k3,k4)
        # temp = (k1+2*k2+2*k3+k4)
        gnew = g + h*(k1+k2*2+k3*2+k4)/6


        #print('här')
        #print(gnew)
        Ulist.append(g)
        g = gnew
        tlist.append(t)
        t += h
        i += 1
    return tlist, Ulist

error(rk4,f,t0,tf,U0,R,xbar,ybar)

plt.figure()

N = 1000
wl =4
wr = 2
b = 1
U = [x,y,o]
U0 = [0,0,0]
#U =[]
t0=0
t, U = rk4(f,t0,tf,U0,N)
U = np.array(U)
plt.plot(U[:,0],U[:,1])
plt.axis('equal')
# plt.show()

# n = [10,20,40,80,160,320,640,1280]

# konv1 = []
# konv2 = []

# for steg in n:
#     t1, U1 = forward_euler(f,t0,tf,U0,steg)
#     t2, U2 = rk4(f,t0,tf,U0,steg)

#     eX = U1[0]
#     eY = U1[1]
#     eT = U1[2]

#     rX = U2[0]
#     rY = U2[1]
#     rT = U2[2]

#     Efel = []
#     rFel = []

#     for i in range(len(eX)):
#         Efel.append(np.abs(np.sqrt(eX[i]-xbar)**2 + (eY[i]-ybar)**2)-R)
#         rFel.append(np.abs(np.sqrt(rX[i]-xbar)**2 + (rY[i]-ybar)**2)-R)

#     konv1.append(max(Efel))
#     konv2.append(max(rFel))

# print(konv1)
# print(konv2)
# plt.loglog(n, konv1)
# plt.loglog(n, konv2)
plt.show()