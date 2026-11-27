from collections import UserList

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

#Uppgift 2
def robotcircle(b,wl,wr,x0,y0,o0):# Där o0 = theta0
    Dw = -(wr-wl)/b
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
    A = 0
    B = 0
    while i < N:
        Ulist.append(g)
        tlist.append(t)
        gnew = g + h*f(t, g)
        g = gnew
        t += h
        i += 1
    Ulist.append(g)
    tlist.append(t)
    return tlist, Ulist, A, B
def x():
    return x
def y():
    return y
def o():
    return o
N = 15
wl =4
wr = 2
b = 1
U = [x,y,o]
U0 = np.array([0.0,0.0,0.0])
#U =[]
t0=0
xbar,ybar,R,tf=robotcircle(b,wl,wr,U0[0],U0[1],U0[2])
# print("här")
# print(xbar)
# print(ybar)
# print(R)
# print(tf)
vl = lambda t: wl
vr = lambda t: wr
f = lambda t, U: odefun(t,U,b,vl,vr)
t, U, A,B = forward_euler(f,t0,tf,U0,N)
U = np.array(U)
plt.plot(U[:,0],U[:,1])
plt.axis('equal')
#plt.show()

def error(function,f,t0,tf,U0,R,xbar,ybar):
    errlist = []
    N = 2
    i = 0
    plotlist = []
    while i < 10:
        t, U, A,B = function(f, t0, tf, U0, N)
        # print(".")
        # print(N)
        # print(t0)
        # print(tf)
        # print(U0)
        # print(len(t))
        # print(U)
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
    kvot = np.log((errvec[:-1]/errvec[1:]))/np.log(2)
    # kvot = (errvec[:-1]/errvec[1:])
    # plt.loglog(plotlist,errlist)
    # plt.axis('equal')
    print(kvot)
# print("här")
# print(xbar)
# print(ybar)
# print(R)
# print(tf)
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
        Ulist.append(g)
        tlist.append(t)
        k1 = f(t, g)
        k2 = f(t+h/2, g + h*k1/2)
        k3 = f(t+h/2, g + h*k2/2)
        k4 = f(t + h, g + h*k3)
        #print(k1,k2,k3,k4)
        # temp = (k1+2*k2+2*k3+k4)
        gnew = g + h*(k1+k2*2+k3*2+k4)/6
        #print('här')
        #print(gnew)
        g = gnew
        t += h
        i += 1
    
    Ulist.append(g)
    tlist.append(t)
    A = g[0]
    B = g[1]
    return tlist, Ulist, A, B
# print("här")
# print(xbar)
# print(ybar)
# print(R)
# print(tf)
error(rk4,f,t0,tf,U0,R,xbar,ybar)

# N = 15
# wl =4
# wr = 2
# b = 1
# U = [x,y,o]
# # U0 = [0,0,0]
# U0 = np.array([0.0,0.0,0.0])
# #U =[]
# t0=0
# t, U, diff = rk4(f,t0,tf,U0,N)
# U = np.array(U)
# plt.plot(U[:,0],U[:,1])
# plt.axis('equal')



"""Uppgift 4"""

def bana(f):
    N = 5
    tf = 10
    b = 1
    U0 =[0,0,0]
    tol = 1e-5
    diff = 1.0

    t, U1, X1, Y1 = rk4(f,t0,tf,U0,N)
    i = 0
    while diff > tol:
        N = N*2
        X2 = X1
        Y2 = Y1
        t, U1, X1, Y1 = rk4(f,t0,tf,U0,N)
        deltaX = np.abs(X1-X2)
        deltaY = np.abs(Y1-Y2)
        if deltaY < deltaX:
            diff = deltaX
        else:
            diff = deltaY
        # print(N)
        # print(diff)
        U1 = np.array(U1)
        i += 1
    plt.figure()
    plt.plot(U1[:,0],U1[:,1])
    plt.axis('equal')

###a)
print("a")
vl = lambda t: 3*t + 5
vr = lambda t: 3*t + 1
f = lambda t, U: odefun(t,U,b,vl,vr)

bana(f)

###b)
print("b")
vl = lambda t: 3*t + 1
vr = lambda t: t + 6
f = lambda t, U: odefun(t,U,b,vl,vr)

bana(f)

###c)
print("c")
k = 2.7
vl = lambda t: k*(np.cos(k*t)-0.75)
vr = lambda t: k*(np.cos(k*t)+0.5)
f = lambda t, U: odefun(t,U,b,vl,vr)

bana(f)

###d)
print("d")
k=0.75
vl = lambda t: k*(0.5-5*np.cos(5*k*t))
vr = lambda t: k*(1.5+5*np.cos(5*k*t))
f = lambda t, U: odefun(t,U,b,vl,vr)

bana(f)



"""Uppgift 5"""

# def interp(X, Y, u):
#     X = np.array(X)
#     Y = np.array(Y)
#     A = []
#     j = 0
#     for i in X:
#         A.append(x**j)
#     e = np.ones(X.shape)




a = 1
def uppg5(a):
    N = 1000
    J0 = sp.special.jv(0,2*a)
    g = lambda t: a*np.sin(3*J0*t)
    tf = 30
    vl = lambda t: (np.cos(g(t)))**2
    vr = lambda t: (np.sin(g(t)))**2 
    f = lambda t, U: odefun(t,U,b,vl,vr)

    t, U, X1, Y1 = rk4(f,t0,tf,U0,N)

    plt.figure()
    U = np.array(U)
    plt.plot(U[:,0],U[:,1])
    plt.axis('equal')

    ###a)
    N = 4
    tol = 1e-6
    diff = 1
    Tn = 0
    while diff > tol:
        T = 0
        j = 0
        t,U,X1,Y1 = rk4(f,t0,tf,U0,N)
        U = np.array(U)
        templist = U[:,0]
        Tg = Tn
        for i in templist:
            if i>=0:
                j+=1
            elif i<0 and T==0:
                T = np.interp(0, [templist[j], templist[j-1]],
                                [t[j], t[j-1]])
                Tn = T

        if Tn != 0:
            diff = np.abs(Tn - Tg)

        N=N*2
        Tk = 2*Tn
        # print(Tk,diff)
    print("Omloppstiden är:", Tk)
    ###b)

    tol = 1e-6
    N = 2
    area = 0
    diff = 1
    while diff > tol:
        Ag = area
        t, U, X1, Y1 = rk4(f, t0, Tk, U0, N)
        U = np.array(U)
        x = U[:, 0]
        y = U[:, 1]
        o = U[:, 2]
        i = 0
        A = 0
        h = (Tk - t0) / N
        for i in range(N + 1):
            af = 1 / 4 * (x[i] * np.sin(o[i]) - y[i] * np.cos(o[i]))
            if i == 0 or i == N:
                A += af / 2
            else:
                A += af

        area = np.abs(h * A)
        diff = np.abs(area-Ag)
        N = N*2
        # print(area,"ae",diff)
    print("Arean är:", area)

uppg5(a)

"""Uppgift 6"""



def compute_area(a, N):
    t0 = 0
    tf = 100
    J0 = sp.special.jv(0,2*a)
    g = lambda t: a*np.sin(3*J0*t)
    vl = lambda t: (np.cos(g(t)))**2
    vr = lambda t: (np.sin(g(t)))**2 
    f = lambda t, U: odefun(t,U,b,vl,vr)

    n = N
    
    tol = 1e-6
    diff = 1
    Tn = 0
    Tg = 0
    while diff > tol:
        T = 0
        j = 0
        t,U,X1,Y1 = rk4(f,t0,tf,U0,N)
        U = np.array(U)
        templist = U[:,0]
        
        Tg = Tn
        for i in templist:
            if i >= 0:
                j+=1
            elif i<0 and T==0:
                T = np.interp(0, [templist[j], templist[j-1]],
                                [t[j], t[j-1]])
                Tn = T
        if Tn != 0:
            diff = np.abs(Tn - Tg)
        N=N*2
        Tk = 2*Tn

    tol = 1e-6
    N = n
    area = 0
    diff = 1
    while diff > tol:
        Ag = area
        t, U, X1, Y1 = rk4(f, t0, Tk, U0, N)
        U = np.array(U)
        x = U[:, 0]
        y = U[:, 1]
        o = U[:, 2]
        i = 0
        A = 0
        h = (Tk - t0) / N
        for i in range(N + 1):
            af = 1 / 4 * (x[i] * np.sin(o[i]) - y[i] * np.cos(o[i]))
            if i == 0 or i == N:
                A += af / 2
            else:
                A += af

        area = np.abs(h * A)
        diff = np.abs(area-Ag)
        N = N*2

    return area

# a = np.linspace(0,1.131,100)
# pltlist = []
# pltlist1 = []
# pltlist2 = []
# pltlist3 = []
# pltlist4 = []
# pltlist5 = []

# for i in a:
#     pltlist.append(compute_area(i, 1)-20)
#     pltlist1.append(compute_area(i, 2)-20)
#     pltlist2.append(compute_area(i, 3)-20)
#     pltlist3.append(compute_area(i, 4)-20)
#     pltlist4.append(compute_area(i, 5)-20)
#     pltlist5.append(compute_area(i, 6)-20)

# plt.figure()
# plt.grid()
# plt.plot(a,pltlist)
# plt.plot(a,pltlist1)
# plt.plot(a,pltlist2)
# plt.plot(a,pltlist3)
# plt.plot(a,pltlist4)
# plt.plot(a,pltlist5)

# plt.show()

N = 6

def sekant(g0,g1):
    tol = 1e-7
    diff = 1
    while diff > tol:
        f1 = compute_area(g1, N) - 20
        # print(f1)
        f0 = compute_area(g0, N) - 20
        # print(f0)
        g2 = g1 - f1*((g1-g0)/(f1-f0))
        print(g2)

        diff = np.abs(g2 - g1)
        g0 = g1
        g1 = g2
        # print ("sek",g2, diff)
    print("a = ", g2)
    print("Arean är",compute_area(g2, N))
    return g2

print("1")
a1 = sekant(1.03,1.04)

print("2")
a2 = sekant(1.117,1.121)

# print("3")
# sekant(1.127,1.219)
uppg5(a1)
uppg5(a2)
uppg5(1)