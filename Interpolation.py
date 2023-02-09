import math

from numpy import double
import numpy as np
import matplotlib.pyplot as plt
from sympy import *

# e**x+2**(-x)+2*cos(x)-6
print("Enter e**x+2**(-x)+2*cos(x)-6 as function")
func = sympify(input("Enter the function: "))
w = 0
while (w == 0):
    UB = double(input("Enter the upper bound"))
    LB = double(input("Enter the lower bound"))
    if (UB > LB):
        w = 1
    else:
        print("Upper Bound should be greater than Lower Bound")

n = int(input("Enter the number of points 'n' to be taken"))
e = symbols('e')
x = symbols('x')

Ny = []
Nx = []

step = (UB-LB)/(n-1)
for i in range(n):
    a = LB+step*i
    Nx.append(a)
    Ny.append(func.evalf(subs={x:a,e:math.e}))

An = Ny
NumberofPoints = 100
CSxStep = (UB-LB)/(NumberofPoints-1)
CSx = []
for i in range(NumberofPoints):
    a = LB+CSxStep*i
    CSx.append(a)

Px = []
Py = []
An = Ny


step = (UB - LB) / (NumberofPoints-1)
for i in range(NumberofPoints):
    a = LB + step * i
    Px.append(a)
    Py.append(func.evalf(subs={x: a, e: math.e}))

Lx = []
for i in range(100):
    a = LB + step * i
    Lx.append(a)

print(Ny)
print(Nx)



def ARowMaker(func,UB,LB,n,R):
    h = (UB-LB)/(n-1)
    if(R==1):
        Row = np.array([1])
        for i in range(n-1):
            Row = np.append([Row],[0])

    if(R==n):
        Row = np.array([0])
        for i in range(n-2):
            Row = np.append([Row],[0])
        Row = np.append([Row],[1])


    if(R!=n and R!=1):
        Row = np.array([0])
        for i in range(n-1):
            Row = np.append([Row],[0])
        Row[R-1] = 4*h
        Row[(R-2)] = h
        Row[(R)] = h
    print(Row)
    return Row

def AMaker(func,UB,LB,n):
    A = ARowMaker(func,UB,LB,n,1)
    for i in range(n-1):
        A = np.vstack([A,ARowMaker(func,UB,LB,n,i+2)])
    return A


def bMaker(func,UB,LB,n):
    h = (UB - LB) / (n - 1)
    b = np.array(0)
    for i in range(n-2):
        T = ((3)/h)*(Ny[i+2]-2*Ny[i+1]+Ny[i])
        b = np.append(b,T)
    b = np.append(b,0)
    print(b.transpose())
    return b.transpose()


def Cmaker(func,UB,LB,n):
    A = AMaker(func,UB,LB,n)
    b = bMaker(func,UB,LB,n)
    Ainv = np.linalg.inv(A)
    Cn = np.matmul(Ainv,b)
    return Cn

def Bmaker(Cn,An,UB,LB,n):
    h = (UB-LB)/(n-1)
    i = 0
    T = (An[i + 1] - An[i]) / h - (h / (3)) * (2 * Cn[i] + Cn[i + 1])
    Bn = np.array(T)
    for i in range(len(An)-1):
        i=i+1
        if(i+1>=len(An)):
            T = (An[0] - An[i]) / h - (h / (3)) * (2 * Cn[i] + Cn[0])
        else:
            T = (An[i+1]-An[i])/h - (h/(3))*(2*Cn[i]+Cn[i+1])
        Bn = np.append(Bn,T)
    return Bn

def Dmaker(Cn,UB,LB,n):
    h = (UB-LB)/(n-1)
    i = 0
    T = T = (Cn[i+1]-Cn[i])/((3)*h)
    Dn = np.array(T)
    for i in range(len(An)-1):
        i=i+1
        if(i+1>=len(An)):
            T = (Cn[0] - Cn[i]) / ((3) * h)
        else:
            T = (Cn[i+1]-Cn[i])/((3)*h)
        Dn = np.append(Dn,T)
    return Dn

def CubicSplineEqs(func,An,UB,LB,n):

    Cn = Cmaker(func,UB,LB,n)
    Bn = Bmaker(Cn,An,UB,LB,n)
    Dn = Dmaker(Cn,UB,LB,n)


    print(Nx)
    Eqs = []
    for i in range(n-1):

      S = Dn[i]*((x-Nx[i])**3)+Cn[i]*((x-Nx[i])**2)+Bn[i]*(x-Nx[i])+An[i]
      Eqs = np.append(Eqs,S)

    for i in range(len(Eqs)):
        print(Eqs[i])
    return Eqs


def CSyMaker(func,An,UB,LB,n,CSx):
    Eqs = CubicSplineEqs(func, An, UB, LB, n)


    CSy = []
    i2 = 0
    for i in range(len(CSx)):
        if(Nx[i2]<=CSx[i] and CSx[i]<=Nx[i2+1]):
            CSy= np.append(CSy,sympify(Eqs[i2]).evalf(subs={x:CSx[i],e:math.e}))

        else:
            #print(str(i2) + " "+ str(i))
            i2 = i2 + 1
            CSy = np.append(CSy, sympify(Eqs[i2]).evalf(subs={x: CSx[i], e: math.e}))

    return CSy

def nmaker(n, N):
    while (n > N):
        n = n - N
    return n


def Lagrange(X, Y):
    L = sympify("0")
    for i in range(len(X)):
        L = L + LagTemp(X, Y, i)
    return L


def LagTemp(X, Y, n):
    LagrangeN = sympify("1")
    LagrangeN = LagrangeN * Y[n]
    for i in range(len(X)):
        if (i != n):
            LagrangeN = LagrangeN * (x - X[i])

    LagrangeD = sympify("1")
    for i in range(len(X)):
        if (i != n):
            LagrangeD = LagrangeD * (X[n] - X[i])

    # print("Numerator: ")
    # print(LagrangeN)
    # print("Denominator: ")
    # print(LagrangeD)
    # print("Lagrange: ")
    # print(LagrangeN/LagrangeD)
    return LagrangeN / LagrangeD


def LagVals(X, Y):
    Ly = []
    for i in range(len(Lx)):
        Ly.append(Lagrange(X, Y).evalf(subs={x: Lx[i], e: math.e}))
    return Ly

def LenMaker(V,Len):

    L = len(str(V))
    d = Len-L
    if(L>Len):
        print("X")
    for i in range (d):
        V = str(V) + " "

    return V
def Write(F,Nx):
    F = F.tolist()
    A = len(F)
    print(A)
    L = (int)((math.sqrt(1+8*A)-1)/2)
    gap = "                         " #25 spaces
    File = open("DividedDifference.txt", "w")

    xline = "x     "
    for i in range(len(Nx)):
        xline = xline + LenMaker(Nx[i],25) + gap
    File.write(xline + "\n")

    for i in range(L):
        if(i==0):
            start = "f(x) "
        else:
            start = "f" + str(i) + "(x)"
        st = ""
        for j in range(L-i):
            V = F.pop(0)
            st =  st + LenMaker(str(V),25) + gap

        for j in range(i):
            st = gap + st
        print(st)
        File.write(start +" "+st+ "\n" + "\n")
    File.close()

def DividedDiffTemp(Nx,Fn):
    d = len(Nx)-len(Fn) +1
    T = np.array([])
    for i in range(len(Fn)-1):

        T = np.append(T,(Fn[i]-Fn[i+1])/(Nx[i]-Nx[i+d]))
    return T

def DividedDiff(Nx ,Ny):


    File = open("DividedDifference.txt", "a")
    X = np.array(Nx)
    F = np.array(Ny)
    Fn = Ny
    for i in range(len(Nx)-1):
        Fn = DividedDiffTemp(Nx,Fn)

        F = np.append(F,Fn)
    Write(F,Nx)
    return 0

Ly = LagVals(Nx, Ny)
print(Ly)

Ey = []
for i in range(len(Px)):
    Ey.append(Py[i] - Ly[i])

print(func)
print(len(An))
print(UB)
print(LB)
print(n)
print(len(CSx))

CSy = CSyMaker(func,An,UB,LB,n,CSx)
DividedDiff(Nx,Ny)

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)

ax1.set_title("f(x)")
ax2.set_title("Lagrange")
ax3.set_title("Cublic Spline")
ax4.set_title("f(x) vs L(x) vs CS(x)")
ax1.plot(Px, Py, color='red', label='f(x)', linewidth=0.5)  # plot f(x)

ax2.plot(Lx, Ly, color='blue', label='L(x)', linewidth=0.5)  # plot lagrange polynomial

ax3.plot(CSx, CSy, color='green', label='CS(x)', linewidth=0.5)  # plot Cubic Spline
ax4.plot(Px, Py, color='red', label='f(x)', linewidth=0.5)  # plot f(x)

ax4.plot(Lx, Ly, color='blue', label='L(x)', linewidth=0.5)  # plot lagrange polynomial

ax4.plot(CSx, CSy, color='green', label='CS(x)', linewidth=0.5)  # plot Cubic Spline
plt.legend()
plt.show()


