import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.legendre import Legendre
from numpy.polynomial.chebyshev import Chebyshev
import scipy.linalg as la
import scipy.integrate as integrate
from pytex import *

#from hw4 import polyinterp

def p1():
    cond = []
    est = []
    ratio = []
    for n in range(5,20):
        H = la.hilbert(n)
        K = np.linalg.cond(H,p=np.inf)
        cond.append(K)
        E = (1+np.sqrt(2))**(4*n) / np.sqrt(n)
        est.append( E )
        ratio.append (K/E)
    latex_table((range(5,len(cond)+5), cond, est, ratio), ('$n$','$K(H_n)$', 'est', '$K(H_n)/$est'))
        
        
def p4():
    deg = 4
    c = np.zeros(deg+1)
    dom = [-np.pi, np.pi]
    for i in range(0,5):
        coef = np.zeros(deg+1)
        coef[i] = 1
        prod = lambda x: Legendre(coef, dom)(x)**2
        den = integrate.quad(prod, -np.pi, np.pi)
        prod = lambda x: Legendre(coef, dom)(x)*np.cos(x)
        num = integrate.quad(prod, -np.pi, np.pi)
        c[i] = num[0]/den[0]
    p = Legendre(c, dom)
    #plot them
    x = np.linspace(-np.pi, np.pi,500)
    yp = p(x)
    ycos = np.cos(x)
    fig, (ax1, ax2) = plt.subplots(2,1, sharex=True)
    ax1.plot(x,ycos, '-k')
    ax1.plot(x,yp, '-c')
    ax2.plot(x, ycos-yp,'-k')
    plt.show()
    

def taylor_series(x):
    sum = 0
    for i in range(11):
        sum += (-1)**i /(2*i+1) * x**(2*i+1)
    return sum

def cheb_series(x):
    sum = 0
    alpha = np.sqrt(2)-1
    coef = np.zeros(11*2+1)
    for i in range(11):
        coef[2*i+1] = 1
        sum += (-1)**i /(2*i+1) * alpha**(2*i+1) * Chebyshev(coef)(x)
        coef[2*i+1] = 0
    return 2*sum
    
def p5b():
    x = np.linspace(-1, 1, 1000)
    y = np.arctan(x)
    yt = taylor_series(x)
    yc = cheb_series(x)
    plt.plot(x, y-yt, 'r-')
    plt.plot(x, y-yc, 'b-')
    plt.show()
    print('Max error for Taylor Series: ', np.max(y-yt))
    print('Max error for Chebyshev Series: ', np.max(y-yc))
    
    
def polyinterp(u, x, y, w=None):
    if w == None:
        w = baryweights(x)

    ret = np.zeros(len(u))
    for i in range(len(ret)):
        if u[i] in x:
            ret[i] = y[np.where(x==u[i])]
        else:
            weights = w /(u[i] - x)
            ret[i] = weights.dot(y)/sum(weights)
    return ret
    
def baryweights(x):
    w = np.ones(len(x))
    for j, xj in enumerate(x):
        for i, xi in enumerate(x): 
            if j != i:
                w[j] /= (xj - xi)
    return w
    
def cheb_points2(n, a, b):
    x = []
    for j in range(0,n):
        x.append( -1* np.cos( (j*np.pi)/(n-1) ) )
    return x
    
def p5c():
    x = cheb_points2(21, -1, 1)
    y = np.arctan(x)
    
    u = np.linspace(-1, 1, 1000)
    p = polyinterp(u,x,y)
    
    f = np.arctan(u)
    
    print('Max error: ', np.max(np.abs(p-f)))
    
    plt.plot(u, f-p, 'r-')
    plt.show()
        
        
        
        