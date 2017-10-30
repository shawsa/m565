import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import CubicSpline

def p1():
    plt.plot( (-100, 100), (0,0), 'k-')
    plt.plot( (0,0), (-100, 100), 'k-')
    plt.plot( (1,3,4), (0,1,0), 'bo' )
    
    b0 = 1
    b1 = -.5
    d0 = -.125
    d1 = .25
    
    x = np.linspace(1,3,100)
    y = b0*(x-1) + d0*(x-1)**3
    plt.plot(x,y, 'r-')
    
    x = np.linspace(3,4,100)
    y = 1 + b1*(x-3) - .75*(x-3)**2 + d1*(x-3)**3
    plt.plot(x,y, 'g-')
    
    plt.xlim(0, 5)
    plt.ylim(0, 2)
    
    plt.show()
    
    
def not_a_knot_interp(xs, ys, us):
    #sort the data
    xs = np.array(xs)
    p = xs.argsort()
    xs = xs[p]
    ys = np.array(ys)[p]
    us = np.sort(us)
     
    n = len(xs)
    hs = xs[1:] - xs[:-1]
    
    delta = (ys[1:] - ys[:-1]) / hs
    
    A = np.zeros( (n,n) )
    A[0,0] = 1
    A[0,n-1] = -1
    for i in range(1, n-1):
        A[i, i-1] = hs[i]
        A[i, i] = 2*(hs[i]+hs[i-1])
        A[i, i+1] = hs[i-1]
    A[n-1, 0] = 2*hs[n-2]
    A[n-1, 1] = hs[n-2]
    A[n-1, n-2] = hs[0]
    A[n-1, n-1] = 2*hs[0]
    
    B = np.zeros( (n, 1) )
    for i in range(1, n-1):
        B[i] = 3*( hs[i]*delta[i-1] + hs[i-1]*delta[i] )
    B[n-1] = 3*( hs[-1]*delta[0] + hs[0]*delta[-1] )
    
    d = np.linalg.solve(A, B)
    d.shape = (n,)
    
    b = (d[:-1] -2*delta + d[1:])/(hs*hs)
    c = (3*delta - 2*d[:-1] - d[1:])/(hs)

    f = np.zeros( len(us) )
    j = 0
    for i, u in enumerate(us):
        while xs[j+1] < u:
            j += 1
        dis = u - xs[j] 
        f[i] = ys[j] + d[j]*dis + c[j]*dis**2 + b[j]*dis**3
    return f
    
def piecewise_interp(xs, ys, us):
    '''
    if type not in ['not-a-knot', 'periodic']:
        raise ValueError("Illegal argument. Type cannot equal " + str(type))
    '''
    #sort the data
    xs = np.array(xs)
    p = xs.argsort()
    xs = xs[p]
    ys = np.array(ys)[p]
    us = np.sort(us)
     
    n = len(xs)
    hs = xs[1:] - xs[:-1]
    
    delta = (ys[1:] - ys[:-1]) / hs
    
    A = np.zeros( (n,n) )
    for i in range(1, n-1):
        A[i, i-1] = hs[i]
        A[i, i] = 2*(hs[i]+hs[i-1])
        A[i, i+1] = hs[i-1]
    A[n-1, 0] = 2*hs[n-2]
    A[n-1, 1] = hs[n-2]
    A[n-1, n-2] = hs[0]
    A[n-1, n-1] = 2*hs[0]
    
    B = np.zeros( (n, 1) )
    for i in range(1, n-1):
        B[i] = 3*( hs[i]*delta[i-1] + hs[i-1]*delta[i] )
    
    '''
    if type == 'not-a-knot':
        A[0,0] = hs[1]**2
        A[0,1] = hs[1]**2 - hs[0]**2
        A[0,2] = -hs[0]**2
        B[0] = 2*(hs[1]**2 * delta[0] - hs[0]**2 * delta[1])
        A[n-1,n-3] = hs[-1]**2
        A[n-1,n-2] = hs[-1]**2 - hs[-2]**2
        A[n-1,n-1] = -hs[-2]**2
        B[n-1] = 2*(hs[-1]**2 * delta[-2] - hs[-2]**2 * delta[-1])
    '''
    #elif type == 'periodic':
    A[0,0] = 1
    A[0,n-1] = -1
    B[n-1] = 3*( hs[-1]*delta[0] + hs[0]*delta[-1] )
    
    d = np.linalg.solve(A, B)
    d.shape = (n,)
    
    b = (d[:-1] -2*delta + d[1:])/(hs*hs)
    c = (3*delta - 2*d[:-1] - d[1:])/(hs)

    f = np.zeros( len(us) )
    j = 0
    for i, u in enumerate(us):
        while xs[j+1] < u:
            j += 1
            if j >= n-1:
                j = n-2
                
                break
        dis = u - xs[j] 
        f[i] = ys[j] + d[j]*dis + c[j]*dis**2 + b[j]*dis**3
    return f
    
def p2c_helper(n=5):
    xs = np.linspace(0,1, n+1)
    ys = np.random.rand(n+1)
    ys[-1] = ys[0]
    
    us = np.linspace(0,1,1000)
    fs = piecewise_interp(xs, ys, us)
    
    plt.plot( (-100, 100), (0,0), 'k-')
    plt.plot( (0,0), (-100, 100), 'k-')
    plt.plot( xs, ys, 'bo')
    plt.plot( us, fs, 'r-')
    
    plt.xlim( (-0.1, 1.1) )
    plt.ylim( (-0.1, 1.1) )
    plt.show()
    
def p2d():
    xs = []
    for i in range(0,6):
        xs.append( np.sin(np.pi * i / 10) - 1 )
    for i in range(6, 11):
        xs.append( 1 - np.sin(np.pi * i / 10) )
    xs = np.array(xs)
    print(xs)
    ys = np.sin( np.pi*(xs - .01) )
    us = np.linspace(-1, 1, 100)
    fs = piecewise_interp(xs, ys, us)
    
    plt.plot( (-100, 100), (0,0), 'k-')
    plt.plot( (0,0), (-100, 100), 'k-')
    plt.plot( xs, ys, 'bo')
    plt.plot( us, fs, 'r-')
    
    #fs = piecewise_interp(xs, ys, us, type='not-a-knot')
    #plt.plot( us, fs, 'g-')
    
    cs = CubicSpline(xs,ys)
    plt.plot(us, cs(us), 'g-')
    
    plt.xlim( (-1.1, 1.1) )
    plt.ylim( (-1.1, 1.1) )
    plt.show()
    
    