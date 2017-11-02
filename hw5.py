import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import CubicSpline
from hw4 import polyinterp

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
    ys = np.sin( np.pi*(xs - 0.1) )
    us = np.linspace(-1, 1, 100)
    fs = piecewise_interp(xs, ys, us)
    
    plt.plot( (-100, 100), (0,0), 'k-')
    plt.plot( (0,0), (-100, 100), 'k-')
    plt.plot( xs, ys, 'bo')
    plt.plot( us, fs, 'r-')
    
    cs = CubicSpline(xs,ys)
    plt.plot(us, cs(us), 'g-')
    
    plt.xlim( (-1.1, 1.1) )
    plt.ylim( (-1.1, 1.1) )
    plt.show()
    
def p3():
    str = open('res/hw5_DotToDot.txt', 'r').read()
    x = []
    y = []
    for l in str.split('\n'):
        xy = l.split('\t')
        x.append( float(xy[0]))
        y.append( float(xy[1]))
    t = np.linspace(0, 1, len(x))
    u = np.linspace(0, 1, len(x)*100)
    fx = piecewise_interp(t, x, u)
    fy = piecewise_interp(t, y, u)
    
    plt.plot( (-100, 100), (0,0), 'k-')
    plt.plot( (0,0), (-100, 100), 'k-')
    plt.plot( x, y, 'bo')
    plt.plot( x, y, 'g-')
    plt.plot( fx, fy, 'r-')
    
    plt.xlim( (0, 1.1) )
    plt.ylim( (0, 1.1) )
    plt.show()
    
def trig_interp(y, u):
    N = len(y)
    if N%2 == 0:
        raise ValueError("y must contain an odd number of values.")
    k = np.linspace(0, 2*np.pi, N, endpoint=False)
    f = []
    sign = np.array( [ (-1)**k for k in range(0,N)] )
    for t in u:
        if t in k:
            p = y[ int(np.argwhere(k==t)) ]
        else:
            denominator = sign / np.sin( (t - k)/2 )
            numerator = denominator * y
            p = np.sum(numerator)/np.sum(denominator)
        f.append(p)
        
    return f

def p6b_helper(t):
    return np.sin(2*t)*np.cos(3*t)
    
def p6b():
    k = np.linspace(0, 2*np.pi, 11, endpoint=False)
    y = p6b_helper(k)
    
    x = np.linspace(0,2*np.pi, 111)
    f = trig_interp(y, x)
    
    plt.plot( (-100, 100), (0,0), 'k-')
    plt.plot( (0,0), (-100, 100), 'k-')
    plt.plot( x, f, 'r-')
    plt.plot( k, y, 'bo')
    plt.plot( x, p6b_helper(x), 'g-')
    plt.xlim( (-.1, 2*np.pi+.1) )
    plt.ylim( (-1.1, 1.1) )
    plt.show()
    
    plt.plot( (-100, 100), (0,0), 'k-')
    plt.plot( (0,0), (-100, 100), 'k-')
    plt.plot( x, f - p6b_helper(x), 'r-')
    plt.xlim( (-.1, 2*np.pi+.1) )
    plt.ylim( (-2e-15, 2e-15) )
    plt.show()
    
def p6c_helper(t):
    return np.sin(20*np.exp(np.cos(t - .5)))
    
def p6c():
    k = np.linspace(0, 2*np.pi, 101, endpoint=False)
    y = p6c_helper(k)
    
    x = np.linspace(0,2*np.pi, 1001)
    f = trig_interp(y, x)
    
    chebyshev = np.array( [ np.pi*(1 - np.cos(j*np.pi/100)) for j in range(101)] )
    y2 = p6c_helper(chebyshev)
    f2 = polyinterp(x, chebyshev, y2)
    
    plt.plot( (-100, 100), (0,0), 'k-')
    plt.plot( (0,0), (-100, 100), 'k-')
    plt.plot( x, f, 'r-')
    plt.plot( x, f2, 'g-')
    plt.plot( k, y, 'bo')
    plt.plot( x, p6c_helper(x), 'g-')
    plt.xlim( (-.1, 2*np.pi+.1) )
    plt.ylim( (-1.1, 1.1) )
    plt.show()
    
    plt.plot( (-100, 100), (0,0), 'k-')
    plt.plot( (0,0), (-100, 100), 'k-')
    plt.plot( x, f - p6c_helper(x), 'r-')
    plt.xlim( (-.1, 2*np.pi+.1) )
    plt.ylim( (-3e-6, 3e-6) )
    plt.show()
    
    plt.plot( (-100, 100), (0,0), 'k-')
    plt.plot( (0,0), (-100, 100), 'k-')
    plt.plot( x, f2 - p6c_helper(x), 'r-')
    print(np.max(f2-p6c_helper(x)))
    plt.xlim( (-.1, 2*np.pi+.1) )
    plt.ylim( (-1e-2, 1e-2) )
    plt.show()
    
def square(t):
    if t%(2*np.pi) < np.pi:
        return 1
    else:
        return -1
    
def p6d():
    np_square = np.vectorize(square)
    k = np.linspace(0, 2*np.pi, 101, endpoint=False)
    y = np_square(k)
    x = np.linspace(0, 6*np.pi, 10001)
    f = trig_interp(y,x)
    
    plt.plot(x, np_square(x), 'b-')
    plt.plot(x, f, 'r-')
    
    plt.xlim(0,6*np.pi)
    plt.show()
    
    print("Max amplitude: %f" % np.max(f))
    
