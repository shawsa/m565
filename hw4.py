import numpy as np
import matplotlib.pyplot as plt
from matplotlib2tikz import save as tikz_save


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
        for xi in x[np.arange(len(x))!=j]: 
            w[j] /= (xj - xi)
    return w

def p4_foo(x):
    return np.sqrt( x**2 + .1)

def p4b(tikz=False):
    samples = []
    labels = []
    x = np.linspace(-1,1,15)
    samples.append(x)
    labels.append('Equispaced points:')
    x = []
    for j in range(0,15):
        x.append( -1* np.cos( (2*j+1)/29 * np.pi))
    samples.append( np.sort( np.array(x) ) )
    labels.append('First kind Chebyshev points:')
    x = []
    for j in range(0,15):
        x.append( -1* np.cos( (j*np.pi)/14 ) )
    samples.append( np.sort( np.array(x) ) )
    labels.append('Second kind Chebyshev points:')
    x = [-0.987992518020485,
            -0.394151347077563,
            0.570972172608539,
            -0.937273392400706,
            -0.201194093997435,
            0.724417731360170,
            -0.848206583410427,
            0,
            0.848206583410427,
            -0.724417731360170,
            0.201194093997435,
            0.937273392400706,
            -0.570972172608539,
            0.394151347077563,
            0.987992518020485]
    samples.append( np.sort( np.array(x) ) )
    labels.append('Legendre points:')
    for x, l, file in zip(samples, labels, range(1,5)):
        y = p4_foo(x)
        u = np.linspace(-1,1,201)
        p = polyinterp(u, x, y)
        f = p4_foo(u)
        print(l)
        print('Eucidian norm: %f' % np.linalg.norm(f-p))
        print('Infinity norm: %f' % np.linalg.norm(f-p,np.inf))
        plt.subplot(2,1,1)
        plt.plot(u,p, 'r-')
        plt.plot(u,f, 'k-')
        plt.plot(x,y, 'bo')
        plt.xlim( (-1, 1) )
        plt.ylim( (.2, 1.2) )
        
        plt.subplot(2,1,2)
        plt.plot(u, p-f, 'k-')
        plt.xlim( (-1, 1) )
        y_range = np.max( np.abs(p-f) )
        plt.ylim( (-y_range, y_range) )
        if tikz:
            tikz_save('images/hw4_figure_' + str(file) + '_tikz.tex')
        plt.show()

def p5a():
    x = np.array( [0,2,4] )
    y = np.array( [0,1,2] )
    
    u = np.linspace(0,4, 200)
    f = polyinterp(u, x, y)
    
    plt.plot(u,f, 'r-')
    
    plt.plot( (0,2,4), (0,1,2), 'bo')
    plt.plot( (-10,10), (2,2), 'k-')
    plt.plot( (-10,10), (0,0), 'k-')
    
    plt.xlim( (-0.5, 4.5) )
    plt.ylim( (-.5, 2.5) )
    
    plt.show()
    
def p5b():
    
    x = np.array( [0,2,4] )
    y = np.array( [0,1,2] )
    
    u = np.linspace(0,4, 200)
    f = .25 * u**2 - 1/16 * u**2 * (u-2)
    
    plt.plot(u,f, 'r-')
    
    plt.plot( (0,2,4), (0,1,2), 'bo')
    plt.plot( (-10,10), (2,2), 'k-')
    plt.plot( (-10,10), (0,0), 'k-')
    
    plt.xlim( (-0.5, 4.5) )
    plt.ylim( (-.5, 2.5) )
    
    plt.show()
    
    