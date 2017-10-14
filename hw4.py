import numpy as np
import matplotlib.pyplot as plt


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

def p4b():
    samples = []
    x = np.linspace(-1,1,15)
    samples.append(x)
    for x in samples:
        y = p4_foo(x)
        u = np.linspace(-1,1,201)
        p = polyinterp(u, x, y)
        f = p4_foo(u)
        print('Eucidian norm: %f' % np.linalg.norm(f-p))
        print('Infinity norm: %f' % np.linalg.norm(f-p,np.inf)) 
        plt.plot(u,p, 'r-')
        plt.plot(u,f, 'k-')
        plt.plot(x,y, 'bo')
        plt.show()
