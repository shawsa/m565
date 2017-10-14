import numpy as np
import matplotlib.pyplot as plt

def polyinterp(u, x, y, w=None):
    n = len(x)
    
    return f
    
def baryweights(x):
    w = np.ones(len(n))
    for j, xj in enumerate(x):
        for xi in x[:j]+[j+1:]
            w[j] /= (xj - xi)
    return w