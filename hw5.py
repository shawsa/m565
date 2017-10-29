import numpy as np
from matplotlib import pyplot as plt

def p1():
    plt.plot( (-100, 100), (0,0), 'k-')
    plt.plot( (0,0), (-100, 100), 'k-')
    plt.plot( (1,3,4), (0,1,0), 'bo' )
    
    b0 = .5
    b1 = 0
    d0 = 0
    d1 = -.25
    
    x = np.linspace(1,3,100)
    y = b0*(x-1) + d0*(x-1)**3
    plt.plot(x,y, 'r-')
    
    x = np.linspace(3,4,100)
    y = 1 + b1*(x-3) - .75*(x-3)**2 + d1*(x-3)**3
    plt.plot(x,y, 'g-')
    
    plt.xlim(0, 5)
    plt.ylim(0, 2)
    
    plt.show()