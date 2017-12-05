import numpy as np
import matplotlib.pyplot as plt
from pytex import *

def trap_int(f, a, b, n):
    h = (b-a)/n
    x = a + h
    T = 0
    for i in range(n-1):
        T += f(x)
        x += h
    T += (f(a)+f(b))/2
    T *= h
    return T
    
def f1(x):
    return (x-1)**2 * np.exp(-x**2)
    
def f2(x):
    return np.exp( np.cos( np.pi*x ) )
    
int1_actual = 1.87259295726583875
int2_actual = 2.53213175550401667
    
def p1a():
    ns = [4, 8, 16, 32]
    T1 = []
    T2 = []
    for n in ns:
        T1.append( abs( trap_int(f1, -1, 1, n) - int1_actual ) )
        T2.append( abs( trap_int(f2, -1, 1, n) - int2_actual ) )
    latex_table( (ns, T1, T2), ('$n$', '$I(f_1)$ error', '$I(f_2)$ error') ) 
    
def df1(x):
    return 2*(x-1)*np.exp(-x**2) - 2*x*(x-1)**2 * np.exp(-x**2)
    
def df2(x):
    return -np.pi*np.sin(np.pi*x)*np.exp(np.cos(np.pi*x))
    
def trap_int_corrected(f, df, a, b, n):
    T = trap_int(f,a,b,n)
    h = (b-a)/n 
    return T - h**2/12*(df(b)-df(a))
  
def p1b():
    int1_actual = 1.87259295726583875
    int2_actual = 2.53213175550401667
    ns = [4, 8, 16, 32]
    T1 = []
    T2 = []
    for n in ns:
        T1.append( abs( trap_int_corrected(f1, df1, -1, 1, n) - int1_actual ) )
        T2.append( abs( trap_int_corrected(f2, df2, -1, 1, n) - int2_actual ) )
    latex_table( (ns, T1, T2), ('$n$', '$I(f_1)$ error', '$I(f_2)$ error') ) 
    
    
def p1c():
    ns = np.arange(40)+20
    #ns = np.array([4, 8, 16, 32])
    T_err = []
    for n in ns:
        T_err.append( abs( trap_int(f1, -1, 1, n) - int1_actual ) )
    plt.plot(2/ns**2, T_err, 'b.')
    plt.xlabel('h^2')
    plt.ylabel('error')
    plt.show()
    
    T_err = []
    for n in ns:
        T_err.append( abs( trap_int_corrected(f1, df1, -1, 1, n) - int1_actual ) )
    plt.plot(2/ns**4, T_err, 'b.')
    plt.xlabel('h^4')
    plt.ylabel('error')
    plt.show()
    
    
    
    
    
    
    
    
    
    
    
    