import numpy as np
from pytex import *
from decimal import *
from matplotlib import pyplot as plt





# Problem 1     ********************************************************************

def secant(f, x0, x1, epsilon=1e-16, n_max=10**5):
    xs = [x0,x1]
    xn = x1
    xn_1 = x0
    fn = f(xn)
    fn_1 = f(xn_1)
    for i in range(n_max):
        if abs(fn) < epsilon:
            return xs
        xn, xn_1 = xn - (xn-xn_1)/(fn-fn_1)*fn, xn
        fn, fn_1 = f(xn), f(xn_1)
        xs.append(xn)
    return None

def newton(f, df, x0, epsilon=1e-16, n_max=10**5):
    x = x0
    xs =[x]
    fx = f(x)
    dfx = df(x)
    for i in range(n_max):
        if abs(fx) < epsilon:
            return xs
        x = x - fx/dfx
        xs.append(x)
        fx, dfx = f(x), df(x)
    return None

def foo_p1(x):
    return np.exp(-1*x)-x

def dfoo_p1(x):
    return -1*np.exp(-1*x)-1

def hw2_p1c():
    root = 0.56714329040978
    sec = secant(foo_p1, 5.5, 5.25, epsilon=1e-14)
    err_rel = []
    err_ratio = ['N/A']
    ns = []
    for i, x in enumerate(sec):
        err_rel.append(abs(x - root)/root)
        if i != 0:
            err_ratio.append(np.log(err_rel[i])/np.log(err_rel[i-1]))
        ns.append(i+1)
    latex_table((ns, sec, err_rel, err_ratio), ("$n$", "$x_n$", "relative err", 'ratio'))
    
    newt = newton(foo_p1, dfoo_p1, 5.25, epsilon=1e-14)
    err_rel = []
    err_ratio = ['N/A']
    ns = []
    for i, x in enumerate(newt):
        err_rel.append(abs(x - root)/root)
        if i != 0:
            err_ratio.append(np.log(err_rel[i])/np.log(err_rel[i-1]))
        ns.append(i+1)
    latex_table((ns, newt, err_rel, err_ratio), ("$n$", "$x_n$", "relative err", "ratio"))
    
# Problem 3     ********************************************************************
def my_sqrt(a, x0, epsilon=1e-15, n_max=10**5):
    xs = [x0]
    x = x0
    for i in range(n_max):
        if abs(x**2 - a) < epsilon:
            return xs
        x = x/2 * (3 - x**2/a)
        xs.append(x)
    return  None
    

def hw2_p3a():
    xs = my_sqrt(2, 1)
    root = np.sqrt(2)
    err_rel = []
    err_ratio = ['N/A']
    ns = []
    for i, x in enumerate(xs):
        err_rel.append(abs(x - root))
        if i != 0:
            err_ratio.append(np.log(err_rel[i])/np.log(err_rel[i-1]))
        ns.append(i+1)
    latex_table((ns, xs, err_rel, err_ratio), ("$n$", "$x_n$", "abs err", 'ratio'))
    
def foo_3(x,a):
    return .5*(x + a/x)
def foo_4(x,a):
    return x/2 * (3-x**2/a)
      
def hw2_p3b():
    xs = np.linspace(-.1,10,1000)
    ys = np.vectorize(foo_3)(xs,2)
    plt.plot(xs,ys, 'b-')
    plt.plot(xs,xs, 'r-')
    plt.xlim(0, 10)
    plt.ylim(0,10)
    #plt.show()
    
    xs = np.linspace(-2,5,1000)
    ys = np.vectorize(foo_4)(xs,2)
    plt.plot(xs,ys, 'b-')
    plt.plot(xs,xs, 'r-')
    plt.xlim(-2, 5)
    plt.ylim(-2,5)
    #plt.show()
    
    #check for convergence
    starts = [.0001, 1, 2, 4, 100]
    ends = []
    epsilon = 1e-13
    for x0 in starts:
        x = x0
        for i in range(10**6):
            x_prev = x
            x = foo_3(x,2)
            if abs(x-x_prev) < epsilon:
                break
        ends.append(x)
    latex_table((starts, ends), ('$x_0$', 'converged'))
    
    starts = [-100, -10, -5, -1*np.sqrt(10/3)-.001, -1*np.sqrt(10/3)+.001,-1, -.01, .01,1, 1*np.sqrt(10/3)-.001, 1*np.sqrt(10/3)+.001, 1.999,2, 2.001 ,3,3.1, 3.15, 3.2, 3.4, 3.5,3.6,3.7,3.84,5, 10, 100]
    ends = []
    for x0 in starts:
        x = x0
        for i in range(10**2):
            x_prev = x
            x = foo_4(x,2)
            if abs(x-x_prev) < epsilon:
                break
            if abs(x) > 5000:
                break
        ends.append(x)
    latex_table((starts, ends), ('$x_0$', 'converged'))
    
def my_inv_root(a, x0, epsilon=1e-15, n_max=10**5):
    xs = [x0]
    x = x0
    for i in range(n_max):
        if abs(x**2 - a) < epsilon:
            return xs
        x = x/2 * (3 - x**2*a)
        xs.append(x)
    return  None
def hw2_p3c():
    xs = my_inv_root(2, .07)
    root = 1/np.sqrt(2)
    err_rel = []
    err_ratio = ['N/A']
    ns = []
    for i, x in enumerate(xs):
        err_rel.append(abs(x - root))
        if i != 0:
            err_ratio.append(np.log(err_rel[i])/np.log(err_rel[i-1]))
        ns.append(i+1)
    latex_table((ns, xs, err_rel, err_ratio), ("$n$", "$x_n$", "abs err", 'ratio'))
    
    
