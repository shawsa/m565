import numpy as np
from pytex import *
from decimal import *
from matplotlib import pyplot as plt
from math import erf, pi, sqrt, exp, ceil

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

def fixed_point(x0, foo, a, n):
    xs  = [x0]
    ys = [-100]
    for i in range(n-1):
        if abs(xs[-1]) > 1000:
            break
        xs.append(xs[-1])
        ys.append(foo(xs[-1],a))
        xs.append(ys[-1])
        ys.append(xs[-1])
    return xs, ys
      
def hw2_p3b():
    xs = np.linspace(-.1,10,1000)
    ys = np.vectorize(foo_3)(xs,2)
    plt.plot(xs,ys, 'b-')
    plt.plot(xs,xs, 'r-')
    fp_xs, fp_ys = fixed_point(.2, foo_3, 2, 10)
    plt.plot(fp_xs, fp_ys, 'g-')
    fp_xs, fp_ys = fixed_point(8.5, foo_3, 2, 10)
    plt.plot(fp_xs, fp_ys, 'y-')
    plt.xlim(0, 10)
    plt.ylim(0,10)
    plt.show()
    
    xs = np.linspace(-10,10,1000)
    ys = np.vectorize(foo_4)(xs,2)
    plt.plot(xs,ys, 'b-')
    plt.plot(xs,xs, 'r-')
    fp_xs, fp_ys = fixed_point(.2, foo_4, 2, 10)
    plt.plot(fp_xs, fp_ys, 'g-')
    fp_xs, fp_ys = fixed_point(3, foo_4, 2, 10)
    plt.plot(fp_xs, fp_ys, 'c-')
    fp_xs, fp_ys = fixed_point(2, foo_4, 2, 10)
    plt.plot(fp_xs, fp_ys, 'm-')
    fp_xs, fp_ys = fixed_point(3.1, foo_4, 2, 10)
    plt.plot(fp_xs, fp_ys, 'y-')
    fp_xs, fp_ys = fixed_point(3.2, foo_4, 2, 10)
    plt.plot(fp_xs, fp_ys, 'k-')
    plt.xlim(-5, 5)
    plt.ylim(-5,5)
    plt.show()
    
    #check for convergence
    xs = []
    converges = []
    for x0 in np.arange(-4,4,.1):
        fp_xs, fp_ys = fixed_point(x0, foo_4, 2, 10000)
        xs.append(x0)
        if abs(fp_xs[-1] - fp_xs[-2]) < 1e-10:
            converges.append(fp_xs[-1])
        else:
            converges.append('diverges')
    latex_table((xs, converges), ('$x0$', 'converges to'))
    
    
    
    
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
    
    starts = [0, .01,1, np.sqrt(2), 2, np.sqrt(6)-.001, 1*np.sqrt(6)+.001]
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
    
    
# Problem 5     ********************************************************************

def foo_6(x):
    return -1*erf(2*(x-1))
def dfoo_6(x):
    return -4/sqrt(pi)*exp(-4*(x-1)**2)

def hw2_p5a():
    xs = [0]
    x = 0
    for i in range(30):
        x = foo_6(x)
        xs.append(x)
    latex_table((range(1,31), xs),('$n$', '$x_n$'))
    
def hw2_p5b():
    getcontext().prec = 30
    xs = [0]
    x = 0
    ns = [0]
    #perform first iteration so that error checking can compute
    g = foo_6(x)
    x = x - (g-x)**2/(foo_6(g) -2*g + x)
    xs.append(x)
    for i in range(1,1000000): 
        if abs(xs[-1] - xs[-2]) < 1e-15:
            break
        g = foo_6(x)
        if abs((foo_6(g) -2*g + x)) < 1e-15:
            break
        print(x)
        x = x - (g-x)**2/(foo_6(g) -2*g + x)
        xs.append(x)
        ns.append(i)
    ns.append(ns[-1]+1)
    xs.append(foo_6(x))
        
    latex_table((ns, xs),('$n$', '$x_n$'))
    #part c
    print(dfoo_6(xs[-1]))

    
# Problem 7     ********************************************************************
 
def foo_7(x):
    return np.exp(x) - 1
    
def fixed_point_7(x0, foo, n):
    xs  = [x0]
    ys = [-100]
    for i in range(n-1):
        if abs(xs[-1]) > 1000:
            break
        xs.append(xs[-1])
        ys.append(foo(xs[-1]))
        xs.append(ys[-1])
        ys.append(xs[-1])
    return xs, ys
    
def hw2_p7a():
    xs = np.linspace(-10,10,1000)
    ys = foo_7(xs)
    plt.plot(xs,ys, 'b-')
    plt.plot(xs,xs, 'r-')
    fp_xs, fp_ys = fixed_point_7(.2, foo_7, 10)
    plt.plot(fp_xs, fp_ys, 'g-')
    fp_xs, fp_ys = fixed_point_7(.1, foo_7, 10)
    plt.plot(fp_xs, fp_ys, 'y-')
    fp_xs, fp_ys = fixed_point_7(-.1, foo_7, 10)
    plt.plot(fp_xs, fp_ys, 'c-')
    fp_xs, fp_ys = fixed_point_7(-2, foo_7, 10)
    plt.plot(fp_xs, fp_ys, 'm-')
    fp_xs, fp_ys = fixed_point_7(-5, foo_7, 10)
    plt.plot(fp_xs, fp_ys, 'k-')
    plt.xlim(-5, 5)
    plt.ylim(-5,5)
    plt.show()
    
def pred_n(epsilon, x0):
    return 2/epsilon - 2/x0
    
def hw2_p7b():
    eps = [1e-4, 1e-5, 1e-6]*3
    x0s = [-.1]*3 + [-.01]*3 + [-.001]*3
    pred = []
    actual = []
    for e, x0 in zip(eps, x0s):
        pred.append( ceil( pred_n(e,x0) ) )
        x = x0
        for i in range(pred[-1]*2):
            if x > -1*e:
                actual.append(i)
                break
            x = foo_7(x)
        else:
            actual.append('>2*predicted')
    rel_errors = []
    for p, a in zip(pred, actual):
        if a == '>2*predicted':
            rel_errors.append(' ')
        else:
            rel_errors.append('{:0.2f}\\%'.format(100*(p-a)/a))
    latex_table((eps, x0s, pred, actual, rel_errors), ('$\epsilon$', '$x_0$', '$n$ predicted', '$n$ actual', 'relative error'))
        
    
    
