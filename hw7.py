import numpy as np
import matplotlib.pyplot as plt
from pytex import *
from sympy import *
#from ghalton import *

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
    
    
def p3b():
    xs = [0.02913447215197205,
        0.1739772133208976,
        0.4117025202849020,
        0.6773141745828204,
        0.8947713610310083]
    ws = [0.2978934717828945,
        0.3497762265132242,
        0.2344882900440524,
        0.09893045951663315,
        0.01891155214319580]
    correct = 0.239811742000565
    
    ys = np.sin(xs)
    
    GQ = np.dot(ws,ys)
    print("Gaussian Quadrature:")
    print(GQ)
    print("Error: %f")
    print((GQ - correct))
    
def p3c():
    xs = [0.04691007703066800,
        0.2307653449471585,
        0.5000000000000000,
        0.7692346550528415,
        0.9530899229693320]
    ws = [0.1184634425280945,
        0.2393143352496832,
        0.2844444444444444,
        0.2393143352496832,
        0.1184634425280945]
    correct = 0.239811742000565
    
    ys = np.sin(xs)*np.log(1.0/np.array(xs))
    
    GQ = np.dot(ws,ys)
    print("Gaussian Quadrature:")
    print(GQ)
    print("Error: %f")
    print((GQ - correct))
    
def p3d():
    xs = [0.01088567092697150,
        0.05646870011595235,
        0.1349239972129753,
        0.2404519353965941,
        0.3652284220238275,
        0.5000000000000000,
        0.6347715779761725,
        0.7595480646034059,
        0.8650760027870247,
        0.9435312998840476,
        0.9891143290730285]
    ws = [0.02783428355808683,
        0.06279018473245231,
        0.09314510546386713,
        0.1165968822959952,
        0.1314022722551233,
        0.1364625433889503,
        0.1314022722551233,
        0.1165968822959952,
        0.09314510546386713,
        0.06279018473245231,
        0.02783428355808683]
    correct = 0.239811742000565
    ys = np.sin(xs)*np.log(1.0/np.array(xs))
    GQ = np.dot(ws,ys)
    print("Gaussian Quadrature:")
    print(GQ)
    print("Error: %f")
    print((GQ - correct))
    
def p3d_print():
    xs = [0.01088567092697150,
        0.05646870011595235,
        0.1349239972129753,
        0.2404519353965941,
        0.3652284220238275,
        0.5000000000000000,
        0.6347715779761725,
        0.7595480646034059,
        0.8650760027870247,
        0.9435312998840476,
        0.9891143290730285]
    ws = [0.02783428355808683,
        0.06279018473245231,
        0.09314510546386713,
        0.1165968822959952,
        0.1314022722551233,
        0.1364625433889503,
        0.1314022722551233,
        0.1165968822959952,
        0.09314510546386713,
        0.06279018473245231,
        0.02783428355808683]
    i = np.arange(11)
    latex_table( (i, xs, ws), ('$i$', '$x_i$', '$w_i$') )
    
def p4f(x):
    return 4/(1+x**2)
 
def p4a():
    N=5
    xs = np.linspace(0,1,2**N+1,endpoint=True)
    fs = p4f(xs)
    fs[0] *= .5
    fs[-1] *= .5
    
    Rom = np.zeros((N+1,N+1))
    for i in range(0,N+1):
        Rom[i][0] = np.sum(fs[::2**(N-i)])/2**i
        
    for m in range(1, N+1):
        for j in range(1, m+1): #j = k-1
            Rom[m,j] = (4**j * Rom[m,j-1] - Rom[m-1,j-1]) / (4**j - 1)
    ns = 2**np.arange(0,N+1)
    Rs = np.array([Rom[i][i] for i in range(0,N+1)])
    err = Rs - np.pi
    latex_table( (ns, Rs, err), ('$n$', '$R_{i,i}$', 'err') )
    
    
def p5a():
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')
    w = Symbol('w')
    a = Symbol('a')
    b = Symbol('b')
    return integrate(sin(pi/2*(x+y+z+w+a+b)), (x,0,1),(y,0,1),(z,0,1),(w,0,1),(a,0,1),(b,0,1))

def p5b():
    N = 10000
    points = np.random.rand(N,6)
    col_sum = np.sum(points, axis=1)
    QMC = sum( np.sin(np.pi/2 * col_sum) ) / N
    return QMC

'''
def p5c():
    N = 10000
    seq = GeneralizedHalton(6)
    seq.reset()
    points = seq.get(N)
    col_sum = np.sum(points, axis=1)
    QMC = sum( np.sin(np.pi/2 * col_sum) ) / N
    return QMC
'''


    
gamma = 0.5772156649015328
    
def p6_f(x):
    return 1/x + np.log(1-1/x)
def p6_f1(x):
    return -x**-2 + (x-1)**-1 - x**-1
def p6_f3(x):
    return -6*x**-4 + 2*(x-1)**-3 - 2*x**-3
def p6_f5(x):
    return -120*x**-6 + 24*(x-1)**-5 - 24*x**-5
def euler_maclaurin_terms(N):
    return p6_f(N)/2 - p6_f1(N)/12 + p6_f3(N)/720 - p6_f5(N)/30240
    
def gamma_approx():
    gamma_approx = 1
    for k in range(2,20):
        gamma_approx += 1/k + np.log(1-1/k)
    gamma_approx += 38*np.arctanh(1/39)-1 + euler_maclaurin_terms(20)
    return gamma_approx
    
def p6b():
    gam = gamma_approx()
    print(gam)
    print(gam - gamma)
    

def p6c():
    approx = 1
    for k in range(2,1003):
        approx += 1/k + np.log(1-1/k)
    print(approx - gamma)
    

    
    
