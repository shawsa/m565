import numpy as np
from scipy.linalg import *

def foo_3b_bad(x,v,A,B,C,T):
    k2 = inv(B.transpose()).dot(x)
    k2 = x.transpose().dot(k2)
    
    k1 = inv(A).dot(x)
    k1 = inv(T).dot(k1)
    k1 = C.dot(k1)
    k1 = inv(B).dot(k1)
    k1 = v.transpose().dot(k1)
    
    return k1 + k2

def foo_3b(x,v,A,B,C,T):
    n = A.shape[0]
    
    P, L, U = lu(B)
    
    k2 = solve_triangular(U.transpose(), x, lower=True)
    k2 = solve_triangular(L.transpose(), k2, lower=False)
    k2 = P.dot(k2)
    k2 = x.transpose().dot(k2)
    
    H = cholesky(A, lower=True)
    
    k1 = solve_triangular(H,x, lower=True)
    k1 = solve_triangular(H.transpose(),k1, lower=False)
    
    #convert T to a banded matrix datatype
    bandedT = np.zeros( (3,T.shape[0]) )
    bandedT[0][1:n] = T.diagonal(1)
    bandedT[1][0:n] = T.diagonal(0)
    bandedT[2][0:n-1] = T.diagonal(-1)
    
    k1 = solve_banded( (1,1), bandedT, k1 )
    
    k1 = C.dot(k1)
    k1 = P.transpose().dot(k1)
    k1 = solve_triangular(L, k1, lower=True)
    k1 = solve_triangular(U, k1, lower=False)
    
    k1 = v.transpose().dot(k1)
    
    return k1 + k2
    
def foo_3_helper(n):
    A = np.random.rand(n,n)
    A = A.dot(A.transpose()) + np.identity(n)
    B = np.random.rand(n,n)
    C = np.random.rand(n,n)
    T = np.diag(np.random.rand(n-1), -1) 
    T += np.diag(np.random.rand(n), 0) 
    T += np.diag(np.random.rand(n-1), 1)
    x = np.random.rand(n,1)
    v = np.random.rand(n,1)
    return x,v,A,B,C,T
    
def p3c():
    for i in range(20):
        x,v,A,B,C,T = foo_3_helper(100)
        print(foo_3b(x,v,A,B,C,T))
        
def p4c_helper(n):
    P = np.zeros( (n,n) )
    P += np.diag(np.random.rand(n))
    P += np.diag(np.random.rand(n-1), 1)
    P += np.diag(np.random.rand(n-1), -1)
    P += np.diag(np.random.rand(n-2), 2)
    P += np.diag(np.random.rand(n-2), -2)

    return P, np.random.rand(n)
    
def penta_solve(P,v):
    n = P.shape[0]
    a = P.diagonal(0)
    b = np.zeros(n)
    b[1:n] = P.diagonal(-1)
    c = P.diagonal(1)
    d = np.zeros(n)
    d[2:n] = P.diagonal(-2)
    e = P.diagonal(2)
    l = np.zeros(n)
    u = np.zeros(n-1)
    w = np.zeros(n-2)
    k = np.zeros(n)
    
    #step 1
    l[0] = a[0]
    u[0] = c[0]/l[0]
    w[0] = e[0]/l[0]
    
    #step 2
    k[1] = b[1]
    l[1] = a[1]-k[1]*u[0]
    u[1] = (c[1]-b[1]*w[0])/l[1]
    w[1] = e[1]/l[1]
    
    #step j
    for j in range(2,n-2):
        k[j] = b[j] - d[j]*u[j-2]
        l[j] = a[j] - d[j]*w[j-2] - k[j]*u[j-1]
        u[j] = (c[j]-k[j]*w[j-1])/l[j]
        w[j] = e[j]/l[j]
     
    #step n-1
    k[n-2] = b[n-2] - d[n-2]*u[n-2-2]
    l[n-2] = a[n-2] - d[n-2]*w[n-2-2] - k[n-2]*u[n-2-1]
    u[n-2] = (c[n-2]-k[n-2]*w[n-2-1])/l[n-2]
    
    #step n
    k[n-1] = b[n-1] - d[n-1]*u[n-1-2]
    l[n-1] = a[n-1] - d[n-1]*w[n-1-2] - k[n-1]*u[n-1-1]
    
    '''
    L = np.zeros( (n,n) )
    U = np.zeros( (n,n) )
    U += np.diag([1]*n,0)
    U += np.diag(u, 1)
    U += np.diag(w, 2)
    L += np.diag(l,0)
    L += np.diag(d[2:n], -2)
    L += np.diag(k[1:n], -1)
    return L, U
    '''
    
    #forward sub L
    y = np.zeros(n)
    #step 1
    y[0] = v[0]/l[0]
    #step 2
    y[1] = (v[1]-k[1]*y[0])/l[1]
    #step j
    for j in range(2,n):
        y[j] = (v[j] - d[j]*y[j-2] - k[j]*y[j-1])/l[j]
        
    #back sub U
    x = np.zeros(n)
    #step n
    x[n-1] = y[n-1]
    #step n-1
    x[n-2] = y[n-2] - u[n-2]*x[n-1]
    #step j
    for j in range(n-3, -1, -1):
        x[j] = y[j] - u[j]*x[j+1] - w[j]*x[j+2]
    x.shape = (n,1)
    return x
    
def p4d_helper(n):
    a = [i for i in range(1,n+1)]
    b = [-(i+1)/3 for i in range(1, n)]
    d = [-(i+2)/6 for i in range(1,n-1)]
    P = np.zeros( (n,n) )
    P += np.diag(a, 0)
    P += np.diag(b, 1)
    P += np.diag(b, -1)
    P += np.diag(d, 2)
    P += np.diag(d, -2)
    f = np.zeros( (n,1) )
    f[0] = .5
    f[1] = 1/6.0
    f[n-2] = 1/6.0
    f[n-1] = .5
    return P, f
    
def p4d():
    for n in [100, 1000]:
        P,f = p4d_helper(n)
        x = penta_solve(P,f)
        r = P.dot(x)-f
        print('residual %f' % r.transpose().dot(r))
        dif = x - solve(P,f)
        print('GE %f' % dif.transpose().dot(dif))
    
    

def p5b():
	A = np.zeros((5,5))
	A[0] = [1/2, 1/3, 1/4, 1/5, 1/6]
	A[1] = [1/3, 1/4, 1/5, 1/6, 1/7]
	A[2] = [1/4, 1/5, 1/6, 1/7, 1/8]
	A[3] = [1/5, 1/6, 1/7, 1/8, 1/9]
	A[4] = [1/6, 1/7, 1/8, 1/9, 1/10]
	
	b = np.array([.882, .744, .618, .521, .447]).transpose()
	
	x = solve(A, b)
	
	print(x)
	
	K_A = np.linalg.cond(A)
	print('Condition number %f' % K_A)

	xh = np.array([-2.1333, 0.6258, 17.4552, -11.8692, -1.4994]).transpose()

	print('norm of left %f' % (np.linalg.norm(x-xh)/np.linalg.norm(x)) )
	print('norm of right %f ' % (K_A * np.linalg.norm(b - A.dot(xh))/np.linalg.norm(b)) )
	
	return x, np.linalg.cond(A)
