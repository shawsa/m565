import numpy as np
from scipy.linalg import *

def foo_3b_bad(x,A,B,C,T):
    k2 = inv(B.transpose()).dot(x)
    k2 = x.transpose().dot(k2)
    
    k1 = inv(A).dot(x)
    k1 = inv(T).dot(k1)
    k1 = C.dot(k1)
    k1 = inv(B).dot(k1)
    k1 = x.transpose().dot(k1)
    
    return k1 + k2

def foo_3b(x,A,B,C,T):
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
    k1 = solve_triangular(U, k1, lower=False)
    k1 = P.dot(k1)
    k1 = solve_triangular(L, k1, lower=True)
    k1 = x.transpose().dot(k1)
    
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
    return x,A,B,C,T

