import numpy as np
from scipy.linalg import *

def foo_3b_bad(x,A,B,C,T):
    K_2 = inv(B.transpose()).dot(x)
    K_2 = x.transpose().dot(K_2)
    
    K_1 = inv(A).dot(x)
    K_1 = inv(T).dot(K_1)
    K_1 = C.dot(K_1)
    K_1 = inv(B).dot(K_1)
    K_1 = x.transpose().dot(K_1)
    
    return K_1 + K_2

def foo_3b(x,A,B,C,T):
    
    P, L, U = lu(B)
    y = 
    
    
    
    return P
    
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

