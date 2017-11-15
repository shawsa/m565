import numpy as np
import scipy.linalg as la
from pytex import *

def p1():
    cond = []
    est = []
    ratio = []
    for n in range(5,20):
        H = la.hilbert(n)
        K = np.linalg.cond(H,p=np.inf)
        cond.append(K)
        E = (1+np.sqrt(2))**(4*n) / np.sqrt(n)
        est.append( E )
        ratio.append (K/E)
    latex_table((range(5,len(cond)+5), cond, est, ratio), ('$n$','$K(H_n)$', 'est', '$K(H_n)/$est'))
        