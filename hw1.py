from pytex import *

from math import sqrt, pi, log, atan, factorial, exp, floor
from decimal import *

import matplotlib.pyplot as plt
import numpy as np

#pi to various degrees of accuracy is hard-coded in order to measure error.

pi_102 = Decimal('3.14159265358979323846264338327950288419716939937510582097494459' + '23078164062862089986280348253421170679')

pi_1002 = Decimal('3.1415926535 8979323846 2643383279 5028841971 6939937510 5820974944 5923078164 0628620899 8628034825 3421170679 8214808651 3282306647 0938446095 5058223172 5359408128 4811174502 8410270193 8521105559 6446229489 5493038196 4428810975 6659334461 2847564823 3786783165 2712019091 4564856692 3460348610 4543266482 1339360726 0249141273 7245870066 0631558817 4881520920 9628292540 9171536436 7892590360 0113305305 4882046652 1384146951 9415116094 3305727036 5759591953 0921861173 8193261179 3105118548 0744623799 6274956735 1885752724 8912279381 8301194912 9833673362 4406566430 8602139494 6395224737 1907021798 6094370277 0539217176 2931767523 8467481846 7669405132 0005681271 4526356082 7785771342 7577896091 7363717872 1468440901 2249534301 4654958537 1050792279 6892589235 4201995611 2129021960 8640344181 5981362977 4771309960 5187072113 4999999837 2978049951 0597317328 1609631859 5024459455 3469083026 4252230825 3344685035 2619311881 7101000313 7838752886 5875332083 8142061717 7669147303 5982534904 2875546873 1159562863 8823537875 9375195778 1857780532 1712268066 1300192787 6611195909 2164201989'.replace(' ',''))

pi_2002 = Decimal('3.14159265358979323846264338327950288419716939937510'+
'58209749445923078164062862089986280348253421170679'+
'82148086513282306647093844609550582231725359408128'+
'48111745028410270193852110555964462294895493038196'+
'44288109756659334461284756482337867831652712019091'+
'45648566923460348610454326648213393607260249141273'+
'72458700660631558817488152092096282925409171536436'+
'78925903600113305305488204665213841469519415116094'+
'33057270365759591953092186117381932611793105118548'+
'07446237996274956735188575272489122793818301194912'+
'98336733624406566430860213949463952247371907021798'+
'60943702770539217176293176752384674818467669405132'+
'00056812714526356082778577134275778960917363717872'+
'14684409012249534301465495853710507922796892589235'+
'42019956112129021960864034418159813629774771309960'+
'51870721134999999837297804995105973173281609631859'+
'50244594553469083026425223082533446850352619311881'+
'71010003137838752886587533208381420617177669147303'+
'59825349042875546873115956286388235378759375195778'+
'18577805321712268066130019278766111959092164201989'+
'38095257201065485863278865936153381827968230301952'+
'03530185296899577362259941389124972177528347913151'+
'55748572424541506959508295331168617278558890750983'+
'81754637464939319255060400927701671139009848824012'+
'85836160356370766010471018194295559619894676783744'+
'94482553797747268471040475346462080466842590694912'+
'93313677028989152104752162056966024058038150193511'+
'25338243003558764024749647326391419927260426992279'+
'67823547816360093417216412199245863150302861829745'+
'55706749838505494588586926995690927210797509302955'+
'32116534498720275596023648066549911988183479775356'+
'63698074265425278625518184175746728909777727938000'+
'81647060016145249192173217214772350141441973568548'+
'16136115735255213347574184946843852332390739414333'+
'45477624168625189835694855620992192221842725502542'+
'56887671790494601653466804988627232791786085784383'+
'82796797668145410095388378636095068006422512520511'+
'73929848960841284886269456042419652850222106611863'+
'06744278622039194945047123713786960956364371917287'+
'46776465757396241389086583264599581339047802759010')



# problem 1 **********************************************************************************


def pi_approx(epsilon):
    pi_approx = []

    a = 1
    b = 1/sqrt(2)
    t = 1/4
    j = 1
    while abs(a-b) >= epsilon:
        y = a
        a , b = (a+b)/2 , sqrt(b*y)
        t -= j*(a-y)**2
        j *= 2
        pi_approx.append(a**2/t)
    return pi_approx
    
def pi_approx_dec(digits):
    getcontext().prec = digits+2
    epsilon = Decimal(10)**(-1*digits)
    pi_approx = []
    a = Decimal(1)
    b = 1 / Decimal(2).sqrt()
    t = Decimal(1/4)
    j = Decimal(1)
    while abs(a-b) >= epsilon:
        y = a
        a , b = (a+b)/2 , (b*y).sqrt()
        t -= j*(a-y)**2
        j *= 2
        pi_approx.append(a**2/t)
    return pi_approx

def print_hw_output_1b(epsilon=1E-14):
    print('(b)')
    pa = pi_approx(epsilon)
    #latex format
    print('\t\\begin{center}')
    print('\t\t\\begin{tabular}{|c|c|c|}')
    print('\t\t\t\\hline')
    #print('itter', 'approx', 'abs error')
    print('\t\t\t$n$ & $\\tilde{\\pi}$ & $\\vert \\tilde{\\pi} - \\pi \\vert$ \\\\ \\hline')
    for i, p in enumerate(pa):
        #print(i+1, p, abs(p - pi))\
        print('\t\t\t' + str(i+1) + ' & ' + str(p) + ' & ' + str(abs(p - pi)) + ' \\\\ \\hline')
    print('\t\t\\end{tabular}')
    print('\t\\end{center}')
    
def print_hw_output_1c(epsilon=1E-14):
    pa = pi_approx(epsilon)
    n = []
    error = []
    for i, p in enumerate(pa):
        n.append(i+1)
        error.append(abs(p - pi))
    plt.semilogy(n, error, 'bo')
    plt.xlabel('n')
    plt.ylabel('absolute error')
    plt.xticks(np.arange(1,len(n)+1,1))
    plt.show()
    
def print_hw_output_1e(epsilon=2000):
    pa = pi_approx_dec(epsilon)
    print('n', 'abs error exponent')
    for i,p in enumerate(pa):
        print(i, (abs(p-pi_2002)).log10().quantize(Decimal('.00001')) )
    #latex format
    print('\t\\begin{center}')
    print('\t\t\\begin{tabular}{|c|c|}')
    print('\t\t\t\\hline')
    #print('itter', 'approx', 'abs error')
    print('\t\t\t$n$ & exponent of absolute error \\\\ \\hline')
    for i, p in enumerate(pa):
        #print(i+1, p, abs(p - pi))\
        exponent_error = (abs(p-pi_2002)).log10().quantize(Decimal('.00001'))
        print('\t\t\t' + str(i+1) + ' & ' + str(exponent_error) + ' \\\\ \\hline')
    print('\t\t\\end{tabular}')
    print('\t\\end{center}')
    
    
# problem 2 **********************************************************************************

def e_taylor(x, n):
    partial_sum = 0
    for i in range(0,n+1):
        partial_sum += x**i / factorial(i)
    return partial_sum
    
def e_better_taylor(x,n):
    ret = e_taylor(x/1024,n)
    for i in range(10):
        ret *= ret
    return ret
    
def e_cos(theta,n):
    if theta < 0:
        theta = theta * -1
        
    pi_sub = floor(theta/(2*pi))
    theta = theta - pi*2*pi_sub
    return e_better_taylor(theta * 1j, n).real
    
def e_sin(theta,n):
    sign = 1
    
    if theta < 0:
        theta = theta * -1
        sign = -1
        
    pi_sub = floor(theta/(2*pi))
    theta = theta - pi*2*pi_sub
    return sign * e_better_taylor(theta * 1j, n).imag
    
def print_hw_output_2a():
    xs = [0.5, -0.5, 30*pi, -30*pi] * 2
    ns = [10]*4 + [40]*4
    
    for i, (x, n) in enumerate(zip(xs,ns)):
        print('{:.4f}'.format(x), n, (exp(x) - e_taylor(x,n)) / exp(x))
        
    x_labels = ['$0.5$', '$-0.5$', '$30\pi$', '$-30\pi$'] * 2  
    rel_errors = [''] * len(x_labels)
    for i, x in enumerate(xs):
        rel_error = abs( (exp(x) - e_taylor(x,ns[i])) / exp(x) )
        rel_errors[i] = "{:.5e}".format(rel_error)
        
    latex_table( (x_labels, ns, rel_errors), ('$x$', '$n$', 'relative error'))
    
def print_hw_output_2b():
    xs = [0.5, -0.5, 30*pi, -30*pi] * 2
    ns = [10]*4 + [40]*4
    
    for i, (x, n) in enumerate(zip(xs,ns)):
        print('{:.4f}'.format(x), n, (exp(x) - e_taylor(x,n)) / exp(x))
        
    x_labels = ['$0.5$', '$-0.5$', '$30\pi$', '$-30\pi$'] * 2
    rel_errors = [''] * len(x_labels)
    for i, x in enumerate(xs):
        rel_error = abs( (exp(x) - e_better_taylor(x,ns[i])) / exp(x) )
        rel_errors[i] = "{:.5e}".format(rel_error)
        
    latex_table( (x_labels, ns, rel_errors), ('$x$', '$n$', 'relative error'))
    
# problem 5 **********************************************************************************
        
        
def print_hw_output_5c(N=10**8):
    print('arctan(N+1) - arctan(N) = ' + str(atan(N+1) - atan(N)) )
    print('arctan(1/(N**2+N+1) = ' + str( atan(1/(N**2+N+1)) ) )
        
        
# problem 6 **********************************************************************************

# failed to reproduce the error in python. Problem 6 is done in C code instead.
def cascsum(p):
    e = 0
    s = p[0]
    for i in range(1, len(p)):
        x = p[i] + s
        z = x - p[i]
        y = (p[i] - (x-z)) + (s-z)
        e += y
        s = x
    s += e
    return s
    

# problem 7 **********************************************************************************

def division(a, x_0=.1, epsilon=10**-15):
    x = x_0
    xs = [x_0]
    n = 0
    ns = [n]
    while abs( (x-(1/a)) *a ) > epsilon:
        n += 1
        x *= (2-a*x)
        ns.append(n)
        xs.append(x)
        
    return (ns,xs)
    
def print_hw_output_7():
    (ns, xs) = division(5)
    latex_table((ns,xs), ('$n$', '$x_n$'))
    





