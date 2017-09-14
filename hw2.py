import numpy as np






# Problem 1     ********************************************************************

def secant(f, x0, x1, epsilon=1e-16, n_max=10**9):
	xs = [x0,x1]
	xn = x1
	xn_1 = x0
	fn = f(xn)
	fn_1 = f(xn_1)
	for i in range(n_max):
		if abs(fn)<epsilon:
			return xn, i
		#print(xn)
		xn, xn_1 = xn - (xn-xn_1)/(fn-fn_1)*fn, xn
		fn, fn_1 = f(xn), f(xn_1)
	return None

def newton(f, df, x0, epsilon=1e-16, n_max=10**9):
	x = x0
	fx = f(x)
	dfx = df(x)
	for i in range(n_max):
		if abs(fx) < epsilon:
			return x, i
		x = x - fx/dfx
		fx, dfx = f(x), df(x)
	return None

def foo_p1(x):
	return np.exp(-1*x)-x

def dfoo_p1(x):
	return -1*np.exp(-1*x)

def hw2_p1c():
	
