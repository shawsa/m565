'''

sudo code for algorithm

a = 1;
b = 1/root(2)
t = 1/4;
j = 1;
while |a − b| ≥ ε do
y = a;
a = (a + b)/2;
b = sqrt(by)
t = t − j(a − y)**2
j = 2j;
pi_approx = a
2/t;
end while
return pi_approx;

'''

from math import sqrt, pi
from decimal import *

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
