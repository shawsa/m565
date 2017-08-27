'''

sudo code for algorithm

a = 1;
b = 1/root(2)
t = 1/4;
j = 1;
while |a − b| ≥ ε do
y = a;
a = (a + b)/2;
b =
√
by
t = t − j(a − y)
2
;
j = 2j;
pi_approx = a
2/t;
end while
return pi_approx;

'''


