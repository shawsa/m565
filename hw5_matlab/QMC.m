format long

rng(11032016);
N = 10000;
points = net(haltonset(6),N);
col_sum = sum(points,2);

I = sum(sin(pi/2*col_sum))/N + 512/pi^6

n = 10;
h = 1/n;
d = 6;
% Trapezoidal rule weights in one dimension
w1 = h*[0.5;ones(n-1,1);0.5];
w = w1;
% Generate the trapezoidal rule weights in 6 dimensions
for j=1:d-1
w = bsxfun(@times,w,reshape(w1,[ones(1,j) n+1]));
end
% Grid of points
[x1,x2,x3,x4,x5,x6] = ndgrid(linspace(0,1,n+1));
g = w.*f(x1,x2,x3,x4,x5,x6);
Qtrap = sum(g(:));

Qtrap + 512/pi^6

function [y] = f(x1,x2,x3,x4,x5,x6)
y = sin(pi/2*(x1+x2+x3+x4+x5+x6));
end
