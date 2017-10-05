n = 12;
A = zeros(n,n);
b = zeros(n,1);
for i=1:n
    b(i,1) = i^(n-1);
    for j=1:n
        A(i,j) = i^(j-1);
    end
end
x = zeros(n,1);
x(n,1) = 1;

xhatcp = gecp(A,b);
xhatpp = gepp(A,b);

'xhatCP residual'
norm(b-A*xhatcp,Inf)
'xhatCP error'
norm(x-xhatcp,Inf)

'xhatPP residual'
norm(b-A*xhatpp,Inf)
'xhatPP error'
norm(x-xhatpp,Inf)

A\b