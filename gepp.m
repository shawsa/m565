function [x] = gepp(A,b)
n = size(A,1);
[LU,p] = lupp(A);
new_b = b;
P = zeros(n,n);
for i=1:n
   P(p(i),i)=1;
end
x = backsub(LU,forsub(LU,P'*b));
end