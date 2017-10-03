function [x] = gecp(A,b)
n = size(A,1)
   
[LU,p,q] = lucp(A);
new_b = b;
P = zeros(n,n);
for i=1:n
   P(p(i),i)=1;
end
Q = zeros(n,n);
for i=1:n
   Q(i,q(i))=1;
end

x = Q'*backsub(LU,forsub(LU,P'*b))

%A\b - Q'*backsub(LU,forsub(LU,P'*b))

end