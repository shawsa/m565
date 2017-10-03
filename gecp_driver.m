A = rand(5,5);
b = rand(5,1);
[LU,p,q] = lucp(A);
new_b = b;
P = zeros(5,5);
for i=1:5
   P(p(i),i)=1;
end
Q = zeros(5,5);
for i=1:5
   Q(i,q(i))=1;
end

A\b - Q'*backsub(LU,forsub(LU,P'*b))