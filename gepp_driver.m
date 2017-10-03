A = rand(5,5);
b = rand(5,1);
[LU,p] = lupp(A);
P = zeros(5,5);
for i=1:5
   P(p(i),i)=1
end
inv(A)*b - backsub(LU,forsub(LU,P'*b))