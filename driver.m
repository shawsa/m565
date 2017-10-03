B = eye(21,21) - tril(ones(21,21),-1)
B(:,21)=1
[A,p,g] = lupp_growth(B)