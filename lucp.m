function [A,p,q] = lucp(A)
%lupp_loop  Computes the LU decomposition of A with partial pivoting
%
%   [A,p] = lupp(A) Computes the LU decomposition of A with partial
%   pivoting using vectorization. The factors L and U are returned in the
%   output A, and the permutation of the rows from partial pivoting are
%   recorded in the vector p. Here L is assumed to be unit lower
%   triangular, meaning it has ones on its diagonal. The resulting
%   decomposition can be extracted from A and p as follows:
%       L = eye(length(LU))+tril(LU,-1);     % L with ones on diagonal
%       U = triu(LU);                        % U 
%       P = p*ones(1,n) == ones(n,1)*(1:n);  % Permutation matrix
%   A is then given as L*U = P*A, or P'*L*U = A.
%
%   Use this function in conjuction with backsub and forsub to solve a
%   linear system Ax = b.

n = size(A,1);
p = (1:n)';
q = (1:n);

for k=1:n-1
    % Find the entry that contains the largest entry in magnitude
    [temp, col_max] = max(abs(A(k:n,k:n)));
    [temp, x_pos] = max(col_max);
    y_pos = col_max(x_pos);
    %pos = pos +1
    %x_pos = mod(pos,n)+1
    %y_pos = floor(pos/n)+1;
    row2swap = k-1+x_pos
    col2swap = k-1+y_pos
    % Swap the rows and columns of A and p (perform complete pivoting)
    A([row2swap, k],:) = A([k, row2swap],:);
    p([row2swap, k]) = p([k, row2swap]);
    A(:,[col2swap,k]) = A(:,[k,col2swap]);
    q(:,[col2swap,k]) = q(:,[k,col2swap]);
    % Perform the kth step of Gaussian elimination
    J = k+1:n;
    A(J,k) = A(J,k)/A(k,k);
    A(J,J) = A(J,J) - A(J,k)*A(k,J);
end

end
