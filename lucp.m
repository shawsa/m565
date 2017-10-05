function [A,p,q] = lucp(A)
n = size(A,1);
p = (1:n)';
q = (1:n);

for k=1:n-1
    % Find the entry that contains the largest entry in magnitude
    [temp, col_max] = max(abs(A(k:n,k:n)));
    [temp, x_pos] = max(col_max);
    y_pos = col_max(x_pos);
    row2swap = k-1+x_pos;
    col2swap = k-1+y_pos;
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
