%RBFFIT Fits the thin plate spline radial basis function (RBF) to the data (x(i),y(i),f(i)).
%    RBFFIT(x,y,f) finds the coefficients of the thin-plate spline RBF that 
%    interpolates the data (x(i),y(i),f(i)), where x and y are the independet
%    variables.  The thin-plate spline is defined as
%                          phi(r) = (r^2)*log(r)
%
%   Inputs: 
%          x : x coordinates of the given nodes.
%          y : y coordinates of the given nodes.
%          f : Scalar function values to interpolate at each (x(i),y(i)).
%
%   Outputs: 
%        lam : the coefficients (weights) for the 2-D thin plate spline
%              interpolant.
%
%  See also: rbfval.
function lam = rbffit(x,y,f)

sz = size(x);

% Flatten x, y, and f;
x = x(:); y = y(:); f = f(:);
n = length(x);

if n ~= length(y) || n ~=length(f)
   error('The sizes of the input vectors must be the same');
end

%
% Compute the pairwise distances between all the nodes.
%

% Prepare the distance squared matrix.
A = zeros(n);

[xd1,xd2] = meshgrid(x);  % x coordinates
A = A + (xd1-xd2).^2;

[xd1,xd2] = meshgrid(y);  % y coordinates
A = A + (xd1-xd2).^2;

clear xd1 xd2;

% Need to apply the TPS radial function to A.  However, the TPS is r^2*log(r)
% and will have numerical trouble when r=0.  It should be 0 when r=0, so we
% explicity make this the case.
id = 1:(n+1):n^2;    % indicies for the diagonal of A (i.e. where r=0).
A(id) = 1;           % make the distance 1 since log(1) = 0.
% Compute the TPS interpolation matrix.  Note that we are dealing with the
% distances squared.  The TPS kernel with distances squared reduces to 
% 0.5*A*log(A)
A = 0.5*A.*log(A);

% Add on the extra polynomial conditions.
ev = ones(n,1);
A = [[A ev x y];[[ev x y]' zeros(3)]];

% Solve for the RBF expansion coefficients
lam = A\[f;zeros(3,1)];


