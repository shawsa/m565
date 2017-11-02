%RBFVAL Evaluates the 2D thin-plate spline (TPS) radial basis function interpolant 
%       specified by the coefficients lambda at the points ux,uy.
%
% s = RBFVAL(lam,x,y,ux,uy) evaluates the 2D thin plate spline interpolant
% specified by the coefficents lam(j) (obtained from rbffit) and the node 
% locations (x(j),y(j)) at the points (ux(i),uy(i)).
%
%   Inputs: 
%
%        lam : The thin-plate expansion coefficients computed from rbffit.
%          x : x coordinates of the given nodes.
%          y : y coordinates of the given nodes.
%         ux : x coordinates of the evaluation locations.
%         uy : y coordinates of the evaluation locations.
%
%  See also: rbffit
function p = rbfval(lam,x,y,ux,uy)

sz = size(ux);

% Flatten x, y, lam, ux, uy;
x = x(:); y = y(:); lam = lam(:)'; ux = ux(:); uy = uy(:);
n = length(x);
m = length(ux);

% The following is not the most computationally efficient way to do the
% evaluation.  Instead it is designed to be efficient in terms of memory
% use.

% Prepare the distance squared matrix.
B = zeros(n+3,1);
p = zeros(m,1);
for j=1:m
   B(1:n,1) = (x-ux(j)).^2 + (y - uy(j)).^2;
   
   % Need to apply the TPS radial function to B.  However, the TPS is r^2*log(r)
   % and will have numerical trouble when r=0.  It should be 0 when r=0, so we
   % explicity make this the case.
   B = 0.5*B.*log(B+eps);
   
   % Polynomial part.
   B(n+1:n+3,1) = [1;ux(j);uy(j)];
   
   % Evaluate the RBF
   p(j,1) = lam*B;
end

p = reshape(p,sz);



