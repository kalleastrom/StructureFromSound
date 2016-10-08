function ysols=toa_trilateration_one_point(d,x);
% ysols=toa_trilateration_one_point(d,x);
% does trilateration using a vector of distances d
% assuming known positions x 
% and calculates y so that 
% Assumes column vectors d and x
% Assumes minimal case, i e sizes of d and x equal
% Gives the two solutions

[x_dim,m] = size(x);
d_dim = size(d,1);

% Assumes x_dim == m
% Assumes d_dim == m

if x_dim ~= d_dim,
  error('toa_trilateration_one_point only solves minimial case 1 size(x,1) = size(d,1)');
end;
if x_dim ~= m,
  error('toa_trilateration_one_point only solves minimial case 2 size(x,1) = size(x,2)');
end;

%keyboard;

d21 = d(1).^2;
d2lin = d(2:end).^2 - d21;

xbasis = x(:,2:end)-repmat(x(:,1),1,m-1);
x0 = x(:,1);
x2basis = sum(x(:,2:end).^2,1);
x20 = sum(x0.^2);

% 1. solve for the x_dim-1 linear equations first
% in subspace spanned by the points in x
% Linear equations d2lin = x2basis'-x20' -2*(xbasis')*y
% Solve first for y = x0 + xbasis*lambda
% -2*(xbasis')*(x0 + xbasis*lambda) = d2lin-(x2basis-x20);
% solve for lambda 
% -2*(xbasis')*xbasis* lambda  = d2lin-(x2basis-x20) + 2*(xbasis')*x0;
A = -2*(xbasis.')*xbasis;
% There could be a problem here if the 
% elements are complex. 
b = d2lin-(x2basis-x20)' + 2*(xbasis.')*x0;
%b = d2lin-(x2basis-x20) + (xbasis.')*x0 + x0'*xbasis;
lambda = A\b;

y0 = x0 + xbasis*lambda;

% 2. Now solve for the only remaining quadratic equation
% d21 = x20 - 2*x0'*y + y'*y
% where y = y0 + xperp*my
xperp = null(xbasis');
% d21 = x20 - 2*x0'*(y0+xperp*my) + (y0+xperp*my)'*(y0+xperp*my)
% 0 = (-d21 + x20 - 2*x0'*y0 + y0'*y0) + 
%     (-2*x0'*xperp + 2*y0'*xperp)*my + 
%     (xperp'*xperp)*my^2
% By clever choice of solutions the linear term i zero
% thus my = plus minus sqrt(
my = sqrt( - (-d21 + x20 - 2*x0.'*y0 + y0.'*y0) / (xperp.'*xperp) );
ysols = [ (y0+xperp*my) (y0 -xperp*my)];
