function [ysols,osols]=tdoa_trilateration_y_one_point(u,x);
% ysols=tdoa_trilateration_y_one_point(d,x);
% does trilateration using a vector of distances d
% assuming known positions x 
% and calculates position y and offset o so that 
% d = |x-y|+o
% Assumes column vectors d and x
% Assumes minimal case, i e sizes of d and x equal
% Gives the two solutions

[x_dim,m] = size(x);
u_dim = size(u,1);

% Assumes x_dim == m
% Assumes d_dim == m

if x_dim+1 ~= u_dim,
  error('tdoa_trilateration_y_one_point only solves minimial case 1 size(x,1)+1 = size(d,1)');
end;
if x_dim+1 ~= m,
  error('tdoa_trilateration_y_one_point only solves minimial case 2 size(x,1)+1 = size(x,2)');
end;

%keyboard;

u21 = u(1).^2;
u2lin = u(2:end).^2 - u21;
x2 = sum( x.^2);

%% linear equations
% (ui-o)^2-(u1-o)^2 = ui^2-2*ui*o - u1^2 + 2*u1*o = 
%   xi^2-2*xi'*y - x1^2 + 2*x1'*y = 
% = (ui.^2-u1.^2-xi'*xi+x1'*x1)*1 + (-2*ui+2*u1)*o + (2*xi' -2*x1')*y

MM = [ (u2lin-x2(2:end)' + x2(1)*ones(x_dim,1)) ...
    (-2*u(2:end)+2*u(1)*ones(x_dim,1)) ...
    (2*x(:,2:end)' - 2*ones(x_dim,1)*x(:,1)')];
AA = MM(:,2:end);
bb = -MM(:,1);

% we now have AA*[o;y]=bb;
% using all but the first equation. 
% this has one free parameter

z0 = AA\bb;
zbase = null(AA);

y0 = z0(2:end);
yb = zbase(2:end);
o0 = z0(1);
ob = zbase(1);
u1 = u(1);
x1 = x(:,1);

% so y = z0 + lambda*zbase

% First equations is (u(1)-o)^2 = (x(:,1)-y)'*(x(:,1)-y) = 
% u(1)^2 - 2*u(1)*o + o^2 = x(:,1)'*x(:,1) - 2*x(:,1)'*y + y'*y = 
% u1^2 - 2*u1*(o0+l*ob) + (o0^2 + 2*l*o0*ob + l^2*ob^2) =
%  (u1^2 - 2*u1*o0 + o0^2) + (-2*u1*ob + 2*o0*ob)*l + (ob^2)*l^2
%   x1'*x1 - 2*x1'*(y0+l*yb) + (y0+l*yb)'*(y0+l*yb) = 
%   (x1'*x1 - 2*x1'*y0 + y0'*y0) + (-2*x1'*yb + 2*yb'*y0)*l + 
%    (yb'*yb)*l ^2

pp = [ (u1^2 - 2*u1*o0 + o0^2) (-2*u1*ob + 2*o0*ob) (ob^2)] - ...
[ (x1'*x1 - 2*x1'*y0 + y0'*y0)  (-2*x1'*yb + 2*yb'*y0) (yb'*yb)];

lambdasol = roots(pp(3:-1:1));

ysols = [ (y0+yb*lambdasol(1)) (y0+yb*lambdasol(2))];
osols = [ (o0+ob*lambdasol(1)) (o0+ob*lambdasol(2))];



