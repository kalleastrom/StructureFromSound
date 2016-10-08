function [res,Dproj]=calcres(D,I,J,x,y,o);

nn = length(D);
m = size(x,2);
n = size(y,2);
V = x(:,I)-y(:,J);
Vt = V';
dd = sqrt(sum(V.^2,1))';
res = dd+o(J)'-D;
Dproj = dd+o(J)';
