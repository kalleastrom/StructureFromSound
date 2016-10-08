function [res,jac]=calcresandjac(D,I,J,x,y,o);

nn = length(D);
m = size(x,2);
n = size(y,2);
V = x(:,I)-y(:,J);
Vt = V';
dd = sqrt(sum(V.^2,1))';
idd = 1./dd;
res = dd+o(J)'-D;
II = (1:length(I))';
JJ1 = (I-1)*3+1;
JJ2 = (I-1)*3+2;
JJ3 = (I-1)*3+3;
JJ4 = (J-1)*4+1+3*m;
JJ5 = (J-1)*4+2+3*m;
JJ6 = (J-1)*4+3+3*m;
JJ7 = (J-1)*4+4+3*m;

VV1 = idd.*Vt(:,1);
VV2 = idd.*Vt(:,2);
VV3 = idd.*Vt(:,3);
VV4 = -idd.*Vt(:,1);
VV5 = -idd.*Vt(:,2);
VV6 = -idd.*Vt(:,3);
VV7 = ones(size(idd));

jac = sparse([II;II;II;II;II;II;II],[JJ1;JJ2;JJ3;JJ4;JJ5;JJ6;JJ7],[VV1;VV2;VV3;VV4;VV5;VV6;VV7],nn,3*m+4*n);
