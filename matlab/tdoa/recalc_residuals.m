function [cres] = recalc_residuals(cmatches,cres);

x = cres.x;
y = cres.y;
o = cres.o;

utimes = cmatches.utimes;
u = cmatches.u;

dd = reshape(sqrt( sum( (kron(ones(1,size(y,2)),x) - ...
    kron(y,ones(1,size(x,2))) ).^2 , 1 ) ),size(x,2),size(y,2)) + repmat(o,size(x,2),1);

[I,J]=find(cmatches.uinliers);
D = u(find(cmatches.uinliers));
[res,Dproj]=calcres(D,I,J,x,y,o);
resm = full(sparse(I,J,res));
cres.res = res;
cres.resm = resm;
