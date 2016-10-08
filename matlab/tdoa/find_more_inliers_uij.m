function [dmatches,dres] = find_more_inliers_uij(rawmatches,rns,cmatches,cres,settings);
% [dmatches,dres] = find_more_inliers_uij(rawmatches,cmatches,cres,settings);

x = cres.x
u = rawmatches.u;
[y,o,inl]=tdoa_trilateration_y(u,x);
[I,J]=find(inl);
D = u(find(inl));
[res,Dproj]=calcres(D,I,J,x,y,o);
resm = full(sparse(I,J,res));
figure(6);
plot(res,'.');
ok = abs(res)<0.05;
okm = full(sparse(I,J,ok));
okj = (sum(okm)>=5);
fokj = find(okj);
% Packa in i dmatches
% x = dres.x;
% y = dres.y;
% o = dres.o;
dmatches.uij = rawmatches.uij(:,fokj);
dmatches.u = rawmatches.u(:,fokj);
dmatches.uindex = rawmatches.uindex(:,fokj);
dmatches.utimes = rawmatches.utimes(:,fokj);
dmatches.uok = rawmatches.uok(:,fokj);
dmatches.uinliers = okm(:,fokj);
inl = dmatches.uinliers;
u = dmatches.u;
[I,J]=find(inl);
D = u(find(inl));
y = y(:,fokj);
o = o(:,fokj);

% Packa in i dres
dres.x = x;
dres.y = y;
dres.o = o;
dres.xinorm = settings.xinorm;
myplot(rns,dmatches,dres,settings);
x = real(dres.x);
y = real(dres.y);
o = real(dres.o);
%[xx4,yy4,oo4,res,jac]=bundletdoa(D,I,J,x,y,o,1,dres.xinorm);
[xx4,yy4,oo4,res,jac]=bundletdoasmooth(D,I,J,x,y,o,1,dres.xinorm);
% bundletdoasmooth(D,I,J,xt,yt,ot,debug,xinorm)
dres.x = xx4;
dres.y = yy4;
dres.o = oo4;
[dres] = recalc_residuals(dmatches,dres);
%myplot(rns,dmatches,dres,settings);
