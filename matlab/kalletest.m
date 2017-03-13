%
if 0,
    matches = result.matches;
    rns = scores;
    cmatches = result.matches;
    cres = result.res;
    settings = result.settings;
    [dmatches,dres] = find_more_inliers_uij(matches,rns,cmatches,cres,settings);
end

load allt20170217
mmm = [ ...
    0 , -0.22105932 , -0.31131638 , -0.67200956 , -1.2055209 , -1.4236101 , -1.8775508 , -1.9329987 ; ...
    0 , -0.01974818 , 0.63091298 , 0.87860974 , 0.78317116 , -0.21526238 , 1.110223e-16 , 0.87813412 ; ...
    0 , -1.3555765 , 2.220446e-16 , -1.1690118 , -0.049042016 , -1.2436546 , -1.110223e-16 , -1.291896 ];
dres.x = mmm;
[dmatches,dres] = find_more_inliers_uij(matches,rns,dmatches,dres,settings);

%%

dmatches.uinliers = dmatches.uinliers & ( abs(dres.resm)<0.02);
inl = dmatches.uinliers;
u = dmatches.u;
[I,J]=find(inl);
D = u(find(inl));
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
figure(1); clf;
plot(res,'.');

myplot(rns,dmatches,dres,result.settings);

%%
myplot(scores,result.matches,result.res,result.settings);
