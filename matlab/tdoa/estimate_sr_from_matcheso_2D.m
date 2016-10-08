function [cmatches,cres] = estimate_sr_from_matcheso_2D(bmatches,o,bres,settings);


cmatches = bmatches;
cres = bres;

% Remove offsets that would cause d = u-o to be negative
%keyboard
badj = find(min(cmatches.u)-o < 0.1);
cmatches.u(:,badj)=[];
cmatches.uok(:,badj)=[];
cmatches.uindex(:,badj)=[];
cmatches.utimes(:,badj)=[];
cmatches.uinliers(:,badj)=[];

cres.A(:,badj-1)=[];
cres.o(badj)=[];
o(badj)=[];

if 1, % Checking for errors
    [d,dc,d2,d2c,d2r,d2rc]=checkequations(cmatches.u,cmatches.uinliers,o,cres.A);
    [norm(d) norm(dc) norm(dc-d); ...
        norm(d2) norm(d2c) norm(d2-d2c); ...
        norm(d2r) norm(d2rc) norm(d2r-d2rc)];
end

dtmp = cmatches.u - repmat(o,size(cmatches.u,1),1);
[sol,stats,bb,bb0] = toa_linear_fix_2D(dtmp,cres.A);
xx = sol.x;
yy = sol.y;
oo = o;
if 1, % Checking for errors
    [d,dc,d2,d2c,d2r,d2rc]=checkequations(cmatches.u,cmatches.uinliers,oo,cres.A,xx,yy);
    [norm(d) norm(dc) norm(dc-d); ...
        norm(d2) norm(d2c) norm(d2-d2c); ...
        norm(d2r) norm(d2rc) norm(d2r-d2rc)]
end

[bb bb0 d2r(2:end,1) d2rc(2:end,1)]
norm(bb+d2r(2:end,1))
norm(bb0+d2rc(2:end,1))

% Choose coordinates to normalise with respect to
if ~isfield(settings,'xinorm'),
    settings.xinorm = [1 7 3];
end
xinorm = settings.xinorm;
bbb = - xx(:,xinorm(1));
xx = xx+repmat(bbb,1,size(xx,2));
yy = yy+repmat(bbb,1,size(yy,2));
AAA = xx(:,xinorm([2 3]));
[qq,rr]=qr(AAA);
xx = qq'*xx;
yy = qq'*yy;


[I,J,blubb]=find(cmatches.uinliers);
TT = cmatches.u(find(cmatches.uinliers));
[res00,TTproj]=calcres(TT,I,J,xx,yy,oo);
%optionsd

[xx2,yy2,oo2,res,jac]=bundletdoadiffdim(TT,I,J,xx,yy,oo);

[xx2,yy2,oo2,res,jac]=bundletdoadiffdim(TT,I,J,xx,yy,oo);
tmpid = find(yy2(3,:)<0);
yy2(3,tmpid)=-yy2(3,tmpid);

[xx3,yy3,oo3,res,jac]=bundletdoadiffdim(TT,I,J,xx2,yy2,oo2);
%[xx3,yy3,oo3,res,jac]=bundletdoa(TT,I,J,xx3,yy3,oo3,1,xinorm);
std(res)

if 1, % Checking for errors
    [d,dc,d2,d2c,d2r,d2rc]=checkequations(cmatches.u,cmatches.uinliers,oo3,cres.A,xx3,yy3);
end


if settings.doplot,
    figure(3);
    plot(res*1000,'b.');
    resm = sparse(I,J,res);
    %figure(4);
    %plot(resm'*1000,'b.');
    figure(4); clf;
    hist(res(:)*settings.sr/settings.v,100);
    title(['Histogram of errors in terms of sampling points']);
    figure(5); clf;
    hist(res(:)*1000,100);
    title(['Histogram of errors in millimeter']);
    %plot(resm'*1000,'b.');
end;


%% Remove outliers after bundling


if 0,
    % remove further outliers
    [maxv,maxi]=max(sum(abs(resm)));
    oo3(maxi)=[];
    yy3(:,maxi)=[];
    cmatches.u(:,maxi)=[];
    cmatches.uok(:,maxi)=[];
    cmatches.uindex(:,maxi)=[];
    cmatches.utimes(:,maxi)=[];
    cmatches.uinliers(:,maxi)=[];
    utimes(maxi)=[];
    uindex(maxi)=[];
    [I,J,blubb]=find(cmatches.uinliers);
    TT = cmatches.u(find(cmatches.uinliers));
    [xx3,yy3,oo3,res,jac]=bundletdoa(TT,I,J,xx3,yy3,oo3,1,xinorm);
    resm = sparse(I,J,res);
    if settings.doplot,
        figure(3); clf;
        plot(res*1000,'b.');
        resm = sparse(I,J,res);
        figure(4); clf;
        hist(res(:)*settings.sr/settings.v,100);
        title(['Histogram of errors in terms of sampling points']);
        figure(5); clf;
        hist(res(:)*1000,100);
        title(['Histogram of errors in millimeter']);
        %plot(resm'*1000,'b.');
    end;
end

%[maxv,maxi]=max(sum(abs(resm)));

dd3 = reshape(sqrt( sum( (kron(ones(1,size(yy3,2)),xx3) - ...
    kron(yy3,ones(1,size(xx3,2))) ).^2 , 1 ) ),size(xx3,2),size(yy3,2)) + repmat(oo3,size(xx3,2),1);

cres.x = xx3;
cres.y = yy3;
cres.o = oo3;
cres.d_calc = dd3;
cres.u_calc = dd3-repmat(dd3(1,:),settings.mm,1);
cres.res = res;
cres.jac = jac;

