function [oopt,Aopt,paramopt,resopt,jacopt]=bundle_rank1(D,I,ot,rk,At,debug);
% [oopt,Aopt,paramopt,resopt,jacopt]=bundle_rank1(D,I,ot,rk,At,debug);
% non-linear least squeares optimization to find
% matrix A with rank rk and offsets o so that
% the matrix Cl'*(D.^2-2*D*diag(offset))*Cr
% is as close to the matrix A in a least squares sense
% The method allows for missing data in D
% The boolean matrix I is used to
% identify which data are present.


if nargin<=5,
    debug = 0;
end;

if nargin<=4,
    % Estimate At from D,ot
    [m,n]=size(D);
    param.m = m;
    param.n = n;
    param.D = D;   
    Cr = [-ones(1,n-1); eye(n-1)];
    Cl = [-ones(1,m-1); eye(m-1)];
    At = Cl'*(D.^2-2*D*diag(ot))*Cr;
    [u,s,v]=svd(At);
    At = u(:,1:rk)*s(1:rk,1:rk)*v(:,1:rk)';
end;

if 0,
    % for debugging
    dz = randn(27,1)*0.001;
    A0 = At,
    o0 = ot;
    p0 = param;
    [o1,A1,p1]=updatexy(o0,A0,dz,p0);
    [res0,jac0]=calcresandjac(o0,A0,p0);
    [res1,jac1]=calcresandjac(o1,A1,p1);
end

if 1,
    %precalculation of certain elements
    [m,n]=size(D);
    param.m = m;
    param.n = n;
    param.D = D;
    
    Cr = [-ones(1,n-1); eye(n-1)];
    Cl = [-ones(1,m-1); eye(m-1)];
    %Tr = Cl'*T*Cr;
    Ir = I(1,1)*ones(m-1,n-1) & ...
        ones(m-1,1)*I(1,2:end) & ...
        I(2:end,1)*ones(1,n-1) & ...
        I(2:end,2:end);
    param.Ir = Ir;
    
    k0 = 0;
    param.i1 = (1:length(ot))+k0;
    k0 = k0+length(ot);
    param.rk = rk;
    [u,s,v]=svd(At);
    param.u0 = u;
    param.v0 = s(1:rk,1:rk)*v(:,1:rk)';
    kkk = 0;
    m1 = m-1;
    n1 = n-1;
    BB = zeros(m1,m1,1);
    BB0 = BB;
    for j = 1:rk;
        for i = (rk+1):m1;
            %[i j]
            kkk = kkk+1;
            BBtmp = BB0;
            BBtmp(i,j)=1;
            BBtmp(j,i)=-1;
            BB(:,:,kkk)=BBtmp;
        end
    end
    param.nBB = kkk;
    param.i2 = (1:param.nBB)+k0;
    k0 = k0+param.nBB;
    param.BB = BB;
    
    EE = eye(size(At,1));
    E0 = EE(:,1:rk);
    param.E0 = E0;
    
    kkk = 0;
    m1 = m-1;
    n1 = n-1;
    CC = zeros(size(param.v0,1),size(param.v0,2),1);
    CC0 = CC;
    for i = 1:size(param.v0,1);
        for j = 1:size(param.v0,2);
            %[i j]
            kkk = kkk+1;
            CCtmp = CC0;
            CCtmp(i,j)=1;
            CC(:,:,kkk)=CCtmp;
        end
    end
    param.nCC = kkk;
    param.i3 = (1:param.nCC)+k0;
    k0 = k0+param.nCC;
    param.CC = CC;
    
    param.o = ot;
    param.A = At;
    param.D = D;
    param.Cr = Cr;
    param.Cl = Cl;
    param.I = I;
    param.Ir = Ir;
    
    %Tr2  = Cl'*(D.^2-2*D*diag(offset))*Cr;
    EEE = eye(n);
    DD0 = Cl'*(D.^2)*Cr;
    DD = zeros(m-1,n-1,n);
    for kk = 1:n;
        DD(:,:,kk) = Cl'*(-2*D*diag(EEE(:,kk)))*Cr;
    end
    Tr3 = DD0;
    for kk = 1:n;
        Tr3 = Tr3 + DD(:,:,kk)*ot(kk);
    end;
    
    param.DD0 = DD0;
    param.DD = DD;
end

%keyboard;
for kkk = 1:10;
    kkk;
    %keyboard;
    [res,jac]=calcresandjac(ot,At,param);
%     if 1,
%         [u,s,v]=svd(jac);
%         log(abs(v(:,end)))'
%         diag(s)'
%         keyboard;
%     end
    %[res,jac]=calcresandjac(D,I,J,xt,yt,ot);
    %dz = -(jac\res);
    %dz = -(jac'*jac+eye(size(jac,2)))\(jac'*res);
    %[u,s,v]=svd(full(jac),0);
    %     [u,s,v]=svd(full(jac),'econ');
    %     u = u(:,1:(end-1));
    %     s = s(1:(end-1),1:(end-1));
    %     v = v(:,1:(end-1));
    %dz = -v*inv(s)*u'*res;
    dz = -jac\res;
    [otn,Atn,paramn]=updatexy(ot,At,dz,param);
    %[xtn,ytn,otn]=updatexy(xt,yt,ot,dz);
    [res2,jac2]=calcresandjac(otn,Atn,paramn);
    aa = [norm(res) norm(res+jac*dz) norm(res2)]
    bb = aa;
    bb=bb-bb(2);
    bb = bb/bb(1);
    cc = norm(jac*dz)/norm(res);
    if debug,
        aa
        cc
    end
    % Check that the error actually gets smaller
    if norm(res)<norm(res2),
        % bad
        % check that the update is sufficiently big
        % otherwise it is maybe just numerical limitations
        if cc>1e-4,
            % OK then it is probably just the linearization that
            % is not good enough for this large step size, decrease
            kkkk = 1 ;
            while (kkkk<50) && (norm(res)<norm(res2)),
                dz = dz/2;
                [otn,Atn,paramn]=updatexy(ot,At,dz,param);
                %[xtn,ytn,otn]=updatexy(xt,yt,ot,dz);
                [res2,jac2]=calcresandjac(otn,Atn,paramn);
                kkkk = kkkk+1;
            end
        end
    end
    if debug,
        aa = [norm(res) norm(res+jac*dz) norm(res2)];
        bb = aa;
        bb=bb-bb(2);
        bb = bb/bb(1);
        cc = norm(jac*dz)/norm(res);
        %keyboard;
        aa
        bb
        cc
    end;
    if norm(res2)<norm(res)
        ot = otn;
        At = Atn;
        param = paramn;
    else
        %disp([num2str(kkk) '  stalled']);
    end
end;

oopt = ot;
Aopt = At;
paramopt = param;
resopt = res;
jacopt = jac;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of main function bundle_rank1.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ony,Any,paramny]=updatexy(o,A,dz,param);

m = param.m;
n = param.n;
dz1 = dz(param.i1);
dz2 = dz(param.i2);
dz3 = dz(param.i3);
ony = o + dz1;
Btmp = zeros(m-1,m-1);
for kk = 1:size(param.BB,3);
    Btmp = Btmp + param.BB(:,:,kk)*dz2(kk);
end
vny = param.v0;
for kk = 1:size(param.CC,3);
    vny = vny + param.CC(:,:,kk)*dz3(kk);
end;
Any = param.u0*expm(Btmp)*param.E0*vny;

paramny = param;
paramny.o = ony;
paramny.Any = Any; 
[u,s,v]=svd(Any);
rk = paramny.rk;
paramny.u0 = u;
paramny.v0 = s(1:rk,1:rk)*v(:,1:rk)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of updatexy


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res,jac]=calcresandjac(o,A,param);

m = param.m;
n = param.n;
o = param.o;
dodz = zeros(length(o),param.i3(end));
dodz(param.i1,param.i1)=eye(length(o));

%res = o;
%jac = dodz;

A = param.u0*param.E0*param.v0;
dAdz = zeros(size(A,1),size(A,2),param.i3(end));
for kk = 1:size(param.BB,3);
    dAdz(:,:,param.i2(kk)) = param.u0*param.BB(:,:,kk)*param.E0*param.v0;
end
for kk = 1:size(param.CC,3);
    dAdz(:,:,param.i3(kk)) = param.u0*param.E0*param.CC(:,:,kk);
end;

%res = reshape(A,35,1);
%jac = reshape(dAdz,35,27)

Tr = param.DD0;
for kk = 1:n;
  Tr = Tr + param.DD(:,:,kk)*o(kk);
end;
dTrdz = zeros(size(A,1),size(A,2),param.i3(end));
for kk = 1:length(o);
  dTrdz(:,:,kk) = param.DD(:,:,kk);
end;

%res = reshape(Tr,35,1);
%jac = reshape(dTrdz,35,27);

res = Tr-A;
res = res(:);
jac = dTrdz-dAdz;
jac = reshape(jac,[size(dAdz,1)*size(dAdz,2) size(dAdz,3)]);
res = res(find(param.Ir(:)));
jac = jac(find(param.Ir(:)),:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of calcresandjac
