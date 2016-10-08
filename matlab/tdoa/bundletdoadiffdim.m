function [xopt,yopt,oopt,res,jac]=bundletdoadiffdim(D,I,J,xt,yt,ot,options,debug);
% [xopt,yopt]=bundletdoadiffdim(D,I,J,xt,yt,ot,debug);
%

if nargin<8,
    debug = 0;
end;

if nargin<7,
    m = size(xt,2);
    n = size(yt,2);
    ex = ones(size(xt));
    ey = ones(size([yt;ot]));
    ex(3,:)=zeros(1,m);
    ex(:,1)=zeros(3,1);
    ex(2,2)=0;
    ej = [ex(:);ey(:)];
    options.cutjac = find(ej);
end;

m = size(xt,2);
n = size(yt,2);

for kkk = 1:30;
    %kkk
    
    %     xt = xtn;
    %     yt = ytn;
    %     ot = otn;
    
    
    [res,jac]=calcresandjac_tdoa_diffdim(D,I,J,xt,yt,ot,options);
    %[res,jac]=calcresandjac(D,I,J,xt,yt);
    %dz = -(jac\res);
    %dz = -(jac'*jac+eye(size(jac,2)))\(jac'*res);
    [u,s,v]=svd(full(jac),0);
    %keyboard;
    nrparam = size(jac,2);
    dof = nrparam-3;
    u = u(:,1:dof);
    s = s(1:dof,1:dof);
    v = v(:,1:dof);
    dz = -v*inv(s+0.00001*eye(size(s)))*u'*res;
    dzfull = zeros(3*m+4*n,1);
    dzfull(options.cutjac)=dz;
    [xtn,ytn,otn]=updatexy(xt,yt,ot,dzfull);
    [res2,jac2]=calcresandjac_tdoa_diffdim(D,I,J,xtn,ytn,otn,options);
    aa = [norm(res) norm(res+jac*dz) norm(res2)];
    bb = aa;
    bb=bb-bb(2);
    bb = bb/bb(1);
    cc = norm(jac*dz)/norm(res);
    % Check that the error actually gets smaller
    if norm(res)<norm(res2),
        % bad
        % check that the update is sufficiently big
        % otherwise it is maybe just numerical limitations
        if cc>1e-4,
            % OK then it is probably just the linearization that
            % is not good enough for this large step size, decrease
            kkkk = 1;
            while (kkkk<50) & (norm(res)<norm(res2)),
                dzfull = dzfull/2;
                [xtn,ytn,otn]=updatexy(xt,yt,ot,dzfull);
                [res2,jac2]=calcresandjac_tdoa_diffdim(D,I,J,xtn,ytn,otn,options);
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
        xt = xtn;
        yt = ytn;
        ot = otn;
    else
        %disp([num2str(kkk) '  stalled']);
    end
end;

xopt = xt;
yopt = yt;
oopt = ot;


function [res]=calcres_tdoa_diffdim(D,I,J,x,y,o)
nn = length(D);
m = size(x,2);
n = size(y,2);
V = x(:,I)-y(:,J);
Vt = V';
dd = sqrt(sum(V.^2,1))';
res = dd+o(J)'-D;

function [res,jac]=calcresandjac_tdoa_diffdim(D,I,J,x,y,o,options)
nn = length(D);
m = size(x,2);
n = size(y,2);
V = x(:,I)-y(:,J);
Vt = V';
dd = sqrt(sum(V.^2,1))';
%dd = sqrt(sum(V.^2,1))'+o(J)';
idd = 1./(dd);
res = dd+o(J)'-D;
%res = dd-D;
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
%keyboard;
jac = sparse([II;II;II;II;II;II;II],...
    [JJ1;JJ2;JJ3;JJ4;JJ5;JJ6;JJ7],...
    [VV1;VV2;VV3;VV4;VV5;VV6;VV7],nn,3*m+4*n);
jac = jac(:,options.cutjac);


function [xny,yny,ony]=updatexy(x,y,o,dz);

m = size(x,2);
n = size(y,2);
dz1 = dz(1:(3*m));
dz2 = dz((3*m+1):end);
xny = x + reshape(dz1,3,m);
yo = [y;o];
yony = yo + reshape(dz2,4,n);
yny = yony(1:3,:);
ony = yony(4,:);

