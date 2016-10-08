function [xopt,yopt,oopt,res,jac]=bundletdoa(D,I,J,xt,yt,ot,debug,xinorm);
% [xopt,yopt,oopt]=bundletdoa(D,I,J,xt,yt,ot,debug,xinorm);


if nargin<=6,
    debug = 0;
end;

if nargin==8,
    m = size(xt,2);
    n = size(yt,2);
    blubb = reshape(1:(3*m),3,m);
    dontmoveindex = [blubb(:,xinorm(1))'  blubb(2:3,xinorm(2))' blubb(3,xinorm(3))'];
    EE = speye(3*m+4*n);
    EE(:,dontmoveindex) = [];
end

for kkk = 1:10;
    kkk;
    [res,jac]=calcresandjac(D,I,J,xt,yt,ot);
    %dz = -(jac\res);
    %dz = -(jac'*jac+eye(size(jac,2)))\(jac'*res);
    %[u,s,v]=svd(full(jac),0);
    %keyboard;
    if 0,
        keyboard;
    elseif 0,
        jac0 = jac*EE;
        dz = -EE*((jac0'*jac0 + 0.01*eye(size(jac0,2)))\(jac0'*res));
        return
        dz = -EE*((jac0'*jac0 + 1*eye(size(jac0,2)))\(jac0'*res));
        %dz = EE*(-jac0\res);
    elseif 0,
        jac0 = jac*EE;
        %keyboard;
        dz = EE*(-jac0\res);
    elseif 0,
        [u,s,v]=svd(full(jac),'econ');
        u = u(:,1:(end-6));
        s = s(1:(end-6),1:(end-6));
        v = v(:,1:(end-6));
        dz = -v*inv(s)*u'*res;
    elseif 0,
        dz = -jac\res;
    else
        dz = -(jac'*jac + 10^(1)*speye(size(jac,2)))\(jac'*res);
    end
    [xtn,ytn,otn]=updatexy(xt,yt,ot,dz);
    [res2,jac2]=calcresandjac(D,I,J,xtn,ytn,otn);
    aa = [norm(res) norm(res+jac*dz) norm(res2)]
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
            kkkk = 1 ;
            while (kkkk<50) & (norm(res)<norm(res2)),
                dz = dz/2;
                [xtn,ytn,otn]=updatexy(xt,yt,ot,dz);
                [res2,jac2]=calcresandjac(D,I,J,xtn,ytn,otn);
                kkkk = kkkk+1;
            end
        end
    end
    if 1, %debug,
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
        disp([num2str(kkk) '  stalled']);
    end
end;

xopt = xt;
yopt = yt;
oopt = ot;


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


function [xny,yny,ony]=updatexy(x,y,o,dz);

m = size(x,2);
n = size(y,2);
dz1 = dz(1:(3*m));
dz2 = dz((3*m+1):end);
dz2 = reshape(dz2,4,n);
xny = x + reshape(dz1,3,m);
yny = y + dz2(1:3,:);
ony = o + dz2(4,:);



