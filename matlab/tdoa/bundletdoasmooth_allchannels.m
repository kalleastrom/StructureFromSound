function [asolopt,res,jac]=bundletdoasmooth_allchannels(data,asol,model,debug,xinorm);
% [xopt,yopt,oopt]=bundletdoa(D,I,J,xt,yt,ot,debug,xinorm);

%keyboard;
if nargin<3,
    model.smeas = 1;
    model.smotion = 1;
end;

if nargin<4,
    debug = 0;
end;

if nargin<5,
    xinorm = [1 2 3];
end

m = size(asol.x,2);
n = size(asol.y,2);
blubb = reshape(1:(3*m),3,m);
dontmoveindex = [blubb(:,xinorm(1))'  blubb(2:3,xinorm(2))' blubb(3,xinorm(3))'];
EE = speye(3*m+3*n);
EE(:,dontmoveindex) = [];
asolt = asol;

for kkk = 1:10;
    kkk;
    [res,jac]=calcresandjac(data,asolt,model);
    %dz = -(jac\res);
    %dz = -(jac'*jac+eye(size(jac,2)))\(jac'*res);
    %[u,s,v]=svd(full(jac),0);
    if 0,
        keyboard;
        %%
        for ki = 1:100,
        litet = 0.0001;
        asol0 = asolt;
        dz = zeros(size(jac,2),1);
        dz(ki)=1;
        [asol1]=updatexy(asol0,litet*dz);
        [res0,jac0]=calcresandjac(data,asol0,model);
        [res1,jac1]=calcresandjac(data,asol1,model);
        [(res1-res0)/litet jac0*dz];
        d1 = (res1-res0)/litet;
        d2 = jac0*dz;
        ki
        [norm(d1) norm(d2) norm(d1-d2)]
        pause
        end
        %%
    end
    if 0,
        keyboard;
    elseif 0,
        jac0 = jac*EE;
        dz = -EE*((jac0'*jac0 + 0.01*eye(size(jac0,2)))\(jac0'*res));
        return
        dz = -EE*((jac0'*jac0 + 1*eye(size(jac0,2)))\(jac0'*res));
        %dz = EE*(-jac0\res);
    elseif 1,
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
        dz = -(jac'*jac + 10^(-4)*eye(size(jac,2)))\(jac'*res);
    end
    [asoltn]=updatexy(asolt,dz);
    [res2,jac2]=calcresandjac(data,asoltn,model);
    aa = [norm(res) norm(res+jac*dz) norm(res2)];
    aa0 = aa;
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
                [asoltn]=updatexy(asolt,dz);
                [res2,jac2]=calcresandjac(data,asoltn,model);
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
        aa
        bb
        cc
        %keyboard;
    end;
    if norm(res2)<norm(res)
        asolt = asoltn;
    else
        disp([num2str(kkk) '  stalled']);
    end
end;

asolopt = asolt;

function [res,jac]=calcresandjac(data,asol,model);

[res1,jac1]=calcresandjac1(data,asol);
[res2,jac2]=calcresandjac2(data,asol);
res = [(1/model.smeas)*res1;(1/model.smotion)*res2];
jac = [(1/model.smeas)*jac1;(1/model.smotion)*jac2];


function [res,jac]=calcresandjac1(data,asol);

nn = size(data,1);
m = size(asol.x,2);
n = size(asol.y,2);
V1 = asol.x(:,data(:,1))-asol.y(:,data(:,3));
V1t = V1';
V2 = asol.x(:,data(:,2))-asol.y(:,data(:,3));
V2t = V2';
dd1 = sqrt(sum(V1.^2,1))';
idd1 = 1./dd1;
dd2 = sqrt(sum(V2.^2,1))';
idd2 = 1./dd2;
res = dd2-dd1-data(:,4);
II = (1:nn)';
JJ1 = (data(:,1)-1)*3+1;
JJ2 = (data(:,1)-1)*3+2;
JJ3 = (data(:,1)-1)*3+3;
JJ4 = (data(:,2)-1)*3+1;
JJ5 = (data(:,2)-1)*3+2;
JJ6 = (data(:,2)-1)*3+3;
JJ7 = (data(:,3)-1)*3+1+3*m;
JJ8 = (data(:,3)-1)*3+2+3*m;
JJ9 = (data(:,3)-1)*3+3+3*m;

VV1 = -idd1.*V1t(:,1);
VV2 = -idd1.*V1t(:,2);
VV3 = -idd1.*V1t(:,3);
VV4 = idd2.*V2t(:,1);
VV5 = idd2.*V2t(:,2);
VV6 = idd2.*V2t(:,3);
VV7 = -idd2.*V2t(:,1)+idd1.*V1t(:,1);
VV8 = -idd2.*V2t(:,2)+idd1.*V1t(:,2);
VV9 = -idd2.*V2t(:,3)+idd1.*V1t(:,3);

jac = sparse([II;II;II;II;II;II;II;II;II],[JJ1;JJ2;JJ3;JJ4;JJ5;JJ6;JJ7;JJ8;JJ9],[VV1;VV2;VV3;VV4;VV5;VV6;VV7;VV8;VV9],nn,3*m+3*n);

function [res,jac]=calcresandjac2(data,asol);

NN = size(asol.y,2);
j1= (1:(NN-2))';
j2 = j1+1;
j3 = j2+1;
%j1 = cc(:,1);
%j2 = cc(:,2);
%j3 = cc(:,3);

res1 = asol.y(1,j1) + asol.y(1,j3) - 2*asol.y(1,j2);
res2 = asol.y(2,j1) + asol.y(2,j3) - 2*asol.y(2,j2);
res3 = asol.y(3,j1) + asol.y(3,j3) - 2*asol.y(3,j2);
nn = length(j1);
m = size(asol.x,2);
n = size(asol.y,2);
res = [res1';res2';res3'];
IIx = (1:(nn))';
IIy = ((nn+1):(2*nn))';
IIz = ((2*nn+1):(3*nn))';
JJ1x = (j1-1)*3+1+3*m;
JJ1y = (j1-1)*3+2+3*m;
JJ1z = (j1-1)*3+3+3*m;
JJ2x = (j2-1)*3+1+3*m;
JJ2y = (j2-1)*3+2+3*m;
JJ2z = (j2-1)*3+3+3*m;
JJ3x = (j3-1)*3+1+3*m;
JJ3y = (j3-1)*3+2+3*m;
JJ3z = (j3-1)*3+3+3*m;

VV  = ones(size(JJ1x));

jac = sparse([IIx;IIx;IIx;IIy;IIy;IIy;IIz;IIz;IIz], ...
    [JJ1x;JJ2x;JJ3x;JJ1y;JJ2y;JJ3y;JJ1z;JJ2z;JJ3z], ...
    [VV;-2*VV;VV;VV;-2*VV;VV;VV;-2*VV;VV],3*nn,3*m+3*n);
res = res;
jac = jac;

function [asolny]=updatexy(asol,dz);

asolny = asol;
m = size(asol.x,2);
n = size(asol.y,2);
dz1 = dz(1:(3*m));
dz2 = dz((3*m+1):end);
dz2 = reshape(dz2,3,n);
asolny.x = asolny.x + reshape(dz1,3,m);
asolny.y = asolny.y + dz2;



