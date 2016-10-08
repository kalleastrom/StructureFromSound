function y=toa_trilateration_one_bundle(d,x,y0,inliers);
% y=toa_trilateration_one_bundle(d,x,y0)
% non-linear least squares optimization of y with y0
% as initial estimate, i e
% minimise
%  min_y  sum_i (d_i - sqrt(sum( (x(:,i)-y).^2 )))

if nargin<4,
    inliers = find(isfinite(d));
end;

x = x(:,inliers);
d = d(inliers,1);

[x_dim,m] = size(x);
d_dim = size(d,1);

if d_dim ~= m,
    error('toa_trilateration_one_bundle assumes size(d,1) = size(x,2)');
end;

if d_dim < x_dim,
    error('to few data');
else
    debug = 0;
    %keyboard;
    
    D = d;
    I = 1:size(D,1);
    J = ones(1,size(D,1));
    xt = x;
    yt = y0;
    
    for kkk = 1:10;
        %kkk
        [res,jac]=calcresandjac(D,I,J,xt,yt);
        %dz = -(jac\res);
        %dz = -(jac'*jac+eye(size(jac,2)))\(jac'*res);
        [u,s,v]=svd(full(jac),0);
        %keyboard;
        %     nrparam = size(jac,2);
        %     dof = nrparam-6;
        %     u = u(:,1:dof);
        %     s = s(1:dof,1:dof);
        %     v = v(:,1:dof);
        dz = -v*inv(s)*u'*res;
        [ytn]=updatey(yt,dz);
        [res2,jac2]=calcresandjac(D,I,J,xt,ytn);
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
                    dz = dz/2;
                    [ytn]=updatey(yt,dz);
                    [res2,jac2]=calcresandjac(D,I,J,xt,ytn);
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
            yt = ytn;
        else
            %disp([num2str(kkk) '  stalled']);
        end
    end;
    
    y = yt;
end;

function [res,jac]=calcresandjac(D,I,J,x,y);

%keyboard;
nn = length(D);
m = size(x,2);
n = size(y,2);
V = x(:,I)-y(:,J);
Vt = V';
dd = sqrt(sum(V.^2,1))';
idd = 1./dd;
res = dd-D;
II = (1:length(I))';
JJ1 = 1*ones(size(II));
JJ2 = 2*ones(size(II));
JJ3 = 3*ones(size(II));

VV1 = -idd.*Vt(:,1);
VV2 = -idd.*Vt(:,2);
VV3 = -idd.*Vt(:,3);

jac = sparse([II;II;II],[JJ1;JJ2;JJ3],[VV1;VV2;VV3],nn,3);


function [yny]=updatey(y,dz);

yny = y + dz;



