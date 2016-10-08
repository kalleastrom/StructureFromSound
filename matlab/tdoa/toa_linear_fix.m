function [sol,stats,bb,bb0] = toa_linear_fix_2D(d,A)
%[sols,stats] = toa_linear(d)
% Solver for overdetermined case TOA node calibration
%  at least  senders and 5 receivers
%  measurments according to
%    d(i,j) = sqrt(sum( (x(:,i)-y(:,j)).^2 ))
%  Goal of solver is to reconstruct x and y from d.
%  There are 42 solutions counted with multiplicity and
%  complex solutions.
% In:
%    d - 5x5 matrix with distances
%    settings - settings for solver use
%               load toa_3D_55_settings
%    xt -
%    yt -
% Out:
%    sols  - cell array with all real and feasible solvers
%    stats - statistics from solver

%keyboard;
[m,n]=size(d);
if n>m,
    d = d';
    A = A';
    [m,n]=size(d);
    flip = 1;
else
    flip = 0;
end

D = 3;

[uu,ss,vv]=svd(A);

xr11 = (uu(:,1:D)*ss(1:D,1:D))';
yr1 = vv(:,1:D)';

xr33 = (uu(:,1:D))';
yr3 = ss(1:D,1:D)*vv(:,1:D)';

if 1, %(ss(1,1)/ss(3,3) > 10)
    xr11 = xr33;
    yr1  = yr3;
end

xtp = [zeros(D,1) xr11];
yt =[zeros(D,1) yr1];
xt = xtp/(-2);


%% construct linear constraints
nrofunknowns = (D*(D+1))/2 + D;
E = eye(nrofunknowns);
for i = 1:nrofunknowns;
    xv(i) = multipol(1,E(i,:)');
end
one  = multipol(1,zeros(nrofunknowns,1));
zero = multipol(0,zeros(nrofunknowns,1));
if D==3,
    Cv = [xv(1) xv(2) xv(3); xv(2) xv(4) xv(5); xv(3) xv(5) xv(6)];
elseif D==2,
    Cv = [xv(1) xv(2); xv(2) xv(3)];
end
% Lv = [xv(1) 0 0 ; xv(2) xv(3) 0; xv(4) xv(5) xv(6)];
if D==3,
    bv = [xv(7) xv(8) xv(9)]';
elseif D==2,
    bv = [xv(4) xv(5)]';
end


%% here comes the linear equations
for i = 2:size(xt,2);
    eqs(i-1) = (-2*xt(:,i)'*bv + xt(:,i)'*Cv*xt(:,i)) - (d(i,1)^2-d(1,1)^2);
end
[cfm_linear0,monlinear]=polynomials2matrix(eqs);

AA = cfm_linear0(:,1:(end-1));
bb = cfm_linear0(:,end);
zz0 = -pinv(AA)*bb;

bb0 = -AA*zz0;

if 1,
    figure(7); clf; hold off;
    plot([bb bb0]);
    figure(8);
    plot([A]);
    
end

%keyboard;

if 0,
    %    nnn = size(AA,1);
    for kk =1:5;
        zz0 = -(AA'*AA)\(AA'*bb);
        bb0 = -AA*zz0;
        ww = abs(bb-bb0);
        W = diag(ones(size(ww))./ww);
        
        
        zz0 = -(AA'*W*AA)\(AA'*W*bb);
        bb0 = -AA*zz0;
        ww = abs(bb-bb0);
        W = diag(ones(size(ww))./ww);
        plot(ww);
        
        [blubb,badj]=max(ww);
        
        AA(badj,:)=[];
        bb(badj)=[];
        
    end;
end

AA = cfm_linear0(:,1:(end-1));
bb = cfm_linear0(:,end);
bb0 = -AA*zz0;


C = evaluate(Cv,zz0);
b = evaluate(bv,zz0);

try
    L    = chol(inv(C));
    ok   = 1;
catch % if not positive-definite, throw away
    %         disp(' C is not positive definite')
    L    = 0;
    ok   = 0;
    %keyboard;
    % Sometimes C is not positive definite. Here comes a hack
    %
    tt = min(eig(C));
    C = C + (-tt+1)*eye(size(C));
    L    = chol(inv(C));
    ok   = 0;
end
sol.C = C;
sol.b = b;
sol.L = L;
stats.ok = ok;
x = inv(L')*(xt);
y = L*(yt+repmat(b,1,n));
if flip,
    sol.x = y;
    sol.y = x;
    %sol.dd = dd';
else
    sol.x = x;
    sol.y = y;
    %sol.dd = dd;
end
stats.ss = ss;

