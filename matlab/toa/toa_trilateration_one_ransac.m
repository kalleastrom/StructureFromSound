function [y,inliers,nr_inliers,err_rms]=toa_trilateration_one_ransac(d,x,ransac_k,ransac_tol,inlierid);
% [y,inliers,nr_inliers,err_rms]=toa_trilateration_one_ransac(d,x,ransac_k,ransac_tol,inlierid);
% trilateration of one point using ransac
% 
%

if nargin<3,
    ransac_k = 10;
end

if nargin<4,
    ransac_tol = Inf;
end

if nargin<5,
    inlierid = find(isfinite(d));
end

d = d(inlierid);
x = x(:,inlierid);

[xdim,m]=size(x);
md = size(d,1);

if md~=m,
    error('wrong dimension');
end

if md==xdim,
    % Minimal problem
    error('Minimal problem - to little data to perform ransac');
elseif md<xdim,
    error('Underconstrained problem - to little data to perform ransac');
else
    % Main part (more data (md) than dimensionality of problem xdim
    
    counter = 0;
    ys = zeros(xdim,2*ransac_k);
    evals = zeros(1,2*ransac_k);
    nr_inliers = zeros(1,2*ransac_k);
    
    for kk = 1:ransac_k
        rp = randperm(md);
        rp = rp(1:xdim);
        ysols=toa_trilateration_one_point(d(rp),x(:,rp));
        for ii = 1:2;
            counter = counter+1;
            yy = ysols(:,ii);
            d_proj = toa_calc_d_from_xy(x,yy);
            ys(:,counter)=yy;
            inlierid = find(abs(d_proj-d) < ransac_tol);
            nr_inliers(1,counter)=length(inlierid);
            evals(1,counter)= sqrt(sum( (d_proj(inlierid)-d(inlierid)).^2 ));
        end
    end
    
    % Select the solution with (1) the most inliers and
    % among these the one with the lowest error score.
    [maxv,maxi]=max(nr_inliers);
    ok1 = find(nr_inliers==maxv);
    [minv,mini]=min(evals(ok1));
    besti = ok1(mini);
    
    y = ys(:,besti);
    d_proj = toa_calc_d_from_xy(x,y);
    inliers = (abs(d_proj-d) < ransac_tol);
    nr_inliers=sum(inliers);
    err_rms = evals(besti);
    
    
    
end

