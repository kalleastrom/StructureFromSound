function [y,o,inliers_ut,nr_inliers_ut,err_rms]=tdoa_trilateration_y_one_ransac(d,x,ransac_k,ransac_tol,inlierid);
% [y,o,inliers_ut,nr_inliers_ut,err_rms]=tdoa_trilateration_y_one_ransac(d,x,ransac_k,ransac_tol,inlierid);
% trilateration of one point using ransac
% 
%

if nargin<3,
    ransac_k = 10;
end

if nargin<4,
    ransac_tol = 0.01;
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

if md==(xdim+1),
    % Minimal problem
    error('Minimal problem - to little data to perform ransac');
elseif md<(xdim+1),
    error('Underconstrained problem - to little data to perform ransac');
else
    % Main part (more data (md) than dimensionality of problem xdim
    
    counter = 0;
    ys = zeros(xdim,2*ransac_k);
    os = zeros(1,2*ransac_k);
    evals = zeros(1,2*ransac_k);
    nr_inliers = zeros(1,2*ransac_k);
    
    for kk = 1:ransac_k
        rp = randperm(md);
        rp = rp(1:(xdim+1));
        %keyboard;
        [ysols,osols]=tdoa_trilateration_y_one_point(d(rp),x(:,rp));
        for ii = 1:2;
            counter = counter+1;
            yy = ysols(:,ii);
            oo = osols(:,ii);
            % Make solutions real or check that they are real???
            yy = real(yy);
            oo = real(oo);
            d_proj = tdoa_calc_u_from_xyo(x,yy,oo);
            ys(:,counter)=yy;
            os(:,counter)=oo;
            inlierid2 = find(abs(d_proj-d) < ransac_tol);
            %keyboard;
            nr_inliers(1,counter)=length(inlierid2);
            evals(1,counter)= sqrt(sum( (d_proj(inlierid2)-d(inlierid2)).^2 ));
        end
    end
    
    % Select the solution with (1) the most inliers and
    % among these the one with the lowest error score.
    [maxv,maxi]=max(nr_inliers);
    ok1 = find(nr_inliers==maxv);
    [minv,mini]=min(evals(ok1));
    besti = ok1(mini);
    
    y = ys(:,besti);
    o = os(:,besti);
    %keyboard;
    d_proj = tdoa_calc_u_from_xyo(x,y,o);
    inlierid2 = find((abs(d_proj-d) < ransac_tol));
    inliers_ut = inlierid(inlierid2);
    nr_inliers_ut=length(inlierid2);
    err_rms = evals(besti);
%     maxv
%     inlierid2
%     [d d_proj]
%     abs(d-d_proj)
%     y
%     o
    %keyboard
end

