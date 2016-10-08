function y=toa_trilateration(d,x,y0,index,inliers);
% y=toa_trilateration_one_bundle(d,x,y0,index,inliers)
%

[m,n]=size(d);
[x_dim,xm]=size(x);

if nargin<3,
    y0 = NaN*ones(x_dim,n);
end

if nargin<4,
    index = find(ones(1,n));
end;

if nargin<5,
    inliers = isfinite(d);
end;

y = y0;

for jj = index,
    yy = toa_trilateration_one_ransac(d(:,jj),x,10,Inf,find(inliers(:,jj)));
    yy = real(yy);
    yy = toa_trilateration_one_bundle(d(:,jj),x,yy,inliers(:,jj));
    y(:,jj) = yy;
end;