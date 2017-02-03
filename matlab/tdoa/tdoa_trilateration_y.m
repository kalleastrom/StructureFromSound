function [y,o,inl]=tdoa_trilateration_y(u,x,y0,o0,index,inliers);
% y=tdoa_trilateration_one_bundle(u,x,y0,index,inliers)
%

[m,n]=size(u);
[x_dim,xm]=size(x);

if nargin<3,
    y0 = NaN*ones(x_dim,n);
end

if nargin<4,
    o0 = NaN*ones(1,n);
end

if nargin<5,
    index = find(ones(1,n));
end;

if nargin<6,
    inliers = isfinite(u);
end;

inl = zeros(size(inliers));

ransac_tol = 0.1;

%keyboard;

kkk = 0;
for jj = index,
    format compact;
    kkk = kkk+1;
    if round(kkk/10)== (kkk/10),
        [kkk length(index) jj]
    end
    %keyboard;
    [yy,oo,inl_one_y,nr_inlierid,err_rms] = tdoa_trilateration_y_one_ransac(u(:,jj),x,30,ransac_tol,find(inliers(:,jj)));
    yy = real(yy);
    oo = real(oo);
%     [yy,oo] = tdoa_trilateration_y_one_bundle(u(:,jj),x,yy,oo,find(inliers(:,jj))); % BUGG. Man ska inte bundla 
%                                                                                     % över de ursprungliga inliers utan de nya
    %keyboard;
    if length(inl_one_y)>3,
    [yy,oo] = tdoa_trilateration_y_one_bundle(u(:,jj),x,yy,oo,inl_one_y); % BUGG. Man ska inte bundla 
    end;                                                                              % över de ursprungliga inliers utan de nya
    y(:,jj) = yy;
    o(:,jj) = oo;
    inl(inl_one_y,jj)=ones(length(inl_one_y),1);
end;

