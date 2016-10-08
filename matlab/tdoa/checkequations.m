function [d,dc,d2,d2c,d2r,d2rc]=checkequations(u,inl,o,A,x,y);
%
%

%
if nargin<5,
    D = 3;
    [uu,ss,vv]=svd(A);
    xr11 = (uu(:,1:D)*ss(1:D,1:D))';
    yr1 = vv(:,1:D)';
    xr33 = (uu(:,1:D))';
    yr3 = ss(1:D,1:D)*vv(:,1:D)';
    if (ss(1,1)/ss(3,3) > 10)
        xr11 = xr33;
        yr1  = yr3;
    end
    xtp = [zeros(D,1) xr11];
    y =[zeros(D,1) yr1];
    x = xtp/(-2);
end

[m,n]=size(u);
Rr = [1 -ones(1,n-1); zeros(n-1,1) eye(n-1)];
Rl = [1 -ones(1,m-1); zeros(m-1,1) eye(m-1)];
d = u-repmat(o,size(u,1),1);
d2 = d.^2;
d2r = Rl'*d2*Rr;
dc = reshape(sqrt( sum( (kron(ones(1,size(y,2)),x) - ...
    kron(y,ones(1,size(x,2))) ).^2 , 1 ) ),size(x,2),size(y,2));
uc = reshape(sqrt( sum( (kron(ones(1,size(y,2)),x) - ...
    kron(y,ones(1,size(x,2))) ).^2 , 1 ) ),size(x,2),size(y,2)) + repmat(o,size(x,2),1);
d2c = dc.^2
d2rc = Rl'*d2c*Rr;

figure(9);
subplot(4,3,1);
colormap(gray);
imagesc(d.*inl);
subplot(4,3,2);
colormap(gray);
imagesc(dc.*inl);
subplot(4,3,3);
colormap(gray);
imagesc((dc-d).*inl);
subplot(4,3,4);
colormap(gray);
imagesc(d2.*inl);
subplot(4,3,5);
colormap(gray);
imagesc(d2c.*inl);
subplot(4,3,6);
colormap(gray);
imagesc((d2c-d2).*inl);
subplot(4,3,7);
colormap(gray);
imagesc(d2r.*inl);
subplot(4,3,8);
colormap(gray);
imagesc(d2rc.*inl);
subplot(4,3,9);
colormap(gray);
imagesc((d2rc-d2r).*inl);
subplot(4,3,10);
colormap(gray);
imagesc(d2r(2:end,2:end).*inl(2:end,2:end));
subplot(4,3,11);
colormap(gray);
imagesc(d2rc(2:end,2:end).*inl(2:end,2:end));
subplot(4,3,12);
colormap(gray);
imagesc((d2rc(2:end,2:end)-d2r(2:end,2:end)).*inl(2:end,2:end));

d2r(:,1:7)
d2rc(:,1:7)
d2r(:,1:7)-d2rc(:,1:7)

