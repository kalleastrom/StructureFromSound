function [matning]=interp1d(f0,x0,a);

if 0,
    [n]=length(f0);
    NN=round(a*10);
    cutx=max(round(x0)-NN,1):min(round(x0)+NN,n);
    cuty=max(round(y0)-NN,1):min(round(y0)+NN,m);
    cutbild=bild(cuty,cutx);
    [x,y]=meshgrid(cutx,cuty);
    filter=exp(- ((x-x0).^2 + (y-y0).^2)/a^2)/(a^2*pi);
    matning=sum(sum(filter.*cutbild));
else
    [n]=length(f0);
    xx = 1:n;
    filter=exp(- ((xx-x0).^2 )/(2*a^2))/sqrt(2*a^2*pi);
    matning=sum(f0.*filter);
end;