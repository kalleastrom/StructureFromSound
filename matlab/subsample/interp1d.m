function [matning]=interp1d(f0,x0,a);
% Performs interpolation in 1d using a Gaussian kernel (approximation for
% sinc followed by Gaussian)
% f0 is the function to be interpolated, x0 is the value that
% the interpolated values should be found for and a in the standard
% deviation of the Gaussian kernel that the interpolation is performed
% with. 




if 0,
    [n]=length(f0);
    NN=round(a*10);
    cutx=max(round(x0)-NN,1):min(round(x0)+NN,n);
    cuty=max(round(y0)-NN,1):min(round(y0)+NN,m);
    cutbild=bild(cuty,cutx);
    [x,y]=meshgrid(cutx,cuty);
    filter=exp(- ((x-x0).^2 + (y-y0).^2)/a^2)/(a^2*pi);
    matning=sum(sum(filter.*cutbild));
elseif 1,
    % try to make it work for whole matrices, so far no cutting
    turn = 0;
    if size(f0,1)<size(f0,2)
        f0 = f0';
        turn = 1;
    end
    [n] = length(f0);
    xx = 1:n;
    [Ts,T] = ndgrid(x0,xx);
    matning = exp(- ((T-Ts).^2 )./(2*a^2))./sqrt(2*a^2*pi)*f0; % test this
    if turn
        matning = matning';
    end
else
    [n]=length(f0);
    xx = 1:n;
    filter=exp(- ((xx-x0).^2 )/(2*a^2))/sqrt(2*a^2*pi);
    matning=sum(f0.*filter);
end;