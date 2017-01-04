function [data,datap]=interp1d_with_derivative_restricted(f0,t0,w,sigma);
% a - assumed to be one row
% t0 - mid points of interwall
% w - half width of sample size (2w+1 samples) (-w, ... w)

% Cut out extra nr of samples on each side
extra = 20;
% Integer part of position
theshift = t0;
theshift_integerpart = round(theshift);
theshift_rest = theshift - theshift_integerpart;
% cut out 2w+1 samples and extra on each side
tmp1 = f0((theshift_integerpart-w-extra):(theshift_integerpart+w+extra));

% Convolution kernels (gaussian and derivative of gaussian)
x0 = [-10:10];
x0 = -30:30;
%sigma = 1.2;
miu = -theshift_rest; % or maybe -theshift_rest
h = (1/sqrt(2*pi)*sigma)*exp(-(x0 - miu).^2/(2*sigma^2));
h = h/sum(h);

% Convolve
tmp2 = conv(tmp1,h,'same');

% remove extra on edges
data = tmp2( (extra+1):(end-extra) );

% If second out argument is needed, calculate derivatives
if nargout>1,
    hprime = (-2*(x0-miu)/(2*sigma^2)).*h;
    tmp2prime = conv(tmp1,hprime,'same');
    datap = tmp2prime( (extra+1):(end-extra) );
end


