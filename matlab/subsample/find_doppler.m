function [dopp, err, f0t, f1t] = find_doppler(f0, f1, a2, thresh)
% finds the doppler shift alpha between two similar function f0 and f1, 
% s.t. f0(t)=f1(alpha*t). The interpolation of the signals is done using a
% Gaussian kernel of width a2 and the initial search starts at alpha=1.
% thresh is a threshold for the error.

if isempty(a2); a2 = 2; end
if length(f0)~=length(f1)
    warning('The two signals are not of the same length')
end

x = 1:length(f0);
xmid = x(100:end-100);

f0t = interp1d(f0,xmid,a2);
f1t = interp1d(f1,xmid,a2);
alpha0 = 1;

alpha = alpha0;
err = norm(f0t-f1t);

%% Use interpolation and derivative to find alpha

nbr_iter = 0;
while (err(end)>thresh) && (nbr_iter<10)
    alpha0 = alpha(end);
    [f1t,f1td] = interp1d_with_derivative(f1,alpha0*xmid,a2);
    delta = -(f1td'\(f0t'-f1t'));
%     delta = - mean((f0t-f1t)./f1td);
%     alpha(end+1) = alpha0 + 0.1*delta;
    alpha(end+1) = alpha0 + delta;
    err(end+1) = norm(f0t-f1t);
    nbr_iter = nbr_iter+1;
end

dopp = alpha(end);
err = err(end);
