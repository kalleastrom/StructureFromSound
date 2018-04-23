function [amp, err, f0t, f1t] = find_amp(f0, f1, a2)
% finds the amplitude difference gamma between two similar funtions f0 and 
% f1, s.t. f0(t)=gamma*f1(t). 

if isempty(a2); a2 = 2; end
if length(f0)~=length(f1)
    warning('The two signals are not of the same length')
end

x = 1:length(f0);
xmid = x(10:end-10);

if size(f0,1)>size(f0,2)
	f0 = f0';
    f1 = f1';
end

f0t = interp1d(f0,xmid,a2);
f1t = interp1d(f1,xmid,a2);
gamma0 = 1;


gamma = gamma0;
err = [];
% err = norm(f0t-f1t);

%% Use interpolation and derivative to find gamma

nbr_iter = 10;
for i = 1:nbr_iter
    gamma0 = gamma(end);
    [f1t,~] = interp1d_with_derivative(f1,xmid,a2);
    gamma_der = f1t;
    f1t = gamma0*f1t;
%     f1td = gamma0*f1td;

    if 0
        figure(10); clf; plot(xmid,f0t,'b'), hold on;
        plot(xmid,f1t,'r--'); hold off;
    end
        
    res = f0t'-f1t';
    err(end+1) = norm(res);
    
    dgamma = gamma_der'\res;
%     delta = -(f1td'\(f0t'-f1t'));
%     alpha(end+1) = alpha0 + delta;
    gamma(end+1) = gamma0+dgamma;
    nbr_iter = nbr_iter+1;
end

[f1t,~] = interp1d_with_derivative(f1,xmid,a2);
err(end+1) = norm(f0t-f1t);

amp = gamma(end);
err = err(end);
