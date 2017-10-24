function [z0, err, f0t, f1t] = find_translation_doppler_amplitude(f0, f1, thresh, a2, tt)
% finds the translation, doppler and constant amplitude between two similar 
% functions f0 and f1 (f1 wrt f0). The translation is estimated with an 
% accuracy given by thresh. The interpolation of the signals is done using 
% a Gaussian kernel of with std a2 and the initial translations search is 
% performed in the interval given by  tt=[minval,maxval]. The estimated 
% translation, doppler and amplitude can be found in z0 and the 
% corresponging error in err(end)

if isempty(thresh); thresh = 0.001; end
if isempty(a2); a2 = 2; end
if isempty(tt); tt = [-10 10]; end

if length(f0)~=length(f1)
    warning('The two signals are not of the same length')
end

tt = tt(1):tt(end);

x = 1:length(f0);
xmid = x(100+1:end-100);
f0t = interp1d(f0,xmid,a2);
for k = 1:length(tt)
    f1t = interp1d(f1,xmid-tt(k),a2);
    err(k) = norm(f0t-f1t);
end

[minv,mini] = min(err);
tau = tt(mini);
err = minv;

z0 =[tau;1;1]; % starting values of unknown paramters 
% translation and doppler

%% Use interpolation and derivative to find correct value

nbr_iter = 0;
J = zeros(length(xmid),3);
D = diag([1 1/1000 1]); % obs! what is D?
while (err(end)>thresh) && (nbr_iter<20)
%     tau0 = tau(end);
    f0t = interp1d(f0,xmid,a2);
    [f1t,f1td] = interp1d_with_derivative(f1,z0(2)*xmid+z0(1),a2); % obs! Where is the amplitude? 
    f1t = z0(3)*f1t; % add the amplitude
    f1td = z0(3)*f1td;
    J(:,1) = f1td'; % derivative wrt translation
    J(:,2) = (f1td.*xmid)'; % derivative wrt doppler
    J(:,3) = (f1t./z0(3))'; % derivative wrt amplitude
%     res = -(f0t'-f1t');
%     dz = -D*((J*D)\res);
    res = (f0t'-f1t');
    dz = D*((J*D)\res);
    if abs(dz(1))>1
        dz = dz*(1/abs(dz(1)));
    end
    znew = z0+dz;
    [f1tnew,~] = interp1d_with_derivative(f1,znew(2)*xmid+znew(1),a2);
    f1tnew = znew(3)*f1tnew;

    [norm(res) norm(res+J*dz) norm(f0t-f1tnew)];

    % if norm(resnew) > norm(res) d� minska steget.   
    z0 = znew;
    % obs! can we add this since we know that the amplitude is positive?
    % does it help?
    if z0(3)<0 
        z0(3) =1;
    end
    % find a new tau
    %tau(end+1) = tau0 - mean((f0t-f1t)./f1td);
    % Here one could use the least squares fit instead
    % This is similar in spirit to the Gauss-Newton iteration
    % f \approx = f1t + f1td*delta.
    % If we want to find delta which minimizes the sum of squared
    % residuals res = (f0t') - (f1t' + f1td'*delta)
    % This gives
    % (f0t'-f1t') = (-f1td')*delta
    % Think b = A*x
    % x = inv(A'*A)*(A'*b) or matlab x = A\b;
    %delta_old = - mean((f0t-f1t)./f1td);
    %delta_new = -(f1td'\(f0t'-f1t'));
    %[delta_old delta_new]
    %log10(abs(delta_new))
    %keyboard;
    %[norm((f0t-f1t)) norm((f0t') - (f1t' - f1td'*delta_old)) ]
    %tau(end+1) = tau0 + delta_new;
    %
    err(end+1) = norm(f0t-f1t);
    nbr_iter = nbr_iter+1;
end

trans = tau;

