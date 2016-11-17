% test 3 simply generates the function loaded by F0() and a translated copy
% of it. It is possible to choose how many decimals the translation should
% have. Then, noise can be added to thje signals (needs further
% development). Then a search for the translation in implemented. The
% number of valid numbers in the solution can be chosen-

% 
% % clear err
% 
% %% generate the function
% x = 1:2000;
% a = 0;
% f0 = F0(x,a);
% 
% % f0 = f0 + 0.001*randn(size(f0));
% 
% nbr_decimals = 2;
% true_translation = round(10^(nbr_decimals+1)*rand(1))/10^(nbr_decimals);
% 
% f1 = F0(x+true_translation,a);
% % f1 = f1+ 0.001*randn(size(f1));
% 
% % plot the functions
% figure(1); clf; 
% subplot(2,1,1); plot(x,f0,'b');
% subplot(2,1,2); plot(x,f1,'r');
% 
% %% Discrete search for initial guess
% 
% thresh = 10^(-nbr_decimals);
% a = 10;
% tt = [-10 10];
% 
% [trans, err] = find_translation(f0,f1,thresh,a,tt);
% 
% % a = 2;
% % 
% % tt = -10:10; 
% % xmid = x(100+1:end-100);
% % f0t = interp1d(f0,xmid,a);
% % for k = 1:length(tt)
% %     f1t = interp1d(f1,xmid-tt(k),a);
% %     err(k) = norm(f0t-f1t);
% % end
% % 
% % [minv,mini] = min(err);
% % tau = tt(mini);
% % err = minv;
% % 
% % %% Use interpolation and derivative to find correct value
% % 
% % threshold = 10^(-nbr_decimals);
% % a = 10;
% % 
% % while err(end)>threshold
% %     tau0 = tau(end);
% %     f0t = interp1d(f0,xmid,a);
% %     [f1t,f1td] = interp1d_with_derivative(f1,xmid-tau0,a);
% %     % find a new tau
% %     tau(end+1) = tau0 - mean((f0t-f1t)./f1td);
% %     err(end+1) = norm(f0t-f1t);
% % end

%% try finding the translation between more than two channels (variable)

nbr_channels = 8;
channels = cell(nbr_channels,1);

x = 1:2000;
a = 0;
channels{1} = F0(x,a);
% channels{1} = channels{1} + 0.001*randn(size(channels{1}));

nbr_decimals = 2; % how many decimals the translation should have
true_translation = round(2*10^(nbr_decimals+1)*rand(1,nbr_channels-1))/10^(nbr_decimals) - 10;
true_translation = [0 true_translation];
true_translation = zeros(size(true_translation));
true_doppler = 1 + 0.01*randn(1,nbr_channels-1);
true_doppler = [1 true_doppler];

% create the channels and plot them
figure(1); clf; 
subplot(nbr_channels,1,1); plot(x,channels{1});
title(['The ' num2str(nbr_channels) ' different channels']);
for i = 2:nbr_channels
    channels{i} = F0((1/true_doppler(i))*(x-true_translation(i)),a);
    channels{i} = channels{i}+0*randn(size(channels{i}));
    subplot(nbr_channels,1,i); plot(x,channels{i});
end

thresh = 10^(-8); % decides how good the translation estimation needs to be
a = 10;
tt = [-15 15]; % the translations to be tried
trans = zeros(nbr_channels,1);
err = zeros(nbr_channels,1);
for i = 2:nbr_channels    
    [curr_z, curr_err] = find_translation_and_doppler(channels{1}, channels{i}, thresh, a, tt);
    z(:,i) = curr_z;
    err(i) = curr_err(end);
end

% find the "original" signal be taking the mean of all signals, translated

xmid = x(100+1:end-100); % make sure that the cut part is larger than any translation
channels_t = cell(nbr_channels,1);
a = 2; % what should a be here?
channels_mean = zeros(1,length(xmid));
figure(2); clf; hold on;
figure(3); clf; hold on;
for i = 1:nbr_channels
    channels_t{i} = interp1d(channels{i},xmid-trans(i),2);
%     subplot(nbr_channels,1,i); plot(x,channels{i});
    figure(2); plot(xmid,channels{i}(xmid),'b');
    figure(3); plot(xmid,channels_t{i},'b');
    % compute the mean signal
    channels_mean = channels_mean+channels_t{i}./nbr_channels;
end

% the mean will be a bit smoothed, while the original signals in Figure 2
% will not
figure(2); plot(xmid,channels_mean,'r--'); hold off; 
title('The original signals and the computed mean')
figure(3); plot(xmid,channels_mean,'r--'); hold off;
title('The translated signals and the computed mean')

[trans true_translation']

