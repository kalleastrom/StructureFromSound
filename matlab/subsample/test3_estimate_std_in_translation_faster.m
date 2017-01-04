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

nbr_channels = 2;
channels = cell(nbr_channels,1);
noise_std = 0.003;

x = 1:2000;
a1 = 0; % Blurr 1
a2 = 2; % Blurr 2
channels{1} = F0(x,a1);
channels{1} = channels{1} + noise_std*randn(size(channels{1}));

nbr_decimals = 2; % how many decimals the translation should have
true_translation = round(2*10^(nbr_decimals+1)*rand(1,nbr_channels-1))/10^(nbr_decimals) - 10;
true_translation = [0 true_translation];
%true_translation = [0 0];

for kk = 1:2;
    
    % create the channels and plot them
    figure(1); clf;
    subplot(nbr_channels,1,1); plot(x,channels{1});
    title(['The ' num2str(nbr_channels) ' different channels']);
    for i = 2:nbr_channels
        channels{i} = F0(x+true_translation(i),a1);
        channels{i} = channels{i}+noise_std*randn(size(channels{i}));
        subplot(nbr_channels,1,i); plot(x,channels{i});
    end
    
    thresh = 10^(-8); % decides how good the translation estimation needs to be
    tt = [-15 15]; % the translations to be tried
    trans = zeros(nbr_channels,1);
    err = zeros(nbr_channels,1);
    for i = 2:nbr_channels
        [curr_trans, curr_err] = find_translation(channels{1}, channels{i}, thresh, a, tt);
        trans(i) = curr_trans(end);
        err(i) = curr_err(end);
    end
    
    % find the "original" signal be taking the mean of all signals, translated
    
    xmid = x(100+1:end-100); % make sure that the cut part is larger than any translation
    channels_t = cell(nbr_channels,1);
    channels_mean = zeros(1,length(xmid));
    figure(2); clf; hold on;
    figure(3); clf; hold on;
    for i = 1:nbr_channels
        channels_t{i} = interp1d(channels{i},xmid-trans(i),a);
        %channels_t{i} = interp1d(channels{i},xmid-true_translation(i),a);
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
    
    [trans true_translation'];
    all_trans(kk)=trans(2);
    figure(10+kk); plot(channels_t{2}-channels_t{1},'.')
    noise_var_estimate = var(channels_t{2}-channels_t{1}) / (2*(1/sqrt(2*pi*a22^2)))
    noise_std_estimate = sqrt(noise_var_estimate)

end

mean(all_trans)
true_translation(2)
std(all_trans)
% formel för vad vi tror std ska bli

%% Estimate error, we first need \int (F')^2 dx

xx = -15:15;
gg = (-2*xx/(2*a^2)).*(1/sqrt(2*pi*a^2)).*exp( - (xx.^2)/(2*a^2) );
tmp = conv2(channels{1},gg,'same');
% Gör det på medelvärdesbildningen
tmp = tmp(100:1900);
EA = 2*sum(tmp.^2);
% re-estreck = 2*noise_std^2 * diskret deltafunction.
% Re-esteck = 2*noise_std^2 * normalfördelad med std a*sqrt(2) ????
a22 = a2*sqrt(2);
xx = -15:15;
geestreck = 2*noise_std^2 * (1/sqrt(2*pi*a22^2)).*exp( - (xx.^2)/(2*a22^2) );
gg = (-2*xx/(2*a2^2)).*(1/sqrt(2*pi*a2^2)).*exp( - (xx.^2)/(2*a2^2) );
tmp = conv2(channels{1},gg,'same');
tmp2 = conv2(tmp,geestreck,'same');
tmp = tmp(100:1900);
tmp2 = tmp2(100:1900);
Vb = 4*sum(tmp.*tmp2);
VX = Vb/(EA^2);
SX = sqrt(VX)
Vb_approx = noise_std^2*4*sum(tmp.*tmp);
VX_approx = noise_std^2 / (sum(tmp.*tmp));
noise_var_estimate = var(channels_t{2}-channels_t{1}) / (2*(1/sqrt(2*pi*a22^2)))
noise_std_estimate = sqrt(noise_var_estimate)





