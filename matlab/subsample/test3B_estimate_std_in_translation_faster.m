% Simulates a signal and a number of translated copies. Also adds Gaussian 
% distributed noise to all channels. Computes the subsample translation for
% each channel with a given number of decimals. Also estimates the std of 
% the translation. Similar to test3_estimate_std_in_translation_faster.m
% but more organised. Prints the translation, translation std and noise std
% (true and estimated) at the end. 2017-08-14


% test 3B simply generates the function loaded by F0() and a translated copy
% of it. It is possible to choose how many decimals the translation should
% have. Then, noise can be added to thje signals (needs further
% development). Then a search for the translation in implemented. The
% number of valid numbers in the solution can be chosen-


%% First generate signals without noise

% settings
nbr_channels = 2;
channels = cell(nbr_channels,1);
noise_std = 0.03;
x = 1:1000;
a1 = 0; % Blurr 1
a2 = 1.5; % Blurr 2 
a22 = a2*sqrt(2);
N = 100; % Nr of runs with different noise realizations

% generate ground truth translations
nbr_decimals = 2; % how many decimals the translation should have
true_translation = round(2*10^(nbr_decimals+1)*rand(1,nbr_channels-1))/10^(nbr_decimals) - 10;
true_translation = [0 true_translation];

%% Run N times

all_trans = zeros(1,N);

tic;
for kk = 1:N;
    
    % Generate data and add noise
    for i = 1:nbr_channels
        channels{i} = F0(x+true_translation(i),a1);
        channels{i} = channels{i}+noise_std*randn(size(channels{i}));
    end
    
    % Estimate the translation
    thresh = 10^(-8); % decides how good the translation estimation needs to be
    tt = [-15 15]; % the translations to be tried
    trans = zeros(nbr_channels,1);
    err = zeros(nbr_channels,1);
    for i = 2:nbr_channels
        if 0,
            % Version 1
            [curr_trans, curr_err] = find_translation(channels{1}, channels{i}, thresh, a2, tt);
        else
            % Version 2
            [curr_trans, curr_err] = find_translation2(channels{1}, channels{i}, thresh, a2, tt);
        end
        trans(i) = curr_trans(end);
        err(i) = curr_err(end);
    end
    
% Skip this for speed for the moment
%     % find the "original" signal be taking the mean of all signals, translated
%     
    xmid = x(100+1:end-100); % make sure that the cut part is larger than any translation
    channels_t = cell(nbr_channels,1);
    channels_mean = zeros(1,length(xmid));
%     figure(2); clf; hold on;
%     figure(3); clf; hold on;
    for i = 1:nbr_channels
        channels_t{i} = interp1d(channels{i},xmid-trans(i),a2);
        %channels_t{i} = interp1d(channels{i},xmid-true_translation(i),a);
        %     subplot(nbr_channels,1,i); plot(x,channels{i});
        %figure(2); plot(xmid,channels{i}(xmid),'b');
        %figure(3); plot(xmid,channels_t{i},'b');
        % compute the mean signal
        channels_mean = channels_mean+channels_t{i}./nbr_channels;
    end
    
    % the mean will be a bit smoothed, while the original signals in Figure 2
    % will not
%     figure(2); plot(xmid,channels_mean,'r--'); hold off;
%     title('The original signals and the computed mean')
%     figure(3); plot(xmid,channels_mean,'r--'); hold off;
%     title('The translated signals and the computed mean')
    
    [trans true_translation'];
    all_trans(kk)=trans(2);
    %figure(10+kk); plot(channels_t{2}-channels_t{1},'.')
    noise_var_estimate = var(channels_t{2}-channels_t{1}) / (2*(1/sqrt(2*pi*a22^2)));
    noise_std_estimate = sqrt(noise_var_estimate);
    all_noise(kk)=noise_std_estimate;
end
toc

% formel f�r vad vi tror std ska bli

%% Estimate error, we first need \int (F')^2 dx

xx = -15:15;
gg = (-2*xx/(2*a2^2)).*(1/sqrt(2*pi*a2^2)).*exp( - (xx.^2)/(2*a2^2) );
tmp = conv2(channels{1},gg,'same');
% G�r det p� medelv�rdesbildningen
tmp = tmp(100+1:end-100);
EA = 2*sum(tmp.^2);
% re-estreck = 2*noise_std^2 * diskret deltafunction.
% Re-esteck = 2*noise_std^2 * normalf�rdelad med std a*sqrt(2) ????
a22 = a2*sqrt(2);
xx = -15:15;
geestreck = 2*noise_std^2 * (1/sqrt(2*pi*a22^2)).*exp( - (xx.^2)/(2*a22^2) );
gg = (-2*xx/(2*a2^2)).*(1/sqrt(2*pi*a2^2)).*exp( - (xx.^2)/(2*a2^2) );
tmp = conv2(channels{1},gg,'same');
tmp2 = conv2(tmp,geestreck,'same');
tmp = tmp(100+1:end-100);
tmp2 = tmp2(100+1:end-100);
Vb = 4*sum(tmp.*tmp2);
VX = Vb/(EA^2);
SX = sqrt(VX)
Vb_approx = noise_std^2*4*sum(tmp.*tmp);
VX_approx = noise_std^2 / (sum(tmp.*tmp));
noise_var_estimate = var(channels_t{2}-channels_t{1}) / (2*(1/sqrt(2*pi*a22^2)));
noise_std_estimate = sqrt(noise_var_estimate);

%%
disp('Translation first empirical then ground truth');
[mean(all_trans) true_translation(2)]
disp('Standard deviation of translation, first empirical then according to theory');
[std(all_trans) SX]
disp('Noise standard deviation of translation, first ground truth then estimated');
[noise_std noise_std_estimate]





