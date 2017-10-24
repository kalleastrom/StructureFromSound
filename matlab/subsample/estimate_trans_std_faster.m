function out = estimate_trans_std_faster(opts,count)

% the same as test3B_estimate_std_in_translation_faster, but as a function
% with opts inparameter

%% First generate signals without noise

% settings
nbr_channels = opts.nbr_channels;
channels = cell(nbr_channels,1);
noise_std = opts.noise_std;
signal_length = opts.signal_length;
x = 1:signal_length;
a1 = opts.a1; % Blurr 1
a2 = opts.a2; % Blurr 2 
a22 = a2*sqrt(2);
N = opts.N; % Nr of runs with different noise realizations
if isfield(opts, 'tt')
    tt = opts.tt;
else
    tt = [-15 15]; % the translations to be tried
end

% generate ground truth translations
nbr_decimals = opts.nbr_decimals;  % how many decimals the translation should have
true_translation = round(2*10^(nbr_decimals+1)*rand(1,nbr_channels-1))/10^(nbr_decimals) - 10;
if isfield(opts,'true_translation')
    true_translation = opts. true_translation; % optional inparameter
end
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
    trans = zeros(nbr_channels,1);
    err = zeros(nbr_channels,1);
    for i = 2:nbr_channels
        if 1,
            % Version 1
            [curr_trans, curr_err] = find_translation(channels{1}, channels{i}, thresh, a2, tt);
        else
            % Version 2
            [curr_trans, curr_err] = find_translation2(channels{1}, channels{i}, thresh, a2, tt);
        end
        trans(i) = curr_trans(end); % if we want to use more than last trans we should have a look in tind_trans
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
SX = sqrt(VX);
Vb_approx = noise_std^2*4*sum(tmp.*tmp);
VX_approx = noise_std^2 / (sum(tmp.*tmp));
noise_var_estimate = var(channels_t{2}-channels_t{1}) / (2*(1/sqrt(2*pi*a22^2)));
noise_std_estimate = sqrt(noise_var_estimate);

%%
out.all_trans = all_trans; % The N different estimated translations
out.est_trans = mean(all_trans); % Translation empirical
out.true_trans = true_translation(2); % Translation ground truth
out.est_trans_std = std(all_trans); % Standard deviation of translation empirical
out.theoretical_trans_std = SX; % Standard deviation of translation according to theory
out.est_noise_std = noise_std; % Noise standard deviation estimated
out.true_noise_std = noise_std_estimate; % Noise standard deviation ground truth


