function out = estimate_trans_dopp_amp_std(opts,plotopt)
% only works for 2 channels so far. Simulates two similar sound channels
% and then finds the translation, doppler and amplitude for the second
% channel to be similar to the first. opts contains the settings and
% signals and out contains the results. So far there is some problem with
% the amplitude estimation in find_translation_doppler_amplitude.m

% settings
nbr_channels = opts.nbr_channels;
channels = cell(nbr_channels,1);
noise_std = opts.noise_std;
signal_length = opts.signal_length;
x = 1:signal_length;
a1 = opts.a1;
a2 = opts.a2;
a22 = a2*sqrt(2);
N = opts.N;
if isfield(opts, 'tt')
    tt = opts.tt;
else
    tt = [-15 15];
end

% generate ground truth translations, dopplers and amplitudes
nbr_decimals = opts.nbr_decimals;  % how many decimals the translation should have
true_translation = round(2*10^(nbr_decimals+1)*rand(1,nbr_channels-1))/10^(nbr_decimals) - 10;
if isfield(opts,'true_translation')
    true_translation = opts.true_translation; % optional inparameter
end
true_translation = [0 true_translation];
true_doppler = 1 + 0.01*randn(1,nbr_channels-1);
if isfield(opts,'true_doppler')
    true_doppler = opts.true_doppler; % optional inparameter
end
true_doppler = [1 true_doppler];
true_amplitude = abs(2*randn(1,nbr_channels-1));
if isfield(opts,'true_amplitude')
    true_amplitude = opts.true_amplitude; % optional inparameter
end
true_amplitude = [1 true_amplitude];

% matrix to save all estimated trans, doppler and amplitudes in
all_z = zeros(3,N);

tic;
for kk = 1:N
    z = zeros(3,2); % matrix to save estimated trans and doppler in
    z(2:3,1) = 1; % the doppler and amplitude for signal 1 is set to 1
    % generate data and add noise 
    for i = 1:nbr_channels
        channels{i} = 1/true_amplitude(i)*F0((1/true_doppler(i))*(x-true_translation(i)),a1);
        channels{i} = channels{i}+noise_std*randn(size(channels{i}));
    end 
    
    % estimate the translation, doppler and amplitude
    thresh = 10^(-8); % decides how good the translation estimation needs to be
    trans = zeros(nbr_channels,1);
    err = zeros(nbr_channels,1);
    for i = 2:nbr_channels    
        [curr_z, curr_err,f0t,f1t] = find_translation_doppler_amplitude(channels{1}, channels{i}, thresh, a2, tt);
        z(:,i) = curr_z;
        err(i) = curr_err(end);
    end
    
    % find the "original" signal be taking the mean of all signals, changed
    xmid = x(100+1:end-100); % make sure that the cut part is larger than any translation
    channels_t = cell(nbr_channels,1);
    channels_mean = zeros(1,length(xmid));
    if plotopt
        figure(2); clf; hold on;
        figure(3); clf; hold on;
        colors = ['b','g'];
    end
    for i = 1:nbr_channels
        channels_t{i} = z(3,i)*interp1d(channels{i},(z(2,i)*xmid+z(1,i)),a2); % obs!
        channels_mean = channels_mean+channels_t{i}./nbr_channels;       
        if plotopt
            figure(2); plot(xmid,channels{i}(xmid),colors(i));
            figure(3); plot(xmid,channels_t{i},colors(i));
        end
    end
    all_z(:,kk) = z(:,2);
    noise_var_estimate = var(channels_t{2}-channels_t{1}) / (2*(1/sqrt(2*pi*a22^2))); % obs! will this be the same?
    noise_std_estimate = sqrt(noise_var_estimate);
    all_noise_std(kk)=noise_std_estimate;
end

% save data to return

out.all_z = all_z; % the N different estimated trans, doppler and amplitude
out.est_z = mean(all_z,2); % empirical trans, doppler, amplitude
out.true_z = [true_translation(2); true_doppler(2); true_amplitude(2)]; % z ground truth
out.est_z_std = [std(all_z(1,:)); std(all_z(2,:)); std(all_z(3,:))]; % std of z empirical
out.theoretical_z_std = []; % to be continued..
% obs! are the following two rows correct? cmp estimate_trans_std_faster.m
out.est_noise_std = mean(all_noise_std); % Noise standard deviation estimated. obs! ok with mean?
out.true_noise_std = noise_std; % Noise standard deviation ground truth

