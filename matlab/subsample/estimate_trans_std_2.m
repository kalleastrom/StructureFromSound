function out = estimate_trans_std_2(opts,plotopt)
% Simaulates a signal and a copy which is translated, has Doppler shift and
% amlitude changes as well as an offset. Noise is asses. Then, the offset 
% is estimated and subtracted first whereupon the translation is estimated.
% For a signal copy without Doppler and amplitude, set opts.true_doppler=1 
% and opts.true_amplitude=1. (Similar to test3_estimate_std_in_translation) 
% 2018-03-22

%% settings

nbr_channels = opts.nbr_channels;
if nbr_channels~=2
    warning('The function only works for two channels')
end
channels = cell(nbr_channels,1);
noise_std = opts.noise_std;
signal_length = opts.signal_length;
x = 1:signal_length;
a1 = opts.a1;
a2 = opts.a2;
a22 = a2*sqrt(2);
N = opts.N;
if isfield(opts, tt)
    tt = opts.tt;
else
    tt = [-15 15];
end

%% generate ground truth translationd, dopplers and amplitudes
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
true_offset = round(5*10*rand(1,nbr_channels))/10;


%% do the estimations
all_trans = zeros(1,N);
all_offset = zeros(nbr_channels,N);
all_noise_std = zeros(1,N);

for kk = 1:N
    trans = zeros(1,nbr_channels); % matric to save estimated trans in 
    
    % generate data and add noise
    for i = 1:nbr_channels
        channels{i} = true_offset(i)+1/true_amplitude(i)*F0((1/true_doppler(i))*(x-true_translation(i)),a1);
        channels{i} = channels{i}+noise_std*randn(size(channels{i}));
    end     
    
    % estimate the translation
    thresh = 10^(-8); % decides how good the translation estimation needs to be
    err = zeros(1,nbr_channels);
    offset = zeros(nbr_channels,1);
    offset(1) = mean(channels{1});
    for i = 2:nbr_channels
        offset(i) = mean(channels{i});
        [curr_trans, curr_err] = find_translation(channels{1}-offset(1), channels{i}-offset(i), thresh, a2, tt);
        trans(i) = curr_trans(end);
        err(i) = curr_err;
    end
    
    % find the "original" signal by taking the mean of all signals
    xmid = x(100+1:end-100); % make sure that the cut part is larger than any translation
    channels_t = cell(nbr_channels,1);
    channels_mean = zeros(1,length(xmid));
    if plotopt
        figure(2); clf; hold on;
        figure(3); clf; hold on;
        colors = ['b','g'];
    end
    for i = 1:nbr_channels
        channels_t{i} = -offset(i)+interp1d(channels{i},xmid-trans(i),a2);
        channels_mean = channels_mean+channels_t{i}./nbr_channels;
        if plotopt
            figure(2); plot(xmid,channels{i}(xmid),colors(i));
            figure(3); plot(xmid,channels_t{i},colors(i));
            if i~=1; pause(); end
        end
    end
    
    all_trans(:,kk) = trans(:,2);
    all_offset(:,kk) = offset;
    noise_var_estimate = var(channels_t{2}-channels_t{1}) / (2*(1/sqrt(2*pi*a22^2)));
    noise_std_estimate = sqrt(noise_var_estimate);
    all_noise_std(kk)=noise_std_estimate;
end


mean(all_trans)
std(all_trans)
true_translation(2)
% formel f�r vad vi tror std ska bli


%% Find the theoretical (and estimated) estimation of std

cut_len = 100;

W = channels_mean;
W_cut = W(cut_len+1:end-cut_len);

xx = -15:15; % for width of Gaussian kernel. Might be changed if needed.
gg = (-2*xx/(2*a2^2)).*(1/sqrt(2*pi*a2^2)).*exp( - (xx.^2)/(2*a2^2) ); % derivative of Gaussian kernel
W_prime = conv2(W,gg,'same'); % "differentiate" W
W_prime_cut = W_prime(cut_len+1:end-cut_len); % only use middle part
% expected value of A
EA = 2*sum(W_prime_cut.^2);

gEEstreck = 2*noise_std^2*(1/sqrt(2*pi*a22^2)).*exp(-(xx.^2)/(2*a22^2)); 
W_prime_r_EEstreck = conv2(W_prime,gEEstreck,'same');
W_prime_r_EEstreck_cut = W_prime_r_EEstreck(cut_len+1:end-cut_len);

VB = 4*sum(W_prime_cut.*W_prime_r_EEstreck_cut);
VX = VB/(EA^2);
SX = sqrt(VX);

%% save data to return

out.cx = CX;
out.all_trans = all_trans;
out.est_trans = mean(all_trans);
out.true_trans = true_translation(2);
out.est_trans_std = std(all_trans);
out.theoretical_trans_std = SX;
out.est_noise_std = mean(all_noise_std);
out.true_noise_std = noise_std;


