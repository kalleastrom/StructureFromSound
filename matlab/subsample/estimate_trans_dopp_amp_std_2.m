function out = estimate_trans_dopp_amp_std_2(opts,plotopt)
% uses offset, in comparison to estimate_trans_dopp_amp_std.m
% only works for 2 channels so far. Otherwise the functions are the same.
% 2017-08-18

%% settings
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

%% generate ground truth translations, dopplers and amplitudes
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
% matrix to save all estimated trans, doppler and amplitudes in
all_z = zeros(3,N);
all_offset = zeros(2,N);

tic;
for kk = 1:N
    z = zeros(3,2); % matrix to save estimated trans and doppler in
    z(2:3,1) = 1; % the doppler and amplitude for signal 1 is set to 1
    % generate data and add noise 
    for i = 1:nbr_channels
        channels{i} = true_offset(i)+1/true_amplitude(i)*F0((1/true_doppler(i))*(x-true_translation(i)),a1);
        channels{i} = channels{i}+noise_std*randn(size(channels{i}));
    end 
    
    % estimate the translation, doppler and amplitude
    thresh = 10^(-8); % decides how good the translation estimation needs to be
    trans = zeros(nbr_channels,1);
    err = zeros(nbr_channels,1);
    offset = zeros(nbr_channels,1);
    offset(1) = mean(channels{1});
    for i = 2:nbr_channels    
        offset(i) = mean(channels{i});
        [curr_z, curr_err,f0t,f1t] = find_translation_doppler_amplitude(channels{1}-offset(1), channels{i}-offset(i), thresh, a2, tt);
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
        channels_t{i} = -offset(i)+z(3,i)*interp1d(channels{i},(z(2,i)*xmid+z(1,i)),a2); % obs!
        channels_mean = channels_mean+channels_t{i}./nbr_channels;       
        if plotopt
            figure(2); plot(xmid,channels{i}(xmid),colors(i));
            figure(3); plot(xmid,channels_t{i},colors(i));
        end
    end
    all_z(:,kk) = z(:,2);
    all_offset(:,kk) = offset;
    noise_var_estimate = var(channels_t{2}-channels_t{1}) / (2*(1/sqrt(2*pi*a22^2))); % obs! will this be the same?
    noise_std_estimate = sqrt(noise_var_estimate);
    all_noise_std(kk)=noise_std_estimate;
end

%% find the theoretical estimation of std. 
% OBS! Not correct yet. Gives negative variance for amplitude and time
% matters. The covariance of b is wrong. So is also the expected value of
% b (not zero) which makes the expression for the covariance of X more
% complicated.

cut_len = 100; 

% W = channels{1}; % use channels_mean instead?
W = channels_mean;

W_cut = W(cut_len+1:end-cut_len);
xx = -15:15; % for width of Gaussian kernel
gg = (-2*xx/(2*a2^2)).*(1/sqrt(2*pi*a2^2)).*exp(-(xx.^2)/(2*a2^2)); % Gaussian kernel. Make i tsum to 1?
W_prime = conv2(W,gg,'same');
W_prime_cut = W_prime(cut_len+1:end-cut_len); % only use middle part
t = 1:length(W); % what should t be? affects CX(2,2), but not CX(1,1), CX(3,3) 
t_cut = t(cut_len+1:end-cut_len);

% find expected value of A
EA = zeros(3); % 3x3 matrix for expected value of A 
EA(1,1) = 2*sum(W_prime_cut.^2);
EA(1,2) = 2*sum(t_cut.*W_prime_cut.^2);
EA(1,3) = 2*sum(W_cut.*W_prime_cut); % channels{1}\approx W. Use channels_mean instead?;
EA(2,1) = EA(1,2);
EA(2,2) = 2*sum(t_cut.^2.*W_prime_cut.^2);
EA(2,3) = 2*sum(t_cut.*W_cut.*W_prime_cut);
EA(3,1) = EA(1,3);
EA(3,2) = EA(2,3);
EA(3,3) = 2*sum(W_cut.^2+noise_std^2); %C[Esteck,Estreck]=noise_std^2

% find kernel for and r_{E-Estreck}
a22 = a2*sqrt(2);
gEEstreck = 2*noise_std^2*(1/sqrt(2*pi*a22^2)).*exp(-(xx.^2)/(2*a22^2)); % make it sum to 1?
W_prime_r_EEstreck = conv2(W_prime,gEEstreck,'same');
W_prime_r_EEstreck_cut = W_prime_r_EEstreck(cut_len+1:end-cut_len);

Y= W_prime.*t; % Y= W_prime(t_2)*t_2;
Y_r_EEstreck = conv2(Y,gEEstreck,'same'); 
Y_r_EEstreck_cut = Y_r_EEstreck(cut_len+1:end-cut_len);

% find covariance of b
CB = zeros(3); % 3x3 matrix for covariance of B
CB(1,1) = 4*sum(W_prime_cut.*W_prime_r_EEstreck_cut);
CB(1,2) = 4*sum(t_cut.*W_prime_cut.*W_prime_r_EEstreck_cut);
CB(1,3) = CB(1,1);
CB(2,1) = CB(1,2);
CB(2,2) = 4*sum(t_cut.*W_prime_cut.*Y_r_EEstreck_cut); %?????
CB(2,3) = CB(1,2);
CB(3,1) = CB(1,1);
CB(3,2) = CB(1,2);
CB(3,3) = CB(1,1);

% compute covariance of X
% CX = -(EA^(-1))^T*(-EA)^(-1)*CB = (EA^T)^(-1)*EA^(-1)*CB
CX = (EA')\(EA\CB);


%% save data to return

out.all_z = all_z; % the N different estimated trans, doppler and amplitude
out.est_z = mean(all_z,2); % empirical trans, doppler, amplitude
out.true_z = [true_translation(2); true_doppler(2); true_amplitude(2)]; % z ground truth
out.est_z_std = [std(all_z(1,:)); std(all_z(2,:)); std(all_z(3,:))]; % std of z empirical
out.theoretical_z_std = [];
% obs! are the following two rows correct? cmp estimate_trans_std_faster.m
out.est_noise_std = mean(all_noise_std); % Noise standard deviation estimated. obs! ok with mean?
out.true_noise_std = noise_std; % Noise standard deviation ground truth

