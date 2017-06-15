function out = estimate_trans_dopp_std(opts,plotopt)
% Similar to estimate_trans_dopp_amp_std_2 this function uses offset.
% Though the amplitude is not estimated, only doppler and translation, such
% that V = \bar{V}(z2*t+z1). Also compares to analytic estimations of the
% std (needs more work). Only for two channels.

%% setting
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
true_offset = round(5*10*rand(1,nbr_channels))/10;

%% do the estimations
% matrix to save all estimated trans and doppler in
all_z = zeros(2,N);
all_offset = zeros(2,N);

tic;
for kk = 1:N
    z = zeros(2,2); % matrix to save estimated trans and doppler in
    z(2,1) = 1; % the doppler for signal 1 is set to 1
    % generate data and add noise 
    for i = 1:nbr_channels
        channels{i} = true_offset(i)+F0((1/true_doppler(i))*(x-true_translation(i)),a1);
        channels{i} = channels{i}+noise_std*randn(size(channels{i}));
    end 
    
    % estimate the translation and doppler
    thresh = 10^(-8); % decides how good the translation estimation needs to be
    trans = zeros(nbr_channels,1);
    err = zeros(nbr_channels,1);
    offset = zeros(nbr_channels,1);
    offset(1) = mean(channels{1});
    for i = 2:nbr_channels    
        offset(i) = mean(channels{i});
        [curr_z, curr_err,f0t,f1t] = find_translation_and_doppler(channels{1}-offset(1), channels{i}-offset(i), thresh, a2, tt);
        z(:,i) = curr_z;
        err(i) = curr_err(end);
    end
    
    % find the "original" signal be taking the mean of all signals
    xmid = x(100+1:end-100); % make sure that the cut part is larger than any translation
    channels_t = cell(nbr_channels,1);
    channels_mean = zeros(1,length(xmid));
    if plotopt
        figure(2); clf; hold on;
        figure(3); clf; hold on;
        colors = ['b','g'];
    end
    for i = 1:nbr_channels
        channels_t{i} = -offset(i)+interp1d(channels{i}, z(2,i)*xmid+z(1,i), a2); 
        channels_mean = channels_mean+channels_t{i}./nbr_channels;       
        if plotopt
            figure(2); plot(xmid,channels{i}(xmid),colors(i));
            figure(3); plot(xmid,channels_t{i},colors(i));
            pause()
        end
    end
    all_z(:,kk) = z(:,2);
    all_offset(:,kk) = offset;
    noise_var_estimate = var(channels_t{2}-channels_t{1}) / (2*(1/sqrt(2*pi*a22^2))); % obs! will this be the same?
    noise_std_estimate = sqrt(noise_var_estimate);
    all_noise_std(kk)=noise_std_estimate;
end


%% find the theoretical estimation of std. 
% OBS! Not checked yet.

cut_len = 100; 
% freq = 96000;
% t_scale = 1/freq;
t_scale = 1;

% W = channels{1}; % use channels_mean below instead
W = channels_mean;

W_cut  = W(cut_len+1:end-cut_len);
xx = -15:15; % for width of Gaussian kernel. Might be changed if needed.
gg = (-2*xx/(2*a2^2)).*(1/sqrt(2*pi*a2^2)).*exp(-(xx.^2)/(2*a2^2)); % derivative of Gaussian kernel. Make it sum to 1?
% gg = gg./sum(gg);
W_prime = conv2(W,gg,'same');
W_prime_cut = W_prime(cut_len+1:end-cut_len); % only use middle part
t = 1:length(W); % time. Is this correct? Samples? Seconds? Other? Does it matter?
t = t*t_scale; % obs?!?!
t_cut = t(cut_len+1:end-cut_len);

% find the expected value of A
EA = zeros(2);
EA(1,1) = 2*sum(W_prime_cut.^2);
EA(1,2) = 2*sum( t_cut .* W_prime_cut.^2);
EA(2,1) = EA(1,2);
EA(2,2) = 2*sum( (t_cut.^2) .* (W_prime_cut.^2) );

% find kernel for r_{E-Estreck}
gEEstreck = 2*noise_std^2*(1/sqrt(2*pi*a22^2)).*exp(-(xx.^2)/(2*a22^2)); % make it sum to 1?
W_prime_r_EEstreck = conv2(W_prime,gEEstreck,'same');
W_prime_r_EEstreck_cut = W_prime_r_EEstreck(cut_len+1:end-cut_len);

Y = W_prime.*t; % for later convolution. Y = W_prime(t_2)*t_2;
Y_r_EEstreck = conv2(Y,gEEstreck,'same');
Y_r_EEstreck_cut = Y_r_EEstreck(cut_len+1:end-cut_len);
YY = W_prime.*t.^2; % YY = W_prime(t_2)*t_2^2
YY_r_EEstreck = conv2(YY,gEEstreck,'same');
YY_r_EEstreck_cut = YY_r_EEstreck(cut_len+1:end-cut_len);

% find covariance of b
CB = zeros(2);
% CB(1,1) = 4* sum(W_prime_cut .* W_prime_r_EEstreck_cut);
% CB(1,2) = 4* sum( (t_cut.*W_prime_cut) .*  Y_r_EEstreck_cut);
% CB(2,1) = CB(1,2);
% CB(2,2) = 4* sum( (t_cut.^2 .* W_prime_cut) .* YY_r_EEstreck_cut );

CB(1,1) = 4* sum(W_prime_cut .* W_prime_r_EEstreck_cut);
CB(1,2) = 4* sum( (W_prime_cut) .*  Y_r_EEstreck_cut); % should be same as CB(2,1)
CB(2,1) = 4* sum( (t_cut.*W_prime_cut) .*  W_prime_r_EEstreck_cut); % should be same as CB(1,2)
CB(2,2) = 4* sum( (t_cut .* W_prime_cut) .* Y_r_EEstreck_cut );

% compute covariance of X
% CX = -(EA^(-1))^T*(-EA)^(-1)*CB = (EA^T)^(-1)*EA^(-1)*CB
% CX = (EA')\(EA\CB);
% CX = EA\CB/EA'; % see p 34 in http://users.isy.liu.se/en/rt/hendeby/files/PhD1161.pdf
CX = EA'\CB/EA;
% CX = EA'\EA\CB; % same as (EA'\EA)\CB, becomes real

% OBS! Haven't checked part below!!!
%% save data to return

out.cx = CX;
out.all_z = all_z; % the N different estimated trans, doppler and amplitude
out.est_z = mean(all_z,2); % empirical trans, doppler
out.true_z = [true_translation(2); true_doppler(2)]; % z ground truth
out.est_z_std = [std(all_z(1,:)); std(all_z(2,:))]; % std of z empirical
out.theoretical_z_std = [sqrt(CX(1,1)); sqrt(CX(2,2))];
% obs! are the following two rows correct? cmp estimate_trans_std_faster.m
out.est_noise_std = mean(all_noise_std); % Noise standard deviation estimated. obs! ok with mean?
out.true_noise_std = noise_std; % Noise standard deviation ground truth


