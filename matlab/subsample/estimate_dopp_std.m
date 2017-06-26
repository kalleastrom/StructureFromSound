function out = estimate_dopp_std(opts,plotopt)
% Similar to estimate_trans_dopp_std, but only estimates the doppler
% effect, i.e. no translation or amplitude. An offset is though
% implemented. V = \bar{V}(alpha*t). Computes both the empirical std and
% the analytic estimation of the std for alpha. Implemented for two
% channels.

%% settings
nbr_channels = opts.nbr_channels;
if nbr_channels~=2
    warning('The function only works for two channels');
end
channels = cell(nbr_channels,1);
noise_std = opts.noise_std;
signal_length = opts.signal_length;
x = 1:signal_length;
a1 = opts.a1; % the blurr from recordings
a2 = opts.a2; % the blurr added in interpolation
a22 = a2*sqrt(2);
N = opts.N;

%% genereate ground truth dopplers
true_doppler = 1 + 0.01*randn(1,nbr_channels-1);
if isfield(opts, 'true_doppler')
    true_doppler = opts.true_doppler;
end
true_doppler = [1 true_doppler];
true_offset = round(5*10*rand(1,nbr_channels))/10;

%% do the estimations
% matrix to save all the estimated doppler and offsets in
all_dopp = ones(2,N);
all_offset = zeros(2,N);

for kk = 1:N
   % generate data and add noise
   for i = 1:nbr_channels
       channels{i} = true_offset(i) + F0(1/true_doppler(i)*x, a1);
       channels{i} = channels{i} + noise_std*randn(size(channels(i)));
   end
   
   % create a matrix to save the estimated doppler shifts in
   dopp = zeros(1,nbr_channels);
   offset = zeros(nbr_channels,1);
   err = zeros(nbr_channels,1);
   
   % estimate the offset and the doppler
    offset(1) = mean(channels{1});

    cut_length = 100;
    xmid = x(cut_length+1:end-cut_length);
    channels_mean_2 = zeros(1,length(xmid));
    thresh = 0.001;
    
    for i = 2:nbr_channels
        offset(i) = mean(channels{i});
        [curr_dopp, curr_err, f0t, f1t] = find_doppler(channels{1}-offset(1), channels{i}-offset(i), a2, thresh);
        dopp(:,i) = curr_dopp;
        err(i) = curr_err(end);
        channels_mean_2 = f0t./2 + f1t./2;
    end
    
    % find the "original" signal by taking the mean of all signals 
    cut_len = 100;
    xmid = x(cut_len+1:end-cut_len);
    channels_t = cell(nbr_channels,1);
    channels_mean = zeros(1,length(xmid));
    if plotopt
        figure(2); clf; hold on;
        figure(3); clf; hold on;
        colors = ['b','g'];
    end
    for i = 1:nbr_channels
        % obs! channels_t below should be the same as f0t and f1t above?
        % save some time by not recalculating!
        channels_t{i} = -offset(i)+interp1d(channels{i},dopp(i)*xmid,a2);
        channels_mean = channels_mean + channels_t{i}./nbr_channels;
        if plotopt
            figure(2); plot(xmid, channels{i}(xmid),colors(i));
            figure(3); plot(xmid, channels_t{i}, colors(i));
            if i~=1; pause(); end
        end
    end
    all_dopp(2,kk) = curr_dopp;
    all_offset(:,kk) = offset;
    noise_var_estimate = var(channels_t{2}-channels_t{1}) / (2*(1/sqrt(2*pi*a22^2)));
    noise_std_estimate = sqrt(noise_var_estimate);
    all_noise_std(kk) = noise_std_estimate;
end

%% find the theoretical estimation of std

cut_len = 100;
t_scale = 1; % converts time from samples to correct unit

W = channels_mean;

W_cut = W(cut_len+1:end-cut_len);
xx = -15:15; % for width of Gaussian kernel. Change if needed.
gg = (-2*xx/(2*a2^2)).*(1/sqrt(2*pi*a2^2)).*exp(-(xx.^2)/(2*a2^2)); % derivative of Gaussian kernel. Make it sum to 1?
W_prime = conv(W,gg,'same'); % derivative of W convolved with Gaussian \approx W convolved with derivative of Gaussian
W_prime_cut = W_prime(cut_len+1:end-cut_len);
t = 1:length(W); % time in samples
t = t*t_scale; % if time should be in other unit than samples
t_cut = t(cut_len+1:end-cut_len);

% find expected value of a in the Taylor expansion
EA = 2*sum((t_cut.^2) .* (W_prime_cut.^2) );

% find kernel for r_{E-Estreck}
gEEstreck  = 2*noise_std^2*(1/sqrt(2*pi*a22^2)).*exp(-(xx.^2)/(2*a22^2)); % make it sum to 1?
% W_prime_r_EEstreck = conv(W_prime,gEEstreck,'same');
% W_prime_r_EEstreck_cut = W_prime_r_EEstreck(cut_len+1:end-cut_len);

Y = W_prime.*t; % for later convolution. Y = W_prime(t_2)*t_2
Y_r_EEstreck = conv(Y,gEEstreck,'same');
Y_r_EEstreck_cut = Y_r_EEstreck(cut_len+1:end-cut_len);

% find the variance of b in the Taylor expansion
VB = 4 * sum( (t_cut.*W_prime_cut) .* Y_r_EEstreck_cut);

% compute the variance of X
VX = EA/(VB^2);
SX = sqrt(VX);

%% save data to return
% should the offset be returned as well?
% check sixe of all_dopp etc

out.all_dopp = all_dopp; % the N different estimated doppler shifts
out.est_dopp = mean(all_dopp,2); % empirical doppler
out.true_dopp = true_doppler(2); % ground truth
out.est_dopp_std = std(all_dopp); % std of dopp empirical
out.theoretical_dopp_std = SX;
% is the following correct?
out.est_noise_std = mean(all_noise_std); % Noise standard deviation estimated. obs! ok with mean?
out.teur_noise_std = noise_std; % Noise standard deviation ground truth


%Gabrielle, continue here!






















    