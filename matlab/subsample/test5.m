% test 5 simulates 2 channels and finds the translation,
% doppler and a constant amplitude for the second channel (first used as 
% reference). It works sometimes, but is not very robust.

clear;
close all;

nbr_channels = 2;
channels = cell(nbr_channels,1);
noise_std = 0.1;

x = 1:2000;
a = 0;
channels{1} = F0(x,a);
channels{1} = channels{1} + noise_std*randn(size(channels{1}));

nbr_decimals = 2; % how many decimals the translation should have
true_translation = round(2*10^(nbr_decimals+1)*rand(1,nbr_channels-1))/10^(nbr_decimals) - 10;
true_translation = [0 true_translation];
% true_translation = [0 0];
true_doppler = 1 + 0.01*randn(1,nbr_channels-1);
true_doppler = [1 true_doppler];
% true_doppler = [1 1];
true_amplitude = abs(2*randn(1,nbr_channels-1));
true_amplitude = [1 true_amplitude];

z = zeros(3,2); % matrix to save estimated trans and doppler in
z(2:3,1) = 1; % the doppler and amplitude for signal 1 is set to 1

% create the channels and plot them
figure(1); clf; 
subplot(nbr_channels,1,1); plot(x,channels{1});
title(['The ' num2str(nbr_channels) ' different channels']);
for i = 2:nbr_channels
    channels{i} = 1/true_amplitude(i)*F0((1/true_doppler(i))*(x-true_translation(i)),a);
    channels{i} = channels{i}+noise_std*randn(size(channels{i}));
    subplot(nbr_channels,1,i); plot(x,channels{i});
end

thresh = 10^(-8); % decides how good the translation estimation needs to be
a = 2;
tt = [-15 15]; % the translations to be tried
trans = zeros(nbr_channels,1);
err = zeros(nbr_channels,1);
for i = 2:nbr_channels    
    [curr_z, curr_err] = find_translation_doppler_amplitude(channels{1}, channels{i}, thresh, a, tt);
    z(:,i) = curr_z;
    err(i) = curr_err(end);
end

% find the "original" signal be taking the mean of all signals, changed

xmid = x(100+1:end-100); % make sure that the cut part is larger than any translation
channels_t = cell(nbr_channels,1);
a = 2; % what should a be here?
channels_mean = zeros(1,length(xmid));
figure(2); clf; hold on;
figure(3); clf; hold on;
colors = ['b','g'];
for i = 1:nbr_channels
%     channels_t{i} = interp1d(channels{i},z(3,i)*(z(2,i)*xmid+z(1,i)),2); % obs!
    channels_t{i} = z(3,i)*interp1d(channels{i},(z(2,i)*xmid+z(1,i)),2); % obs!
%     subplot(nbr_channels,1,i); plot(x,channels{i});
    figure(2); plot(xmid,channels{i}(xmid),colors(i));
    figure(3); plot(xmid,channels_t{i},colors(i));
    % compute the mean signal
    channels_mean = channels_mean+channels_t{i}./nbr_channels;
end

% the mean will be a bit smoothed, while the original signals in Figure 2
% will not
figure(2); plot(xmid,channels_mean,'r--'); hold off; 
title('The original signals and the computed mean')
figure(3); plot(xmid,channels_mean,'r--'); hold off;
title('The translated signals and the computed mean')

disp(' ')
disp('Estimated translation, doppler and amplitude')
disp(z)
disp('True translation, doppler and amplitude')
disp([true_translation; true_doppler; true_amplitude]);
