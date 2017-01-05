% script to run the function estimate_trans_std_faster (essentially same as
% test3B_estimate_std_in_translation_faster) for several different
% settings. Add several for-loops etc to generate results and plot-


close all;

%% To vary a2

clear opts;

% settings
opts.nbr_channels = 2;
opts.noise_std = 0.03;
opts.signal_length = 1000;
opts.a1 = 0; % Blurr 1
opts.a2 = 2; % Blurr 2
opts.N = 100; % Nr of runs with different noise realizations
opts.nbr_decimals = 2; % how many decimals the translation should have

% OBS! optional parameters
% opts.true_translation = 3.63; % the translation can be used as inparamter to test around the same value several times
opts.tt = [-20 20]; % the translations to be tried. add in opts?

% out = estimate_trans_std_faster(opts);

many_a2 = linspace(0.3,5,48);
all_est_trans_std = zeros(1,length(many_a2));
all_th_trans_std = zeros(1,length(many_a2));

% alpha is used to create the confidence interval and count outliers. Choose
% the correct alpha according to the width of the interval
% alpha = 1.96; % 95%
alpha = 2.58; % 99%
for i = 1:length(many_a2)
    opts.noise_std = 0.03; % same through all iterations
    opts.a2 = many_a2(i); % Blurr 2
    out = estimate_trans_std_faster(opts);
    all_est_trans_std(i) = out.est_trans_std;
    all_th_trans_std(i) = out.theoretical_trans_std;
    
    % create a 2-sided confidence interval for the empirical values
    trans_mean_e = out.est_trans;
    lower_e = trans_mean_e - alpha*out.est_trans_std;
    higher_e = trans_mean_e + alpha*out.est_trans_std;
    % compute the quote outliers/total number
    outliers_e(i) = sum(out.all_trans <= lower_e | ...
                    out.all_trans >= higher_e) / numel(out.all_trans);
    
    % create a 2-sided confidence interval for the theoretical values
    trans_mean_t = out.true_trans;
    lower_t = trans_mean_t - alpha*out.theoretical_trans_std;
    higher_t = trans_mean_t + alpha*out.theoretical_trans_std;
    % compute the quote outliers/total number
    outliers_t(i) = sum(out.all_trans <= lower_t | ...
                    out.all_trans >= higher_t) / numel(out.all_trans);

end


% plot the theoretical and empirical std for the translation as a function
% of the std of the noise (a2)
figure(10); clf; 
subplot(2,1,1)
plot(many_a2, all_th_trans_std, 'r*', many_a2, all_est_trans_std, 'bx');
% title('std of the translation wrt a2')
xlabel('a2'); ylabel('translation std')
legend('Theoretical std', 'Empirical std')
subplot(2,1,2)
plot(many_a2(3:end), all_th_trans_std(3:end), 'r*', ...
            many_a2(3:end), all_est_trans_std(3:end), 'bx');
axis([0 many_a2(end) 0.01 0.02])
% title('std of the translation wrt a2')
xlabel('a2'); ylabel('translation std')
legend('Theoretical std', 'Empirical std')

% plot the ratio of outliers (size of confidence inteval chosen above) as a
% function of the std of the noise (a2)
figure(11); clf;
plot(many_a2,outliers_t,'r*',many_a2,outliers_e,'bx');
title('ratio of outliers wrt a2')
xlabel('a2'); ylabel('ratio of outliers')
legend('Using theoretical values', 'Using empirical values')



%% To vary noise std
clear opts;
clear all;

% settings
opts.nbr_channels = 2;
opts.signal_length = 1000;
opts.a1 = 0; % Blurr 1
opts.a2 = 2; % Blurr 2
opts.N = 100; % Nr of runs with different noise realizations
opts.nbr_decimals = 2; % how many decimals the translation should have

% OBS! optional parameters
opts.true_translation = 3.63; % the translation can be used as inparamter to test around the same value several times
% opts.tt = [-20 20]; % the translations to be tried. add in opts?

% out = estimate_trans_std_faster(opts);

many_noise_std = linspace(0,3,25);
all_est_trans_std = zeros(1,length(many_noise_std));
all_th_trans_std = zeros(1,length(many_noise_std));

% alpha is used to create the confidence interval and count outliers. Choose
% the correct alpha according to the width of the interval
% alpha = 1.96; % 95%
alpha = 2.58; % 99%
for i = 1:length(many_noise_std)
    opts.noise_std = many_noise_std(i); 
    out = estimate_trans_std_faster(opts);
    all_est_trans_std(i) = out.est_trans_std;
    all_th_trans_std(i) = out.theoretical_trans_std;
    
    % create a 2-sided confidence interval for the empirical values
    trans_mean_e = out.est_trans;
    lower_e = trans_mean_e - alpha*out.est_trans_std;
    higher_e = trans_mean_e + alpha*out.est_trans_std;
    % compute the quote outliers/total number
    outliers_e(i) = sum(out.all_trans <= lower_e | ...
                    out.all_trans >= higher_e) / numel(out.all_trans);
    
    % create a 2-sided confidence interval for the theoretical values
    trans_mean_t = out.true_trans;
    lower_t(i) = trans_mean_t - alpha*out.theoretical_trans_std;
    higher_t(i) = trans_mean_t + alpha*out.theoretical_trans_std;
    % compute the quote outliers/total number
    outliers_t(i) = sum(out.all_trans <= lower_t(i) | ...
                    out.all_trans >= higher_t(i)) / numel(out.all_trans);
                
    % do this if the same true_translation is used throughout the script
    if isfield(opts,'true_translation')
        figure(51); 
        h(1) = plot(many_noise_std(i)*ones(1,opts.N),out.all_trans,'r.'); 
        hold on;
    end
end


% plot the different estimated translations for all iterations together
% with a line for the true translation. Optional to add the values used to
% find outliers.
if isfield(opts, 'true_translation')
    figure(51); hold on
    h(2) = plot([0 many_noise_std(end)], ...
                [opts.true_translation opts.true_translation],'-k');
    % uncomment below if you want to plot a line between the different
    % conficence intervals (different for different a2)
    h(3) = plot(many_noise_std,lower_t,'.k'); 
    plot(many_noise_std,higher_t,'.k'); 
    hold off;
    title('Estimated and true translation wrt noise std')
    xlabel('noise std'); ylabel('translation')
    legend(h,'Estimated translation','True translation','Confidence interval')
end

% plot the theoretical and empirical std for the translation as a function
% of the std of the noise (a2)
figure(12); clf; 
plot(many_noise_std, all_th_trans_std, 'r*', many_noise_std, all_est_trans_std, 'bx');
% title('std of the translation wrt noise std')
xlabel('noise std in signal'); ylabel('translation std')
legend('Theoretical std', 'Empirical std')


% plot the ratio of outliers (size of confidence inteval chosen above) as a
% function of the std of the noise (a2)
figure(13); clf;
plot(many_noise_std,outliers_t,'r*',many_noise_std,outliers_e,'bx');
title('ratio of outliers wrt noise std')
xlabel('noise std in signal'); ylabel('ratio of outliers')
legend('Using theoretical values', 'Using empirical values')

