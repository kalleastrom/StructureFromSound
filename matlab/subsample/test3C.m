% script to run the function estimate_trans_std_faster (essentially same as
% test3B_estimate_std_in_translation_faster) for several different
% settings. Add several for-loops etc to generate results and plot-

% settings
opts.nbr_channels = 2;
opts.noise_std = 0.03;
opts.signal_length = 1000;
opts.a1 = 0; % Blurr 1
opts.a2 = 2; % Blurr 2
opts.N = 100; % Nr of runs with different noise realizations
opts.nbr_decimals = 2; % how many decimals the translation should have
% opts.tt = [-15 15]; % add this? then remove in function

out = estimate_trans_std_faster(opts);

% todo
% - try changing some of the variables (noise_std, a2)
% - create nice plots using that
% 
many_a2 = linspace(0.3,3.5,33);
all_est_trans_std = zeros(1,length(many_a2));
all_th_trans_std = zeros(1,length(many_a2));
for i = 1:length(many_a2)   
    opts.noise_std = 0.03; % same through all iterations
    opts.a2 = many_a2(i); % Blurr 2
    out = estimate_trans_std_faster(opts);
    all_est_trans_std(i) = out.est_trans_std;
    all_th_trans_std(i) = out.theoretical_trans_std;
    
    % create a 2-sided 95% confidence interval for the empirical values
    trans_mean_e = out.est_trans;
    lower_e = trans_mean_e - 1.96*out.est_trans_std;
    higher_e = trans_mean_e + 1.96*out.est_trans_std;
    % compute the quote outliers/total number
    outliers_e(i) = sum(out.all_trans <= lower_e | ...
                    out.all_trans >= higher_e) / numel(out.all_trans);
    
    % create a 2-sided 95% confidence interval for the theoretical values
    trans_mean_t = out.true_trans;
    lower_t = trans_mean_t - 1.96*out.theoretical_trans_std;
    higher_t = trans_mean_t + 1.96*out.theoretical_trans_std;
    % compute the quote outliers/total number
    outliers_t(i) = sum(out.all_trans <= lower_t | ...
                    out.all_trans >= higher_t) / numel(out.all_trans);
end

close all;

figure(1); clf; 
subplot(2,1,1)
plot(many_a2, all_th_trans_std, 'r*', many_a2, all_est_trans_std, 'bx');
% title('std of the translation wrt noise std')
xlabel('noise std in signal'); ylabel('translation std')
legend('Theoretical std', 'Empirical std')
subplot(2,1,2)
plot(many_a2(3:end), all_th_trans_std(3:end), 'r*', ...
            many_a2(3:end), all_est_trans_std(3:end), 'bx');
axis([0 3.5 0.01 0.025])
% title('std of the translation wrt noise std')
xlabel('noise std in signal'); ylabel('translation std')
legend('Theoretical std', 'Empirical std')


