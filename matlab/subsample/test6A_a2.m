% Similar to test3C.m, but uses two or three variables (trans+dopp or
% trans+dopp+amp). Estimates the std of both the variable estimations and
% the noise std and plots the empirically and theoretically estimated std:s
% of the variables as a function of the blurring parameter a2. Also

% close all;
%% To vary a2

clear opts;
nbr_vars = 2;

% settings
opts.nbr_channels = 2;
opts.noise_std = 0.03;
opts.signal_length = 1000;
opts.a1 = 0; % Blurr 1
opts.a2 = 2; % Blurr 2
opts.N = 1000; % Nr of runs with different noise realizations
opts.nbr_decimals = 2; % how many decimals the translation should have

% OBS! optional parameters
opts.tt = [-20 20]; % the translations to be tried. add in opts?

many_a2 = linspace(0.3,4,25);
% many_a2 = linspace(0.3,5,5);

all_est_z_std = zeros(nbr_vars,length(many_a2));
all_th_z_std = zeros(nbr_vars,length(many_a2));


% alpha is used to create the confidence interval and count outliers. Choose
% the correct alpha according to the width of the interval
% alpha = 1.96; % 95%
alpha = 2.58; % 99%
outliers_e = zeros(nbr_vars, length(many_a2));
outliers_t = zeros(nbr_vars, length(many_a2));
for i = 1:length(many_a2)
    i
    opts.a2 = many_a2(i); % Blurr 2
    if nbr_vars == 2
        out = estimate_trans_dopp_std(opts, 0); % , 1 for plots
    elseif nbr_vars == 3
        out = estimate_trans_dopp_amp_std_2(opts, 0); % 2 uses offset, 1 doesn't
    end
    all_est_z_std(:,i) = out.est_z_std;
    all_th_z_std(:,i) = out.theoretical_z_std;
    
    % create a 2-sided confidence interval for the empirical values
    mean_e = out.est_z; % estimated mean z 
    lower_e = mean_e - alpha*out.est_z_std; % estimated lower interval boundery for z
    higher_e = mean_e + alpha*out.est_z_std; % estimated higher intervl boundery for z
    % compute the quote outliers/total number
    outliers_e(:,i) = sum(out.all_z <= lower_e | ...
                        out.all_z >= higher_e, 2) / size(out.all_z,2);
    
    % create a 2-sided confidence interval for the theoretical values
    mean_t = out.true_z; % true "mean" z 
    lower_t = mean_t - alpha*out.theoretical_z_std; % theoretical lower interval boundery for z
    higher_t = mean_t + alpha*out.theoretical_z_std; % theoretical higher intervl boundery for z
    % compute the quote outliers/total number
    outliers_t(:,i) = sum(out.all_z <= lower_t | ...
                        out.all_z >= higher_t, 2) / size(out.all_z,2);
end


var_names = {'trans', 'dopp', 'amp'};

for i = 1:nbr_vars
    % plot the theoretical and empirical std for the  as a function
    % of the std of the noise (a2)
    figure(9+2*i); clf;
    subplot(2,1,1)
    plot(many_a2, all_th_z_std(i,:), 'r*', many_a2, all_est_z_std(i,:), 'bx');
    title(['std of the ' var_names{i} ' wrt a2'])
    xlabel('a2'); ylabel([var_names{i} ' std'])
    legend('Theoretical std', 'Empirical std')
    subplot(2,1,2)
    plot(many_a2(3:end), all_th_z_std(i,3:end), 'r*', ...
        many_a2(3:end), all_est_z_std(i,3:end), 'bx');
    axis([0 many_a2(end) 0.0 0.05])
    % title('std of the translation wrt a2')
    xlabel('a2'); ylabel([var_names{i} ' std'])
    legend('Theoretical std', 'Empirical std')
    
    % plot the ratio of outliers (size of confidence inteval chosen above) as a
    % function of the std of the noise (a2)
    figure(10+2*i); clf;
    plot(many_a2,outliers_t(i,:),'r*',many_a2,outliers_e(i,:),'bx');
    title(['ratio of ' var_names{i} ' outliers wrt a2'])
    xlabel('a2'); ylabel('ratio of outliers')
    legend('Using theoretical values', 'Using empirical values')
    
end

