%% To vary noise std
clear opts;
clear all;

nbr_vars = 2;

% settings
opts.nbr_channels = 2;
opts.signal_length = 1000;
opts.a1 = 0; % Blurr 1
opts.a2 = 2; % Blurr 2
opts.N = 1000; % Nr of runs with different noise realizations
opts.nbr_decimals = 2; % how many decimals the translation should have

% OBS! optional parameters
opts.true_translation = 3.63; % the translation can be used as inparamter to test around the same value several times
opts.true_doppler = 1.02;
opts.true_z = [opts.true_translation; opts.true_doppler];
if nbr_vars == 3
    opts.true_amplitude = 2.3;
    opts.true_z = [opts.true_z;opts.true_amplitude];
end


many_noise_std = linspace(0,3,25);
% many_noise_std = linspace(0,3,5);
all_est_z_std = zeros(nbr_vars,length(many_noise_std));
all_th_z_std = zeros(nbr_vars,length(many_noise_std));

% alpha is used to create the confidence interval and count outliers. Choose
% the correct alpha according to the width of the interval
% alpha = 1.96; % 95%
alpha = 2.58; % 99%
outliers_e = zeros(nbr_vars, length(many_noise_std));
outliers_t = zeros(nbr_vars, length(many_noise_std));

lower_t = zeros(nbr_vars, length(many_noise_std));
higher_t = zeros(nbr_vars, length(many_noise_std));
for i = 1:length(many_noise_std)
    i
    opts.noise_std = many_noise_std(i); 
    if nbr_vars == 2
        out = estimate_trans_dopp_std(opts, 0);
    elseif nbr_vars == 3
        out = estimate_trans_dopp_amp_std_2(opts, 0);
    end
    all_est_z_std(:,i) = out.est_z_std;
    all_th_z_std(:,i) = out.theoretical_z_std;
    
    % create a 2-sided confidence interval for the empirical values
    mean_e = out.est_z;
    lower_e = mean_e - alpha*out.est_z_std;
    higher_e = mean_e + alpha*out.est_z_std;
    % compute the quote outliers/total number
    outliers_e(:,i) = sum(out.all_z <= lower_e | ...
                    out.all_z >= higher_e, 2) / size(out.all_z,2);
    
    % create a 2-sided confidence interval for the theoretical values
    mean_t = out.true_z;
    lower_t(:,i) = mean_t - alpha*out.theoretical_z_std;
    higher_t(:,i) = mean_t + alpha*out.theoretical_z_std;
    % compute the quote outliers/total number
    outliers_t(:,i) = sum(out.all_z <= lower_t(:,i) | ...
                    out.all_z >= higher_t(:,i), 2) / size(out.all_z,2);
                
    % do this if the same true_translation is used throughout the script
    if isfield(opts,'true_translation')
        figure(51); 
        h(1,1) = plot(many_noise_std(i)*ones(1,opts.N),out.all_z(1,:),'r.'); 
        hold on;
        figure(52); 
        h(2,1) = plot(many_noise_std(i)*ones(1,opts.N),out.all_z(2,:),'r.'); 
        hold on;
        if nbr_vars == 3
            figure(53); 
            h(3,1) = plot(many_noise_std(i)*ones(1,opts.N),out.all_z(3,:),'r.'); 
            hold on;        
        end
    end
end


var_names = {'trans', 'dopp', 'amp'};

% plot the different estimated variables (trans+dopp or tran+dopp+anmp)
% together with a line for the true values. Optional to add the values used
% to find outliers
if isfield(opts, 'true_translation')
    for i = 1:nbr_vars
        figure(50+i); hold on;
        h(i,2) = plot([0 many_noise_std(end)],...
                    [opts.true_z(i) opts.true_z(i)],'-k');
        % uncomment below if you want to plot a line between the different
        % confidence intervals (different for different noise std)
        h(i,3) = plot(many_noise_std, lower_t(i,:),'.k');
        plot(many_noise_std, higher_t(i,:),'.k');
        hold off;
        title(['Estimated and true ' var_names{i} ' wrt noise std']);
        xlabel('noise std'); ylabel(var_names{i})
        legend(h(i,:),['Estimated ' var_names{i}],['True' var_names{i}],...
                    'Confidence interval');
    end
end

% plot the theoretical and empirical std for the vars as a function of the
% std of the noise
for i = 1:nbr_vars
    figure(15+2*i); clf;
    plot(many_noise_std, all_th_z_std(i,:),'r*', ...
            many_noise_std, all_est_z_std(i,:),'bx');
    title(['std of ' var_names{i} ' wrt noise std']);
    xlabel('noise std in signal'); ylabel([var_names{i} ' std']);
    legend('Theoretical std', 'Empirical std');
end

% plot the ratio of outliers (size of confidence interval chosen above) as
% a function of the std of the noise
for i = 1:nbr_vars
    figure(16+2*i); clf; 
    plot(many_noise_std,outliers_t(i,:),'r*',...
            many_noise_std,outliers_e(i,:),'bx');
    title(['Ratio of ' var_names{i} ' outliers wrt noise std']);
    xlabel('noise std in signal'); 
    ylabel(['ratio of ' var_names{i} ' outliers'])
    legend('Using theoretical values', 'Using empirical values');
end


%
% % plot the different estimated translations for all iterations together
% % with a line for the true translation. Optional to add the values used to
% % find outliers.
% if isfield(opts, 'true_translation')
%     figure(51); hold on
%     h(2) = plot([0 many_noise_std(end)], ...
%                 [opts.true_translation opts.true_translation],'-k');
%     % uncomment below if you want to plot a line between the different
%     % conficence intervals (different for different a2)
%     h(3) = plot(many_noise_std,lower_t,'.k'); 
%     plot(many_noise_std,higher_t,'.k'); 
%     hold off;
%     title('Estimated and true translation wrt noise std')
%     xlabel('noise std'); ylabel('translation')
%     legend(h,'Estimated translation','True translation','Confidence interval')
% end
% 
% % plot the theoretical and empirical std for the translation as a function
% % of the std of the noise (a2)
% figure(12); clf; 
% plot(many_noise_std, all_th_trans_std, 'r*', many_noise_std, all_est_trans_std, 'bx');
% % title('std of the translation wrt noise std')
% xlabel('noise std in signal'); ylabel('translation std')
% legend('Theoretical std', 'Empirical std')
% 
% 
% % plot the ratio of outliers (size of confidence inteval chosen above) as a
% % function of the std of the noise (a2)
% figure(13); clf;
% plot(many_noise_std,outliers_t,'r*',many_noise_std,outliers_e,'bx');
% title('ratio of outliers wrt noise std')
% xlabel('noise std in signal'); ylabel('ratio of outliers')
% legend('Using theoretical values', 'Using empirical values')
% 
