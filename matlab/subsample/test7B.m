% test 7B
% Test to se if the estimations are better or worse if we also estimate
% Doppler and amplitude. Here, the signals have a translation and Doppler
% difference.

close all;

% settings
opts.nbr_channels = 2;
opts.noise_std = 0.02;
opts.signal_length = 1000;
opts.a1 = 0;
opts.a2 = 2;
opts.N = 100;
opts.nbr_decimals = 2;
opts.tt = [-20 20];
plotopt = 1;

% set the parameters
opts.true_translation = 3.63;
opts.true_doppler = 1.02;
opts.true_amplitude = 1.3;


%% Estimate only translation. Also std.

out_t = estimate_trans_std_2(opts,plotopt);

err_t = out_t.all_trans - repmat(out_t.true_trans,[1 opts.N]);


%% Estimate translation, Doppler. Also std.

out_td = estimate_trans_dopp_std_2(opts,plotopt);

err_td = out_td.all_z - repmat(out_td.true_z,[1 opts.N]);

%% Estimate translation, Doppler, amplitude. No std.

% To be continued

%% Write comparison
disp('Test where we had trans and dopp')
disp('[true trans, trans_t, est_std_t, th_est_t, trans_td, est_std_td th_std_td]')
disp(([opts.true_translation out_t.est_trans out_t.est_trans_std ...
    out_t.theoretical_trans_std out_td.est_z(1) ...
    out_td.est_z_std(1) out_td.theoretical_z_std(1)]))