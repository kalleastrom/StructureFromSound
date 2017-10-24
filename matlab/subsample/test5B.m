% Very similar to test5.m, but instead of having all settings inside the
% script the function estimate_trans_dopp_amp_std.m does the computations
% and takes all settings in the paramter opts. 
% See estimate_trans_dopp_amp_std for and test5.m further information.
% 2017-08-17

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

out = estimate_trans_dopp_amp_std(opts,plotopt);
err =out.all_z - repmat(out.true_z,[1 opts.N]);
mean(err,2)
max(err,[],2)
