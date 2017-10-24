% calls estimate_trans_dopp_std which simulates 2 channels and finds the
% translation and doppler for the second channel to be similar to the
% first. opts contains the settings and out the results.

close all;

% settings
opts.nbr_channels = 2;
opts.noise_std = 0.1;
opts.signal_length = 1000;
opts.a1 = 0;
opts.a2 = 2;
opts.N = 100;
opts.nbr_decimals = 2;
opts.tt = [-20 20];
plotopt = 0;

out = estimate_trans_dopp_std(opts,plotopt);
err = out.all_z - repmat(out.true_z,[1 opts.N]);
