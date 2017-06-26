% test_doppler
% calles estimare_dopp_std which simulates 2 channels and finds the doppler
% shift for the second channel. opts contains the settings, out the
% results.

close all;

% settings
opts.nbr_channels = 2;
opts.noise_std = 0.03;
opts.signal_length = 1000;
opts.a1 = 0;
opts.a2 = 2;
opts.N = 10;
opts.nbr_decimals = 2;
plotopt = 0;

out = estimate_dopp_std(opts,plotopt);
% err = out.all_z - repmat(out.true_z,[1 opts.N]);