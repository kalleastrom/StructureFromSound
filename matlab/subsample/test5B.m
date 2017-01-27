close all;

% settings
opts.nbr_channels = 2;
opts.noise_std = 0.02;
opts.signal_length = 1000;
opts.a1 = 0;
opts.a2 = 2;
opts.N = 1000;
opts.nbr_decimals = 2;
opts.tt = [-20 20];
plotopt = 0;

out = estimate_trans_dopp_amp_std(opts,plotopt);
err =out.all_z - repmat(out.true_z,[1 opts.N]);
mean(err,2)
max(err,[],2)
