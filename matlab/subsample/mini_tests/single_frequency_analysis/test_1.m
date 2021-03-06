% to get a feeling for how single frequency analysis can be done

addpath('/home/gabbi/github/StructureFromSound/matlab/subsample')

tt = 1:2000;
freqHz = 50;
% s1 = 3*exp(-(tt-500).^2/(2*sigma^2));
% s2 = 3*exp(-(tt-530).^2/(2*sigma^2));

s1 = 3*sin(1/freqHz*tt);
s2 = 3*sin(1/freqHz*(tt+30));

a = [s1;s2];

figure(1); subplot(2,1,1);
plot(tt,s1); subplot(2,1,2); plot(tt,s2);
subplot(2,1,1); title('Original signals 1')
% calculate the DFT of the two signals
% S11 = fft(s1,2*length(s1)-1);
S1a = fft(s1);
S1b = fft(s1,length(s1));
S1c = fft(s1,1);
S1d = fft(s1,2*length(s1));

figure(10); subplot(2,2,1); plot(abs(S1a));
subplot(2,2,2); plot(abs(S1b));
subplot(2,2,3); plot(abs(S1b));
subplot(2,2,4); plot(abs(S1d));
% S21 = fft(s2,2*length(s2)-1);
S2 = fft(s2,length(s2));

% S1(2+k) seems to be the complex conjugate of S2(end-k)

freq_nbr = 200;
tmp = zeros(size(S1));
tmp(2+freq_nbr) = 1;
tmp(end-freq_nbr) = 1;
S1_part = S1.*tmp;
S2_part = S2.*tmp;

% how to plot phase and amplitude diagrams?
figure(10); subplot(2,1,1); plot(abs(S1));
subplot(2,1,2); plot(angle(S1))


s1_part = ifft(S1_part);
s2_part = ifft(S2_part);

figure(2); subplot(2,1,1);
plot(tt,s1_part); subplot(2,1,2); plot(tt,s2_part);
subplot(2,1,1); title('Signals after')

% now, find translation
thresh = [];
a2 = 2;
tt = [20 40];
[trans, err,f0t,f1t] = find_translation2(s1_part, s2_part, thresh, a2, tt);

figure(3);
plot(1:length(f0t),f0t); hold on; plot(1:length(f1t),f1t,'r--');

% problems with folding?
