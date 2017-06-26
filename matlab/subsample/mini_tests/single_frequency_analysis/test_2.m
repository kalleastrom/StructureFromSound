% test_2

% to get a feeling for how single frequency analysis can be done

addpath('/home/gabbi/github/StructureFromSound/matlab/subsample')

amp = 3;
freqHz = 1/50; % frequency of signal
fsHz = 46000; % sampling frequency
dt = 1/fsHz; % seconds between two samples
 
% tt = hmm?
tt = -1000:1000;
% tt = 1:dt:2000;
% tt = 0:2000;

h = 30;

% s1 = amp*sin(1/freqHz*tt);
s1 = amp*sin(2*pi*freqHz*tt);
s2 = amp*sin(2*pi*freqHz*tt+h);

figure(1); subplot(2,1,1);
plot(tt,s1); subplot(2,1,2); plot(tt,s2);
subplot(2,1,1); title('Original signals 1')
% calculate the DFT of the two signals
% S11 = fft(s1,2*length(s1)-1);
% S1a = fft(s1);
% S1b = fft(s1,length(S2));
% S1c = fft(s1,1);
% S1d = fft(s1,2*length(S2));

S1 = fft(s1,length(s1))/length(s1);
S2 = fft(s2,length(s2))/length(s2);

% figure(10); subplot(2,2,1); plot(abs(S1a));
% subplot(2,2,2); plot(abs(S1b));
% subplot(2,2,3); plot(abs(S1b));
% subplot(2,2,4); plot(abs(S1d));
% % S21 = fft(s2,2*length(s2)-1);
% S2 = fft(s2,length(s2));

% S1(2+k) seems to be the complex conjugate of S2(end-k)

freq_nbr = 200;
tmp = zeros(size(S1));
tmp(2+freq_nbr) = 1;
tmp(end-freq_nbr) = 1;
S1_part = S1.*tmp;
S2_part = S2.*tmp;

% how to plot phase and amplitude diagrams?
figure(10); subplot(2,1,1); plot(tt,fftshift(abs(S1)));
subplot(2,1,2); plot(tt,fftshift(angle(S1)))


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
