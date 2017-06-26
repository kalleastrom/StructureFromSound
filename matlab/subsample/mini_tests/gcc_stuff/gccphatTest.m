% to get a feeling for GCC-PHAT

t = linspace(0,50,1000);
x1 = sin(t);
x2 = sin(t+2);

figure(1); subplot(2,1,1); plot(t,x1); subplot(2,1,2); plot(t,x2);

% compute FFT of the signals
X1 = fft(x1);
X2 = fft(x2);

figure(2); subplot(2,1,1); plot(t,x1)
subplot(2,1,2); plot(t,abs(X1))

tmp= X1.*conj(X2);

gcc = ifft(tmp);
gccphat1 = ifft(tmp/norm(tmp));
gccphat2 = ifft(tmp./abs(tmp));

figure(3); subplot(3,1,1); plot(t,abs(gcc)); 
subplot(3,1,2); plot(t,gccphat1);
subplot(3,1,3); plot(t,gccphat2);
