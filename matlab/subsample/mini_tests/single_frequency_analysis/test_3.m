% test_3_notmine

amp = 1;
freqHz = 12000;
fsHz = 30000; % sampling frequency?
dt = 1/fsHz; %?
index = 1;
tt = 0:dt:1-dt;
yy = amp*sin(2*pi*freqHz*tt);

YY = fft(yy,fsHz)/fsHz;
magYY = abs(YY);
faxis = linspace(-fsHz/2,fsHz/2,fsHz);
plot(faxis/1000,magYY);
figure(); plot(faxis/1000,fftshift(magYY))