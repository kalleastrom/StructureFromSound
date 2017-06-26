% my gcc-test
%% with a gaussian

tt = 1:2000;
sigma = 50;
s1 = 3*exp(-(tt-500).^2/(2*sigma^2));
s2 = 3*exp(-(tt-530).^2/(2*sigma^2));

a = [s1;s2];

figure(1); subplot(2,1,1);
plot(tt,s1); subplot(2,1,2); plot(tt,s2);
subplot(2,1,1); title('Original signals 1')


%% Correlation:
wf1 = @(x) 1; %weighting function for gcc
wf2 = @(x) 1./(abs(x)+(abs(x)<5e-3)); %weighting function for gcc-phat

sw = 800;

scores1 = my_simple_gcc(s1,s2,wf1,sw);
scores2 = my_simple_gcc(s1,s2,wf2,sw);

tt2 = -floor(length(scores1)/2):floor(length(scores1)/2);

figure(2); subplot(2,1,1); plot(tt2,scores1); 
subplot(2,1,2); plot(tt2,scores2);
subplot(2,1,1); title('GCC and GCC-PHAT 1')

%% with a gaussian plus noise

tt = 1:2000;
sigma = 50;
noise_std = 0.1;
s1 = 3*exp(-(tt-500).^2/(2*sigma^2)) + noise_std*randn(size(tt));
s2 = 3*exp(-(tt-530).^2/(2*sigma^2)) + noise_std*randn(size(tt));

a = [s1;s2];

figure(3); subplot(2,1,1);
plot(tt,s1); subplot(2,1,2); plot(tt,s2);
subplot(2,1,1); title('Original signals 2')


%% Correlation:
wf1 = @(x) 1; %weighting function for gcc
wf2 = @(x) 1./(abs(x)+(abs(x)<5e-3)); %weighting function for gcc-phat

sw = 800;

scores1 = my_simple_gcc(s1,s2,wf1,sw);
scores2 = my_simple_gcc(s1,s2,wf2,sw);

tt2 = -floor(length(scores1)/2):floor(length(scores1)/2);

figure(4); subplot(2,1,1); plot(tt2,scores1); 
subplot(2,1,2); plot(tt2,scores2);
subplot(2,1,1); title('GCC and GCC-PHAT 2')
%% with a sine

tt = 1:2000;
sigma = 50;
s1 = 3*sin(1/sigma*tt);
s2 = 3*sin(1/sigma*(tt+30));   

a = [s1;s2];

figure(5); subplot(2,1,1);
plot(tt,s1); subplot(2,1,2); plot(tt,s2);
subplot(2,1,1); title('Original signals 3')


%% Correlation:
wf1 = @(x) 1; %weighting function for gcc
wf2 = @(x) 1./(abs(x)+(abs(x)<5e-3)); %weighting function for gcc-phat

sw = 800;

scores1 = my_simple_gcc(s1,s2,wf1,sw);
scores2 = my_simple_gcc(s1,s2,wf2,sw);

tt2 = -floor(length(scores1)/2):floor(length(scores1)/2);

figure(6); subplot(2,1,1); plot(tt2,scores1); 
subplot(2,1,2); plot(tt2,scores2);
subplot(2,1,1); title('GCC and GCC-PHAT 3')


