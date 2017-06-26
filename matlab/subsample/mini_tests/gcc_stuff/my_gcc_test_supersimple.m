% my gcc-test
%% with a gaussian

s1 = [1:10 9:-1:1 0 0 0];
s2 = [ 0 0 0 1:10 9:-1:1];
tt = 1:length(s1);

figure(1); subplot(2,1,1);
plot(tt,s1); subplot(2,1,2); plot(tt,s2);
subplot(2,1,1); title('Original signals 1')


%% Correlation:
wf1 = @(x) 1; %weighting function for gcc
wf2 = @(x) 1./(abs(x)+(abs(x)<5e-3)); %weighting function for gcc-phat

sw = 8;

scores1 = my_simple_gcc(s1,s2,wf1,sw);
scores2 = my_simple_gcc(s1,s2,wf2,sw);

tt2 = -floor(length(scores1)/2):floor(length(scores1)/2);

figure(2); subplot(2,1,1); plot(tt2,scores1); 
subplot(2,1,2); plot(tt2,scores2);
subplot(2,1,1); title('GCC and GCC-PHAT 1')

