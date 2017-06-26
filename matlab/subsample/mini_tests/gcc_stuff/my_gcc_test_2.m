%% my_gcc_test_2

addpath /home/gabbi/github/StructureFromSound/matlab/gcctracking_rex

% First time load in sound data a and settings, but these are quite
% big files
%load spacecurea; % load a real sound experiment
%load spacecuresettings; % load settings
tt = 1:2000;
sigma = 50;
noise_std = 0.1;
s1 = 3*exp(-(tt-500).^2/(2*sigma^2)) + noise_std*randn(size(tt));
s2 = 3*exp(-(tt-530).^2/(2*sigma^2)) + noise_std*randn(size(tt));

figure(1); subplot(2,1,1);
plot(tt,s1); subplot(2,1,2); plot(tt,s2);
subplot(2,1,1); title('Original signals 1')

a = [s1;s2];

% xtid = settings.xtid;
% ymeter = settings.ymeter;
% Calcualte GCC-PHAT scores
settings.v = 1;                %speed of sound
settings.mm = 2;         %number of microphones
settings.channels = 1:2;         %channels to read
settings.refChannel = 1;         %reference channel
settings.nbrOfSamples = length(a);    
%% Correlation: GCC-PHAT

settings.firstSamplePoint = 1; %center sample point of first frame
settings.frameSize = 4000;     %width of frame in sample points
settings.dx = 4000;            %distance between frames in sample points
settings.frameOverlap = settings.frameSize-settings.dx; %overlap between frames
settings.sw = 800;             %clipping of search window

wf1 = @(x) 1; %weighting function for gcc
wf2 = @(x) 1./(abs(x)+(abs(x)<5e-3)); %weighting function for gcc-phat

settings.wf = wf1; %weighting function

scores1 = gccscores(a,settings);

settings.wf = wf2; %weighting function

scores2 = gccscores(a,settings);

tt2 = -floor(length(scores1{1,2})/2):floor(length(scores1{1,2})/2);

figure(2); subplot(2,1,1); plot(tt2,scores1{2,1}(:,1)); 
subplot(2,1,2); plot(tt2,scores2{2,1}(:,1));
subplot(2,1,1); title('GCC and GCC-PHAT 1')



%% with a sine

tt = 1:2000;
sigma = 50;
s1 = 3*sin(1/sigma*tt);
s2 = 3*sin(1/sigma*(tt+30));   

a = [s1;s2];

figure(3); subplot(2,1,1);
plot(tt,s1); subplot(2,1,2); plot(tt,s2);
subplot(2,1,1); title('Original signals 2')


%% Correlation:

settings.wf = wf1; %weighting function

scores1 = gccscores(a,settings);

settings.wf = wf2; %weighting function

scores2 = gccscores(a,settings);

figure(4); subplot(2,1,1); plot(tt2,scores1{2,1}(:,1)); 
subplot(2,1,2); plot(tt2,scores2{2,1}(:,1));
subplot(2,1,1); title('GCC and GCC-PHAT 2');


