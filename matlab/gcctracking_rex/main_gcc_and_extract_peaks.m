function [settings,scores,u] = main_gcc_and_extract_peaks(a,settings)
% [settings,scores,u] = main_gcc_and_extract_peaks(a,settings)

%% Read audio files:
settings.v = 340;                %speed of sound
settings.mm = size(a,1);         %number of microphones
settings.channels = 1:8;         %channels to read
settings.refChannel = 1;         %reference channel
%[340,8,1:8,1]

%[fileNameBase,dataDir,fileExt] = selectsoundfiles();
%[a,settings.sr] = readaudio([dataDir fileNameBase],fileExt,settings.mm,...
%    settings.channels);
settings.nbrOfSamples = length(a);

%% Correlation: GCC-PHAT
settings.wf = @(x) 1./(abs(x)+(abs(x)<5e-3)); %weighting function
settings.firstSamplePoint = 1; %center sample point of first frame
settings.frameSize = 2048;     %width of frame in sample points
settings.dx = 1000;            %distance between frames in sample points
settings.frameOverlap = settings.frameSize-settings.dx; %overlap between frames
settings.sw = 800;             %clipping of search window
%Default: [@(x) 1./(abs(x)+(abs(x)<5e-3)),1,2048,1048,800]

scores = gccscores(a,settings);

settings.nbrOfFrames = size(scores{1,1},2);
settings.xtid = (1:settings.nbrOfFrames)*settings.dx/settings.sr;

%Additional settings for compatibility with tdoasystem_v5
settings.swstep = 1;
settings.dsel = (-settings.sw):settings.swstep:settings.sw;
settings.ymeter = settings.dsel*settings.v/settings.sr;


%% Delays: Find highest peaks
settings.nbrOfPeaks = 4;       %max number of peaks
settings.minPeakHeight = 0.01; %min value of local maxima
%Default: [4,0.01]

u = getdelays(scores,settings);
%result.matchings.u = u;

