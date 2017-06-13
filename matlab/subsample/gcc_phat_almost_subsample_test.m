%% "gcc-phat subsample"
% takes the gcc-phat "image", fits a polynomial of degree 2 to each of the
% different tie steps and then finds the maximum of these polynomials, to
% achieve a subsaample accuracy for the time delay.


addpath /home/gabbi/github/StructureFromSound/matlab/gcctracking_rex

if ~exist('a'),
    % First time load in sound data a and settings, but these are quite
    % big files
    %load spacecurea; % load a real sound experiment
    %load spacecuresettings; % load settings
    load blubba;
    load blubbsettings;
    xtid = settings.xtid;
    ymeter = settings.ymeter;
    % Calcualte GCC-PHAT scores
    settings.v = 340;                %speed of sound
    settings.mm = size(a,1);         %number of microphones
    settings.channels = 1:2;         %channels to read
    settings.refChannel = 1;         %reference channel
    settings.nbrOfSamples = length(a);    
    %% Correlation: GCC-PHAT
    settings.wf = @(x) 1./(abs(x)+(abs(x)<5e-3)); %weighting function
    settings.firstSamplePoint = 1; %center sample point of first frame
    settings.frameSize = 2048;     %width of frame in sample points
    settings.dx = 1000;            %distance between frames in sample points
    settings.frameOverlap = settings.frameSize-settings.dx; %overlap between frames
    settings.sw = 800;             %clipping of search window
    scores = gccscores(a,settings);
    
    %% Delays: Find highest peaks
    settings.nbrOfPeaks = 1;       %max number of peaks
    settings.minPeakHeight = 0.01; %min value of local maxima
    %Default: [4,0.01]
    
    u = getdelays(scores,settings);
    
end

% To supress the warning that the polynomial in polyfit below is bad
% condotioned.
id = 'MATLAB:polyfit:RepeatedPointsOrRescale';
warning('off',id)

scores_tmp = scores{1,2}; % only look at channel 1 and 2
u_tmp = u{1,2}; % the gcc-phat initial points
[m,n] = size(scores_tmp);

u_sub = nan(1,n);
poly_deg = 2;
range = 5; % how many points to use below + above in polyfit

indices = find(~isnan(u_tmp));

for i = indices
    middle = u_tmp(i)+settings.sw;
    middle_cut = scores_tmp(middle-range:middle+range,i)';
    fitted_poly = polyfit(middle-range:middle+range,middle_cut,poly_deg);
    if 1 
        figure(1); plot(scores_tmp(:,i)); hold on;
        plot(middle-15:0.1:middle+15,polyval(fitted_poly,middle-15:0.1:middle+15),'r--');
        hold off;
        pause();
    end
    sub_middle = roots(polyder(fitted_poly)); % change way of finding max?
    u_sub(i) = sub_middle-settings.sw;
end
 