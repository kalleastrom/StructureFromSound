function [settings,scores,matches,rns,result] = main_gcctracking_rex(a,settings)
    % [settings] = main_gcctracking_rex();

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
    %% Delays: Find highest peaks
    settings.nbrOfPeaks = 4;       %max number of peaks
    settings.minPeakHeight = 0.01; %min value of local maxima
    %Default: [4,0.01]

    u = getdelays(scores,settings);
    %result.matchings.u = u;
    
    %% Clean up delays: Using channel-consistency constraints
    settings.binSize = 3;         %inlier threshold
    settings.minNbrOfInliers = 3; %min number of matching equations
    %Default values: [3,3]

    uref = matchingdelays(u,settings);

    %% RANSAC:
    settings.RANSACnbrOfIterations = 350;     %number of RANSAC iterations
    settings.RANSACframeSize = 21;            %line width
    settings.RANSACframeOverlap = 1;          %overlap of lines
    settings.RANSACmaxNbrOfGroups = 5;        %max nbr of lines
    settings.RANSACminNbrOfInliers = 6;       %min nbr of inliers
    settings.RANSACinlierThreshold = 1;     %max distance to line
    settings.RANSACsharedPointsThreshold = 2; %max nbr of shared points
    settings.RANSACmaxSlope = 3.5;            %max derivative
    %Default valuesOld: [350,21,1,5,5,1.5,2,3.5]
    %Default valuesOld: [350,21,1,5,7,1.5,2,3.5]
    %Default values: [350,21,1,5,6,1,2,3.5]

    [delays,lines,ind] = fitdelayswithransac(uref,settings);

    %% Connect lines:
    settings.lineDistanceThreshold = 10; %max vertical distance between lines
    %Default value: [10]

    [delaysegments,linesegments] = connectlines(delays,lines,ind,settings);

    %% Connect line-segments:
    settings.linesInlierThreshold = 10;
    settings.linesOverlap = 5;
    settings.linesInlierRatio = 0.15;
    [newdelaysegments,newlinesegments] = connectsegments(delaysegments,linesegments,ind,uref,settings);
    %Default values: [10,5,0.1]


    %% Smoothing:
    settings.smoothingDegree = 0.01; %degree of smoothing of uref
    settings.smoothingDistance = 2;  %set to 0 to keep smoothed curve
    %Default values: [0.01,2]

    urefsmooth = smoothdelays(newdelaysegments,newlinesegments,uref,settings);

    %% Clip data:
    uout = clipdata(urefsmooth,settings);

    % %% Save data:
    % if ~exist([fileNameBase 'data.mat'],'file')
    %     save([fileNameBase 'data'],'settings','u','uref','delays','lines','ind',...
    %         'delaysegments','linesegments','newdelaysegments','newlinesegments',...
    %         'urefsmooth','uout')
    % end

    %% Export data:

    %Additional settings for compatibility with tdoasystem_v5
    settings.in = length(uout);
    settings.swstep = 1;
    settings.wlist = [settings.frameSize];
    settings.scorefunction = 'score1';
    settings.xinorm = [1 2 4];
    settings.doplot = 1;
    settings.doverbose = 0;
    settings.doprint = 0;
    settings.nn = size(a,2);
    settings.tt = settings.nn/settings.sr;
    settings.isel = floor((1:(settings.in-1))*settings.nn/settings.in);
    settings.dsel = (-settings.sw):settings.swstep:settings.sw;
    settings.xtid = settings.isel/settings.sr;
    settings.ymeter = settings.dsel*settings.v/settings.sr;

    %Renaming:
    [~,leftLim] = ind2sub(size(uout),find(~isnan(uout),1,'first'));
    [~,rightLim] = ind2sub(size(uout),find(~isnan(uout),1,'last'));
    matches.uij = uout(:,leftLim:rightLim);
    matches.u = matches.uij*settings.v/settings.sr;
    matches.uindex = leftLim:rightLim;
    matches.utimes = matches.uindex*length(a)/settings.in/settings.sr;
    matches.uok = isfinite(matches.uij);

    rns = scores(settings.refChannel,:);
    
end

