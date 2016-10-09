function [result]=sfs_system_v1(a,settings,gt);
% [result]=sfs_system_v1(a,settings,gt);
% INPUT:
% a - raw sound data - matrix of size mxn, 
%     where m is number of channels
%     and n is the number of sample points
% settings - 
%     Required information in settings
%     settings.mm = 8;         %number of microphones, not really required
%     settings.channels = 1:8; %channels to read
%     settings.refChannel = 1; %reference channel
% gt - optional ground truth information for benchmarking?????
% OUTPUT:
% result.matches - raw time-difference measurements and matching
% result.res     - 3D reconstruction of sound events and microphones
% result.scores  - gcc score images
% result.rns     - a subset of scores
% result.debuginfo - not used yet
%

%% Is ground truth provided. If so also do benchmarking against ground truth

if nargin<3,
    gt = NaN;
    hasgt = 0;
else
    hasgt = 1;
end;

%% More settings after the file has been read

% settings from tdoasystem intro
settings.v = 340; % Speed of sound
%settings.in = 1000; % nr of search positions (x-axis in rns-plot)
%settings.sw = 500; % Search width in sampling points (y-axis in rns-plot)
%settings.swstep = 3; % Search steps
%settings.wlist = [2000]; % Window size used in matching
%settings.scorefunction = 'score1';
settings.xinorm = [1 2 4];
% also from tdoasystem intro
settings.nn = size(a,2);
settings.tt = settings.nn/settings.sr;
%settings.isel = floor((1:(settings.in-1))*settings.nn/settings.in);
%settings.dsel = (-settings.sw):settings.swstep:settings.sw;
%settings.xtid = settings.isel/settings.sr;
%settings.ymeter = settings.dsel*340/settings.sr;


%% Here comes the system (with settings and a as input)

%% Run GCC-phat and matching (PART A)

% Saving the scores gives unreasonably large save files. Skip this. 
% partasavename = [settings.saveDir settings.fileNameBase 'A.mat'];
% if ~exist(partasavename),
%     [settings,scores,rawmatches] = main_gcctracking_rex(a,settings);
%     rns = scores(settings.refChannel,:);
%     eval(['save ' partasavename ' settings scores rawmatches rns']);
% else
%     eval(['load ' partasavename]);
% end
% Always recalculate scores.
[settings,scores,matches] = main_gcctracking_rex(a,settings);
rns = scores(settings.refChannel,:);

%settings.isel = floor((1:(settings.in-1))*settings.nn/settings.in);
%settings.dsel = (-settings.sw):settings.swstep:settings.sw;
%settings.xtid = settings.isel/settings.sr;
%settings.ymeter = settings.dsel*340/settings.sr;


%% Estimate offsets from raw matches (PART B)
%settings.in = 1999;
partbsavename = [settings.saveDir settings.fileNameBase 'B.mat'];
if ~exist(partbsavename),
    [bmatches,bres] = estimate_o_from_matches(matches,settings);
    eval(['save ' partbsavename ' bmatches bres']);
else
    eval(['load ' partbsavename]);
end

%% Estimate s and r (PART C)
partcsavename = [settings.saveDir settings.fileNameBase 'C.mat'];
if ~exist(partcsavename),
    [cmatches,cres] = estimate_sr_from_matcheso(bmatches,bres.o,bres,settings);
    eval(['save ' partcsavename ' cmatches cres']);
else
    eval(['load ' partcsavename]);
end
%cres = myplot(rns,cmatches,cres,settings);

%% Try to find more inliers among the remaining columns in uij

partdsavename = [settings.saveDir settings.fileNameBase 'D.mat'];
if ~exist(partdsavename),
    [dmatches,dres] = find_more_inliers_uij(matches,rns,cmatches,cres,settings);
     eval(['save ' partdsavename ' dmatches dres']);
else
    eval(['load ' partdsavename]);
end
%dres = myplot(rns,dmatches,dres,settings);

%%

%% Sub-pixel refinement (PART E)
% partesavename = [settings.saveDir settings.fileNameBase 'E.mat'];
% if ~exist(partesavename),
%     [ematches,eres] = sub_pixel_refinement(a,dmatches,dres,settings);    
%     eval(['save ' partesavename ' ematches eres']);
% else
%     eval(['load ' partesavename]);
% end
%
%eres = myplot(rns,ematches,eres,settings);
%%

[speglingar,spegelplan]=finding_multipath_components(dres,dmatches,scores,settings);


result.settings = settings;
result.matches = dmatches;
result.res = dres;
result.scores = scores;
result.rns = rns;
result.speglingar = speglingar;
result.spegelplan = spegelplan;
result.debuginfo = [];

%eval(['save ' finalsavename ' result']);

