%% Main System Script
%

%% set system specific paths in a separate script
% Don't do this in this script. Run the script separately or in startup.m
% so that this script works for everyone.
% set_localpaths

tdoasystemsettings.multipolpath = '/Users/kalle/Documents/projekt/github/multipol/';
tdoasystemsettings.datapath = '/Users/kalle/Documents/projekt/github/StructureFromSound/data/';

%% add paths

addpath(tdoasystemsettings.multipolpath);
addpath([tdoasystemsettings.datapath filesep 'sfsdb']);
addpath gcctracking_rex
addpath toa
addpath tdoa
addpath sfs_systems
addpath tools

%% Read audio files:
% Some settings are needed for file read and for setting up where temporary
% files will be saved.

[a,settings]=read_from_sfsdb(tdoasystemsettings,9); % Choose dataset here.
savedirname = 'tmp3';                    % Remember to change savefolderpath, so
% that you don't
% overwrite older
% results
settings.saveDir = [tdoasystemsettings.datapath 'savefiles' filesep savedirname filesep];
settings.gtDir = [tdoasystemsettings.datapath 'sfsgt' filesep settings.foldername filesep];
settings.gtFileName = [settings.gtDir 'gt.mat'];
if ~exist(settings.saveDir,'dir'),
    mkdir(settings.saveDir);
end
%% Read Ground Truth file (if there is one)
settings.hasgt = 0;
if exist(settings.gtDir,'dir'),
    if exist(settings.gtFileName,'file'),
        gt = load(settings.gtFileName);
        settings.gt = gt.gt;
        clear gt
        settings.hasgt = 1;
    end
end

% Maybe we should supply additional information here
% in the settings
% such as
% are the microphones in a plane or on a line
% are the sound events in a plane or on a line
% are the sound events continuous or discrete
% are the sound sources point-like or spread out
% is the room anaechoic or not


%% Run a system on the data
%settings for output
settings.doplot = 1;
settings.doverbose = 1;
settings.doprint = 0;

%[result]=sfs_system_v1(a,settings);

% Calculate gcc-phat score plots
[settings,scores,matches] = main_gcctracking_rex(a,settings);
rns = scores(settings.refChannel,:);
for k = 1:8,
    figure(1); clf;
    colormap(gray)
    imagesc(settings.xtid,settings.ymeter,scores{1,2})
    hold on
    plot(matches.utimes,matches.u(2,:),'.');
end;

% Load an initial estimate of
partcsavename = [settings.saveDir settings.fileNameBase 'C.mat'];
eval(['load ' partcsavename]);
[bmatches,bres] = estimate_o_from_matches(matches,settings);
[cmatches,cres] = estimate_sr_from_matcheso(bmatches,bres.o,bres,settings);
myplot(rns,cmatches,cres,settings);
[dmatches,dres] = find_more_inliers_uij(matches,rns,cmatches,cres,settings);
myplot(rns,dmatches,dres,settings);
%[speglingar,spegelplan]=finding_multipath_components(dres,dmatches,scores,settings);
result.settings = settings;
result.matches = dmatches;
result.res = dres;
result.scores = scores;
result.rns = rns;
%result.speglingar = speglingar;
%result.spegelplan = spegelplan;
result.debuginfo = [];
plot_result(result);

%% 
% Rensa bort lite av början

save scia_data dres dmatches

%%

[dmatches,dres] = find_more_inliers_uij(matches,rns,dmatches,dres,settings);
%[speglingar,spegelplan]=finding_multipath_components(dres,dmatches,scores,settings);
result.settings = settings;
result.matches = dmatches;
result.res = dres;
result.scores = scores;
result.rns = rns;
%result.speglingar = speglingar;
%result.spegelplan = spegelplan;
result.debuginfo = [];
plot_result(result);


%%

myplot(rns,dmatches,dres,settings);
inl = dmatches.uinliers;
u = dmatches.u;
[I,J]=find(inl);
D = u(find(inl));
x = real(dres.x);
y = real(dres.y);
o = real(dres.o);
%[xx4,yy4,oo4,res,jac]=bundletdoa(D,I,J,x,y,o,1,dres.xinorm);
[xx4,yy4,oo4,res,jac]=bundletdoasmooth(D,I,J,x,y,o,1,dres.xinorm);
% bundletdoasmooth(D,I,J,xt,yt,ot,debug,xinorm)
dres.x = xx4;
dres.y = yy4;
dres.o = oo4;
[dres] = recalc_residuals(dmatches,dres);

%%

myplot(rns,dmatches,dres,settings);
inl = dmatches.uinliers;
u = dmatches.u;
[I,J]=find(inl);
D = u(find(inl));
x = real(dres.x);
y = real(dres.y);
o = real(dres.o);
%
dfix = u - repmat(o,8,1);

x=toa_trilateration(dfix',y);

%[xx4,yy4,oo4,res,jac]=bundletdoa(D,I,J,x,y,o,1,dres.xinorm);
[xx4,yy4,oo4,res,jac]=bundletdoasmooth(D,I,J,x,y,o,1,dres.xinorm);
% bundletdoasmooth(D,I,J,xt,yt,ot,debug,xinorm)
dres.x = xx4;
dres.y = yy4;
dres.o = oo4;
[dres] = recalc_residuals(dmatches,dres);

y=toa_trilateration(d,x,y0,index,inliers);

%% Handcrafted tracker pieces from main_gcctracking
settings.nbrOfPeaks = 1;       %max number of peaks
settings.minPeakHeight = 0.01; %min value of local maxima
%Default: [4,0.01]
u = getdelays(scores,settings);
%
settings.binSize = 3;         %inlier threshold
settings.minNbrOfInliers = 3; %min number of matching equations
%Default values: [3,3]
uref = matchingdelays(u,settings);
xtid = [0 settings.xtid];
skalf = settings.v/settings.sr;
%
settings.RANSACnbrOfIterations = 350;     %number of RANSAC iterations
settings.RANSACframeSize = 21;            %line width
settings.RANSACframeOverlap = 1;          %overlap of lines
settings.RANSACmaxNbrOfGroups = 5;        %max nbr of lines
settings.RANSACminNbrOfInliers = 6;       %min nbr of inliers
settings.RANSACinlierThreshold = 3; %1    %max distance to line
settings.RANSACsharedPointsThreshold = 2; %max nbr of shared points
settings.RANSACmaxSlope = 7;            %max derivative
%Default valuesOld: [350,21,1,5,5,1.5,2,3.5]
%Default valuesOld: [350,21,1,5,7,1.5,2,3.5]
%Default values: [350,21,1,5,6,1,2,3.5]
[delays,lines,ind] = fitdelayswithransac(uref,settings);
%
% Extract those that were ok with ransac
urefransac = uref;
for k = 2:8;
    tmp = NaN*uref{k};
    for j = 1:size(delays{k},3);
        tmp(ind(:,j)')=delays{k}(1,:,j);
    end
    urefransac{1,k}=tmp;
end
% generate matches structure
uu = zeros(8,length(urefransac{2}));
for k = 2:8;
    uu(k,:)=urefransac{k};
end
tmp = sum(isfinite(uu));
uindex = find(tmp>=6);
matches.uij = uu(:,uindex);
matches.u = skalf*uu(:,uindex);
matches.uindex = uindex;
matches.utimes = xtid(uindex);
matches.uok = isfinite(uu(:,uindex));

%%
for k = 2:8,
    figure(k); clf;
    colormap(gray)
    imagesc(xtid,settings.ymeter,scores{1,k})
    hold on
    %plot(xtid,skalf*u{1,k},'y.');
    %plot(xtid,skalf*uref{1,k},'y.');
    plot(xtid,skalf*urefransac{1,k},'y.');
    %plot(xtid,skalf*urefsmooth{k},'y.');
    pause
end;


%% Visualize result and compare with ground truth

if settings.hasgt,
    [errs]=plot_and_compare_result(result);
else
    plot_result(result);
end

%   Copyright 2013-2017, Kalle Åström