%% Main System Script

%% set system specific paths in a separate script
% Don't do this in this script. Run the script separately or in startup.m
% so that this script works for everyone.
%set_localpaths


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

[a,settings]=read_from_sfsdb(tdoasystemsettings,3); % Choose dataset here.
savedirname = 'data_3_system_2';                    % Remember to change savefolderpath, so
                                                    % that you don't
                                                    % overwrite olderr
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

% H�r skulle man nog s�tta saker som
% �r mikrofonerna i 3D eller i ett plan
% �r det ett kontinuerligt ljud eller klapp
% �r det en ljudk�lla eller flera olika samtidigt
% �r det en distinkt ljudk�lla eller utbredd
% �r det ett ekofritt rum eller inte

%% Run a system on the data
%settings for output
settings.doplot = 1;
settings.doverbose = 1;
settings.doprint = 0;

[result]=sfs_system_v1(a,settings);

%% Visualize result and compare with ground truth

if settings.hasgt,
    [errs]=plot_and_compare_result(result);     
else
    plot_result(result);
end

