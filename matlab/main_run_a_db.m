%% Main script for evaluating against a (small) database of examples
% with ground truth

systems = {'sfs_system_v1'};
n_systems = length(systems);
examples = [1 2 3 8];
n_examples = length(examples);


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

%% The evaluation part.

for example_k = 1:n_examples;     % Loop through all n_examples
    for system_k = 1:n_systems;   % Loop through all n_systems
        
        % Read audio files:
        % Some settings are needed for file read and for setting up where temporary
        % files will be saved.
        [a,settings]=read_from_sfsdb(tdoasystemsettings,example_k); % Choose dataset here.
        
        
        % Remember to change savedirname, so
        % that you don't
        % overwrite older
        % results
        savedirname = ['run_20160719_data_' num2str(example_k) '_system_' num2str(system_k)];
        
        % set up the save directory
        settings.saveDir = [tdoasystemsettings.datapath 'savefiles' filesep savedirname filesep];
        if ~exist(settings.saveDir,'dir'),
            mkdir(settings.saveDir);
        end
        
        % Read Ground Truth file (if there is one)
        settings.gtDir = [tdoasystemsettings.datapath 'sfsgt' filesep settings.foldername filesep];
        settings.gtFileName = [settings.gtDir 'gt.mat'];
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
        
        %settings for output
        settings.doplot = 1;
        settings.doverbose = 1;
        settings.doprint = 0;
        
        % Run a system on the data
        try
            result = feval(systems{system_k},a,settings);
            %[result]=sfs_system_v1(a,settings);
            % Visualize result and compare with ground truth
            if settings.hasgt,
                [errs]=plot_and_compare_result(result);
            else
                plot_result(result);
            end
            
        catch
            result = NaN
            errs = NaN;
        end
        
        allerrs(system_k,example_k).errs = errs;
    end
end
