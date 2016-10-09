function [a,settings]=read_from_sfsdb(tdoasystemsettings,id,full_dirname)
%[a,settings]=read_from_sfsdb(tdoasystemsettings,1);
if nargin < 3
    full_dirname = '';
end

switch id,
    case 1
        foldername = 'bassh1';
        settings.foldername = foldername;
        settings.fileNameBase = 'bassh1-';
        settings.dataDir = [tdoasystemsettings.datapath 'sfsdb' filesep foldername filesep];
        settings.fileExt = '.aiff';
        settings.mm = 8;         %number of microphones
        settings.channels = 1:8; %channels to read
        settings.micDim = 3; % The dimensionality of the span of the microphone positions
        settings.soundDim = 3; % The dimensionality of the span of the microphone positions
        [a,settings.sr] = readaudio([settings.dataDir settings.fileNameBase],settings.fileExt,settings.mm,settings.channels);
    case 2
        foldername = 'bassh2';
        settings.foldername = foldername;
        settings.fileNameBase = 'bassh2-';
        settings.dataDir = [tdoasystemsettings.datapath 'sfsdb' filesep foldername filesep];
        settings.fileExt = '.aiff';
        settings.mm = 8;         %number of microphones
        settings.channels = 1:8; %channels to read
        settings.micDim = 3; % The dimensionality of the span of the microphone positions
        settings.soundDim = 3; % The dimensionality of the span of the microphone positions
        [a,settings.sr] = readaudio([settings.dataDir settings.fileNameBase],settings.fileExt,settings.mm,settings.channels);
    case 3
        foldername = 'bassh3';
        settings.foldername = foldername;
        settings.fileNameBase = 'bassh3-';
        settings.dataDir = [tdoasystemsettings.datapath 'sfsdb' filesep foldername filesep];
        settings.fileExt = '.aiff';
        settings.mm = 8;         %number of microphones
        settings.channels = 1:8; %channels to read
        settings.micDim = 3; % The dimensionality of the span of the microphone positions
        settings.soundDim = 3; % The dimensionality of the span of the microphone positions
        [a,settings.sr] = readaudio([settings.dataDir settings.fileNameBase],settings.fileExt,settings.mm,settings.channels);
    case 4
        foldername = 'grieg1';
        settings.foldername = foldername;
        settings.fileNameBase = 'grieg1-';
        settings.dataDir = [tdoasystemsettings.datapath 'sfsdb' filesep foldername filesep];
        settings.fileExt = '.aiff';
        settings.mm = 8;         %number of microphones
        settings.channels = 1:8; %channels to read
        settings.micDim = 3; % The dimensionality of the span of the microphone positions
        settings.soundDim = 3; % The dimensionality of the span of the microphone positions
        [a,settings.sr] = readaudio([settings.dataDir settings.fileNameBase],settings.fileExt,settings.mm,settings.channels);
    case 5
        foldername = 'grieg2';
        settings.foldername = foldername;
        settings.fileNameBase = 'grieg2-';
        settings.dataDir = [tdoasystemsettings.datapath 'sfsdb' filesep foldername filesep];
        settings.fileExt = '.aiff';
        settings.mm = 8;         %number of microphones
        settings.channels = 1:8; %channels to read
        settings.micDim = 3; % The dimensionality of the span of the microphone positions
        settings.soundDim = 3; % The dimensionality of the span of the microphone positions
        [a,settings.sr] = readaudio([settings.dataDir settings.fileNameBase],settings.fileExt,settings.mm,settings.channels);
    case 6
        foldername = 'grieg3';
        settings.foldername = foldername;
        settings.fileNameBase = 'grieg3-';
        settings.dataDir = [tdoasystemsettings.datapath 'sfsdb' filesep foldername filesep];
        settings.fileExt = '.aiff';
        settings.mm = 8;         %number of microphones
        settings.channels = 1:8; %channels to read
        settings.micDim = 3; % The dimensionality of the span of the microphone positions
        settings.soundDim = 3; % The dimensionality of the span of the microphone positions
        [a,settings.sr] = readaudio([settings.dataDir settings.fileNameBase],settings.fileExt,settings.mm,settings.channels);
    case 7
        foldername = 'grieg4';
        settings.foldername = foldername;
        settings.fileNameBase = 'grieg4-';
        settings.dataDir = [tdoasystemsettings.datapath 'sfsdb' filesep foldername filesep];
        settings.fileExt = '.aiff';
        settings.mm = 8;         %number of microphones
        settings.channels = 1:8; %channels to read
        settings.micDim = 3; % The dimensionality of the span of the microphone positions
        settings.soundDim = 3; % The dimensionality of the span of the microphone positions
        [a,settings.sr] = readaudio([settings.dataDir settings.fileNameBase],settings.fileExt,settings.mm,settings.channels);
    case 8
        foldername = 'grieg5';
        settings.foldername = foldername;
        settings.fileNameBase = 'grieg5-';
        settings.dataDir = [tdoasystemsettings.datapath 'sfsdb' filesep foldername filesep];
        settings.fileExt = '.aiff';
        settings.mm = 8;         %number of microphones
        settings.channels = 1:8; %channels to read
        settings.micDim = 3; % The dimensionality of the span of the microphone positions
        settings.soundDim = 3; % The dimensionality of the span of the microphone positions
        [a,settings.sr] = readaudio([settings.dataDir settings.fileNameBase],settings.fileExt,settings.mm,settings.channels);
    case 9
        foldername = 'spacecure1';
        settings.foldername = foldername;
        settings.fileNameBase = 'space_cure-';
        settings.dataDir = [tdoasystemsettings.datapath 'sfsdb' filesep foldername filesep];
        settings.fileExt = '.aiff';
        settings.mm = 8;         %number of microphones
        settings.channels = 1:8; %channels to read
        settings.micDim = 3; % The dimensionality of the span of the microphone positions
        settings.soundDim = 3; % The dimensionality of the span of the microphone positions
        [a,settings.sr] = readaudio([settings.dataDir settings.fileNameBase],settings.fileExt,settings.mm,settings.channels);
    case 10
        foldername = 'whywereyouawayayearroy';
        settings.foldername = foldername;
        settings.fileNameBase = 'whywereyouawayayearroy-';
        settings.dataDir = [tdoasystemsettings.datapath 'sfsdb' filesep foldername filesep];
        settings.fileExt = '.aiff';
        settings.mm = 8;         %number of microphones
        settings.channels = 1:8; %channels to read
        settings.micDim = 3; % The dimensionality of the span of the microphone positions
        settings.soundDim = 3; % The dimensionality of the span of the microphone positions
        [a,settings.sr] = readaudio([settings.dataDir settings.fileNameBase],settings.fileExt,settings.mm,settings.channels);
    case 11
        foldername = 'axelf';
        settings.foldername = foldername;
        settings.fileNameBase = 'axelfA-';
        settings.dataDir = [tdoasystemsettings.datapath 'sfsdb' filesep foldername filesep];
        settings.fileExt = '.aiff';
        settings.mm = 8;         %number of microphones
        settings.channels = 1:8; %channels to read
        settings.micDim = 3; % The dimensionality of the span of the microphone positions
        settings.soundDim = 3; % The dimensionality of the span of the microphone positions
        [a,settings.sr] = readaudio([settings.dataDir settings.fileNameBase],settings.fileExt,settings.mm,settings.channels);
    case 12
        foldername = 'gone';
        settings.foldername = foldername;
        settings.fileNameBase = 'gone-';
        settings.dataDir = [tdoasystemsettings.datapath 'sfsdb' filesep foldername filesep];
        settings.fileExt = '.aiff';
        settings.mm = 8;         %number of microphones
        settings.channels = 1:8; %channels to read
        settings.micDim = 3; % The dimensionality of the span of the microphone positions
        settings.soundDim = 3; % The dimensionality of the span of the microphone positions
        [a,settings.sr] = readaudio([settings.dataDir settings.fileNameBase],settings.fileExt,settings.mm,settings.channels);
    otherwise
%         dirname = ensure_dirname(full_dirname,'sound files');
%         fileExtension = '.aiff';
%         sound_files = dir([dirname '*' fileExtension]);
%         settings.mm = length(sound_files);
%         settings.channels = 1:settings.mm;
%         
%         for i = 1:1
%             sound_file = sound_files(i);
%             ainfo = audioinfo(sound_file.name);
%             settings.sr = ainfo.SampleRate;
%             a = zeros(settings.mm,ainfo.TotalSamples);
%             a(i,:) = audioread(sound_file.name);
%         end
%         
%         for i = 2:settings.mm
%             sound_file = sound_files(i);
%             a(i,:) = audioread(sound_file.name);
%         end
end
        
end

