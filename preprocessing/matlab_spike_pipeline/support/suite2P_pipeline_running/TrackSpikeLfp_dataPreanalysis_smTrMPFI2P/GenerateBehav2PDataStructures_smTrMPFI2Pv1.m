function [Spike Clu Track Laps xml] = GenerateBehav2PDataStructures_smTrMPFI2Pv1(TwoPhotonPath, BehaviorPath, FileNameBase, isInt)

% function [Data Track Laps Spike xml] = GenerateBehavElectroDataStructures_MPFI(FileNameBase);
%
% This function reads the *.lfp  *.clu  *.res and *.whl files and
% produces data-structures with different sorting of behavioral and
% electrophysiol data:
% ouput structures:
%
%

%% 

% LOAD the 2P imaging .mat file here
% Behavior folder: 
% 'Z:\Xiaoliang\2P_free_moving\Zhao_data_2P\20211111\582\A582-20211111-04T\'
% 2P Data:
% 'Z:\Xiaoliang\2P_free_moving\Zhao_data_2P\20211111\582\2\suite2p\plane0\Fall.mat'
%
% GenerateBehav2PDataStructures_smTrMPFI2Pv1('Z:\Xiaoliang\2P_free_moving\Zhao_data_2P\20211111\582\2\suite2p\plane0\Fall.mat','Z:\Xiaoliang\2P_free_moving\Zhao_data_2P\20211111\582\A582-20211111-04T\','A582-20211111-04',0)

data_2p = load(TwoPhotonPath);
%% 
% analyze the recording

cd(BehaviorPath); 
    
%     load([FileNameBase '-param.mat']);
    
    % get recording file information from the xml file
    % if uncommmenting the following line, make sure to change method
    % signature to include FileNameBaseList
    % fileinfo = dir([FileNameBaseList '.dat']);
nChannels = 0;
        % total number of recording channels    
SampleRate = data_2p.ops.fs; %???
        % recording sampling rate (default: 20,000 Hz)
lfpSampleRate = 500;
        % local field potential sampling rate (default: 1250 Hz)
    
    % load behavior file
fprintf(['\nLoad behavior file and align the event time between '...
            'behavior and recording\n']);
check = Arduino2RecordT_smTr2P(FileNameBase,SampleRate,lfpSampleRate,data_2p);

if(check == 1)
    return;
end

% extract clusters and calculate corrected calcium fluorescence
fprintf('\n Extract clusters and calculate corrected calcium fluorescence... \n')
calCorrectedFluo(FileNameBase,data_2p,isInt);

% align tracking data with .dat file
fprintf('\nAlign tracking data with .dat file...\n')
aligntsp2dat_smTr2P(FileNameBase,SampleRate,lfpSampleRate,nChannels, data_2p);
        
% get tracking related data
fprintf('\nGet Track and Laps...\n')
SortTrials_3armMaze_smTr2P(FileNameBase,lfpSampleRate,SampleRate,isInt);

fprintf('\nSeperating Trials....\n')
GetTrials_smTrMPFI2P(TwoPhotonPath, BehaviorPath, FileNameBase);

fprintf('\nFinished!\n')
clearvars -except recordFolder path
    
end
