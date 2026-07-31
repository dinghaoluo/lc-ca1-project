function [Spike Clu Track Laps xml] = GenerateBehavElectroDataStructures_MPFIv1(path, FileNameBaseList, varargin)

% function [Data Track Laps Spike xml] = GenerateBehavElectroDataStructures_MPFI(FileNameBase);
%
% This function reads the *.lfp  *.clu  *.res and *.whl files and
% produces data-structures with different sorting of behavioral and
% electrophysiol data:
% ouput structures:
%
%
% GenerateBehavElectroDataStructures_JF_v2(...
% '/groups/pastalkova/pastalkovalab/data/eva_recordings/A129/A129-20111104',
% 'A129-20111104','A129-20111104-01','A129-20111104-02')
% GenerateBehavElectroDataStructures_JF_v2(...
% '/groups/pastalkova/pastalkovalab/data/eva_recordings/A129/A129-20111104',
% 'A129-20111104-01')

%% 
% if list of input files = concatenated recording => copy concatenated data 
% to the folders with the original recordings
if ~isempty(varargin)    
    concatFolder = FileNameBaseList;
    for nf = 1 : length(varargin)
        recordFolder{nf} = varargin{nf};
    end
    nTotSampl = redistrClu(path, concatFolder, recordFolder);    
else
    concatFolder = [];
    recordFolder{1} = FileNameBaseList;   
    cd([path '/' FileNameBaseList '/']);
    if(exist([FileNameBaseList '.dat.orig'],'file') == 2)
        system(['rm ' FileNameBaseList '.dat']);
        system(['ln -s ' FileNameBaseList '.dat.orig ' ...
                    FileNameBaseList '.dat']);
    end  
end

cd(path);

%% 
% analyze data recordings one by one
for nf = 1 : size(recordFolder,2)
    
    FileNameBase = recordFolder{nf}; 
    cd([path '/' FileNameBase]); 
    
%     load([FileNameBase '-param.mat']);
    
    % load .xml file
    if exist([FileNameBase '.xml'], 'file') == 2
        [xml] = LoadXml_e([FileNameBase '.xml']);
    else fprintf('\n Cannot find %s - file -> analysis aborted.\n', ...
                [FileNameBase '.xml']);
        return;
    end
    
    % get recording file information from the xml file
    fileinfo = dir([FileNameBaseList '.dat']);
    nChannels = xml.nChannels;
        % total number of recording channels    
    SampleRate = xml.SampleRate;
        % recording sampling rate (default: 20,000 Hz)
    lfpSampleRate = xml.lfpSampleRate;
        % local field potential sampling rate (default: 1250 Hz)
    
    % load behavior file
    fprintf(['\nLoad behavior file and align the event time between',...
            'behavior and recording\n']);
    Arduino2RecordT_mpfi(FileNameBase,SampleRate,lfpSampleRate,nChannels);
    
    % load stimulation file
    fprintf(['\nLoad stimulation file and align the event time between'...
            'stimulation and recording\n']);
    Stim2RecordT_mpfi(FileNameBase,SampleRate,lfpSampleRate,nChannels);
    
    % align tracking data with .dat file
    fprintf('\nAlign tracking data with .dat file...\n')
    aligntsp2dat_mpfi(FileNameBase,SampleRate,lfpSampleRate,nChannels);
        
    % get tracking related data
    fprintf('\nGet Track and Laps...\n')
    purpose = 0;
    load([FileNameBase '-param.mat'],'epochGroup');
    if(~isempty(strfind(epochGroup.protocol.purpose,'Treadmill')))
        purpose = 1; % Treadmill
    end
    SortTrials_3armMaze_mpfi(FileNameBase, lfpSampleRate, purpose);
    
    % get spikes related data
    fprintf('\nGet Spike and Clu....\n')
    getSpikes_mpfi(FileNameBase, SampleRate, lfpSampleRate); 
    
    % calculate CCG
    fprintf('\nCalculate CCG\n');
    calACG_mpfi(FileNameBase,lfpSampleRate);
    
    % extract theta phase
    fprintf('\nExtract theta phase...\n')
    % Generate _thetaEst_Ref.dat if necessary, if there is no artifact in
    % the recording, then use .dat file directly
%     if(exist('ChListRef.txt', 'file') == 2)
%         datFileName = [FileNameBase '-thetaEst_Ref.dat'];
%         if exist(datFileName, 'file') ~= 2
%             subtractThetaEst_subsetCh(FileNameBase,'ChListRef',1,'_Ref');
%         end
%         strParam = 'thetaEst_Ref';
%     else
%         strParam = [];
%     end
    MakeTheta_mpfi(FileNameBase,SampleRate,lfpSampleRate,xml.ElecGp,nChannels);
   
    strExt = [];
    % get theta phase
    fprintf('\nGet theta phase and track information');
    getTheta_mpfi(FileNameBase,lfpSampleRate,strExt); 
   
    % get theta phase and Track information for Spikes 
    fprintf('\nGet theta phase and track information for spikes');
    getSpikeInfo_mpfi(FileNameBase,lfpSampleRate,strExt);
    
    % get Stim information
    fprintf('\nGet stim information');
    getStim_mpfi(FileNameBase, lfpSampleRate);
    
    fprintf('\nxml saved into the output structures file: %s\n',...
            [FileNameBase '_BehavElectrDataLFP.mat']);
    save([FileNameBase '_BehavElectrDataLFP.mat'], ...
        'xml', '-append');
    
    clearvars -except recordFolder path
    
end
