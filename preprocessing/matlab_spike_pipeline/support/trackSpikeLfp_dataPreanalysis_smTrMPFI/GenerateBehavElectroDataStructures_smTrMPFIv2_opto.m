function [Spike Clu Track Laps xml] = GenerateBehavElectroDataStructures_smallTMPFIv2_opto(path, FileNameBaseList, varargin)

% function [Data Track Laps Spike xml] = GenerateBehavElectroDataStructures_MPFI(FileNameBase);
%
% This function reads the *.lfp  *.clu  *.res and *.whl files and
% produces data-structures with different sorting of behavioral and
% electrophysiol data:
% ouput structures:
%
%
% e.g.: GenerateBehavElectroDataStructures_smTrMPFIv1
%  ('Z:\Raphael_tests\mice_expdata\ANM001\A001-20180928', 'A001-20180928-01')
% 

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
    
    % in case this is a pixel rig concatenated recording
    if exist([FileNameBase '-concat.dat'], 'file')==2
        [xml] = LoadXml_e([FileNameBase '-concat.xml']);
    end
    
    % get recording file information from the xml file
    fileinfo = dir([FileNameBaseList '.dat']);
    
    % in case this is a pixel rig concatenated recording
    if exist([FileNameBaseList '-concat.dat'], 'file')==2
        fileinfo = dir([FileNameBaseList '-concat.dat']);
    end
    
    nChannels = xml.nChannels;
        % total number of recording channels    
    SampleRate = xml.SampleRate;
        % recording sampling rate (default: 20,000 Hz)
    lfpSampleRate = xml.lfpSampleRate;
        % local field potential sampling rate (default: 1250 Hz)
    
%     load behavior file
    fprintf(['\nLoad behavior file and align the event time between',...
            'behavior and recording\n']);
    check = Arduino2RecordT_smTr_opto(FileNameBase,SampleRate,lfpSampleRate,nChannels);
    if(check == 1)
        return;
    end
    
    % load stimulation file
%     fprintf(['\nLoad stimulation file and align the event time between'...
%             'stimulation and recording\n']);
%     Stim2RecordT_mpfi(FileNameBase,SampleRate,lfpSampleRate,nChannels);
    
    % align tracking data with .dat file
    fprintf('\nAlign tracking data with .dat file...\n')
    aligntsp2dat_smTr(FileNameBase,SampleRate,lfpSampleRate,nChannels);
        
    % get tracking related data
    fprintf('\nGet Track and Laps...\n')
    SortTrials_3armMaze_smTr(FileNameBase, lfpSampleRate);
    
    % get spikes related data
    fprintf('\nGet Spike and Clu....\n')
    getSpikes_mpfi(FileNameBase, SampleRate, lfpSampleRate); 
    
    % calculate CCG
    fprintf('\nCalculate CCG\n');
    calACG_mpfi(FileNameBase,SampleRate);
    
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
    fprintf('\nGet tagging stim information');
    getStim_mpfi_smTr(FileNameBase, lfpSampleRate);
    
    fprintf('\nxml saved into the output structures file: %s\n',...
            [FileNameBase '_BehavElectrDataLFP.mat']);
    save([FileNameBase '_BehavElectrDataLFP.mat'], ...
        'xml', '-append');
    
    clearvars -except recordFolder path
    
end
