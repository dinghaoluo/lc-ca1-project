function ProcessingMice_smTr(path, fileName, onlyRun)
% process the data after getting the DataStructure file
%
% e.g.: ProcessingMice_smTr('./','xzvr_PR1-20170804-01_DataStructure_mazeSection1_TrialType1',1)
%
% by Yingxue: 08/24/2017
    
    if(nargin == 2)
        onlyRun = 1;
    end
    spaceBin = 20; % mm
    timeBin = 0.0096; % s
    methodTheta = 1;
    methodTheta = 0;    
    % get the basic information
    disp('Get basic information');
    BasicInfo_smTr(path,fileName);
    
    % get spikes during the running period
    disp('Get spikes during the running period');
    SpikeDuringRun(path,fileName);
    
    % get spikes during the immobile periodedit
    disp('Get spikes during the immobile period');
    SpikeImmobile(path,fileName);
        
    % calculate running speed and extract licking information
    disp('Calculate running speed and extract licking information');
    RunSpeed(path,fileName,0);
     
    % calculate mean firing rate
    disp('Calculate mean firing rate');
    MeanFiringRate(path,fileName,onlyRun);
    
    % get the depth of each neuron relative to the layer center
    % MAKE SURE THAT THE BAD CHANNELS ARE NOT DELETED DURING SPIKE SORTING
    disp('Get the relative depth of each neuron')
    GetNeuRelativeDepth(path,fileName,onlyRun);

    disp('Smooth spike trains over time');
    ConvSpikeTrainTimePar_smTr(path,fileName,timeBin,onlyRun,1);
         
    % smooth the spike trains
    disp('Smooth spike trains over dist');
    ConvSpikeTrainDistPar_smTr(path,fileName,spaceBin,onlyRun,1);
    
    % calculate peak firing rate
    disp('Calculate peak firing rate');
    PeakFiringRate_smTr(path,fileName,spaceBin,1,onlyRun,0);
    
    % calculate theta phase
    disp('Calculate theta phase Hilbert');
    methodTheta = 0;  
    ThetaPhaseP2P(path,fileName,methodTheta,onlyRun,1);
    
    disp('Calculate theta phase Linear');
    methodTheta = 1;
    ThetaPhaseP2P(path,fileName,methodTheta,onlyRun,1);

    % calculate CCG for each subsession
    disp('Concatenating spikes from each trials during run');
    ConcatenateSpikes_smTr(path,fileName,onlyRun);

    % calculate CCG for each subsession
    disp('Calculate CCG');
    CalCCG(path,fileName,onlyRun);

    % calculate theta modulation
    disp('Calculate theta modulation');
    ThetaModulation(path,fileName,onlyRun)
    
    % detect interneurons
    disp('Detect interneurons');
    DetectInt(path,fileName, onlyRun);
    
    % calculate spatial information
    disp('Calculate spatial information');
    GetFRMapInfo(path,fileName,onlyRun);

    % identify place fields
    disp('Identify place field');
    FieldWidthLR(path,fileName,spaceBin,methodTheta,2,0,onlyRun);
    close all;
    
    % extract bursting spikes
    disp('Extract bursting spikes Hilbert');
    methodTheta = 0; 
    BurstAll(path,fileName,methodTheta,0,onlyRun);
    
    disp('Extract bursting spikes Linear');
    methodTheta = 1;
    BurstAll(path,fileName,methodTheta,0,onlyRun);
    
%     % order trials based on behavior parameters
%     disp('Order trials based on behavior parameters');
%     OrderTrials(path,fileName,onlyRun);
%     
%     % calculate LFP spectrum
%     LFPSpectrumPerTrial(path,fileName,0); 
%     
%     % identify place fields for each session
%     fileNameInfo = [fileName '_Info.mat'];
%     fullPath = [path fileNameInfo];
%     load(fullPath,'beh');
%     disp('Field identification for each session');
%     for Sess = 1:length(beh.mazeSess)
%         if(beh.goodSess(Sess) == 1)
%             fprintf('\nSession %d\n',Sess);
%             FieldWidthSession(path,fileName,Sess,spaceBin,methodTheta,...
%                                 0,0,onlyRun);
% %             disp('Pause and check the detected field.');
% %             pause;
%         end
%     end  

    clearvars;
end
