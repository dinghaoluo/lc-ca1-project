function ProcessingMice_smTr_GoodTr(path, fileName, onlyRun, mazeSess)
% process the data after running ProcessingMice_smTr, separating good and
% bad trials (non-stimulated ctrl) during this set of analysis
%
% e.g.: ProcessingMice_smTr_GoodTr('./','xzvr_PR1-20170804-01_DataStructure_mazeSection1_TrialType1',1)
%
% by Yingxue: 12/18/2020

    if(nargin == 2)
        onlyRun = 1;
    end
    spaceBin = 20; % mm
    timeBin = 0.0096; % s
    methodTheta = 1;
    methodTheta = 0;    

    % calculate mean firing rate
    disp('Calculate mean firing rate for good and bad trials');
    MeanFiringRateGoodTr(path,fileName,onlyRun);
    
    % calculate peak firing rate
    disp('Calculate peak firing rate for good and bad trials');
    PeakFiringRate_smTr_GoodTr(path,fileName,spaceBin,1,onlyRun);
    
    % calculate theta phase
    disp('Calculate theta phase Hilbert for good and bad trials');
    methodTheta = 0;  
    ThetaPhaseP2P_GoodTr(path,fileName,methodTheta,onlyRun);
    
    disp('Calculate theta phase Linear for good and bad trials');
    methodTheta = 1;
    ThetaPhaseP2P_GoodTr(path,fileName,methodTheta,onlyRun);
    
    % calculate CCG for each subsession
    disp('Concatenating spikes from each trials during run for good and bad trials');
    ConcatenateSpikes_smTr_GoodTr(path,fileName,onlyRun);
    
    % calculate CCG for each subsession
    disp('Calculate CCG for good and bad trials');
    CalCCG_GoodTr(path,fileName,onlyRun);
    
    % calculate theta modulation
    disp('Calculate theta modulation for good and bad trials');
    ThetaModulation_GoodTr(path,fileName,onlyRun);
    
    % calculate spatial information
    disp('Calculate spatial information for good and bad trials');
    GetFRMapInfo_GoodTr(path,fileName,onlyRun);

    % extract bursting spikes
    disp('Extract bursting spikes Hilbert for good and bad trials');
    methodTheta = 0; 
    BurstAll_GoodTr(path,fileName,methodTheta,onlyRun);
    
    disp('Extract bursting spikes Linear for good and bad trials');
    methodTheta = 1;
    BurstAll_GoodTr(path,fileName,methodTheta,onlyRun);
    
    disp('Single neuron correlation distance')
    neuronSpikeCorrDist_GoodTr(path,fileName,onlyRun,spaceBin);
    meanNeuronSpikeCorrDist_GoodTr(path,fileName,onlyRun);
    
    disp('Firing field detection based on correlation distance and spatial information')
    FieldDetection_GoodTr(path,fileName,onlyRun,spaceBin,0,mazeSess); % mazeSess is only used for plotting purpose
    
    clearvars;
