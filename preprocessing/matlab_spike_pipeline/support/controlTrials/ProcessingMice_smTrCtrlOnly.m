function ProcessingMice_smTrCtrlOnly(path, fileName, onlyRun, mazeSess)
% process the data after getting the DataStructure file
%
% e.g.: ProcessingMice_smTr('./','xzvr_PR1-20170804-01_DataStructure_mazeSection1_TrialType1',1)
%
% by Yingxue: 08/24/2017
    
    if(nargin == 3)
        onlyRun = 1;
    end
    spaceBin = 20; % mm
    timeBin = 0.0096; % s
    methodTheta = 1;
    methodTheta = 0;    

    % calculate mean firing rate
    disp('Calculate mean firing rate for ctrl trials');
    MeanFiringRateCtrlOnly(path,fileName,onlyRun,mazeSess);
    
    % calculate peak firing rate
    disp('Calculate peak firing rate for ctrl trials');
    PeakFiringRate_smTrCtrlOnly(path,fileName,spaceBin,1,onlyRun,mazeSess);
    
    % calculate theta phase
    disp('Calculate theta phase Hilbert for ctrl trials');
    methodTheta = 0;  
    ThetaPhaseP2PCtrlOnly(path,fileName,methodTheta,onlyRun,mazeSess);
    
    disp('Calculate theta phase Linear for ctrl trials');
    methodTheta = 1;
    ThetaPhaseP2PCtrlOnly(path,fileName,methodTheta,onlyRun,mazeSess);
    
    % calculate CCG for each subsession
    disp('Concatenating spikes from each trials during run for ctrl trials');
    ConcatenateSpikes_smTr_CtrlOnly(path,fileName,onlyRun,mazeSess);
    
    % calculate CCG for each subsession
    disp('Calculate CCG for ctrl trials');
    CalCCG_CtrlOnly(path,fileName,onlyRun,mazeSess);
    
    % calculate theta modulation
    disp('Calculate theta modulation for ctrl trials');
    ThetaModulation_CtrlOnly(path,fileName,onlyRun,mazeSess);
    
    % calculate spatial information
    disp('Calculate spatial information for ctrl trials');
    GetFRMapInfo_CtrlOnly(path,fileName,onlyRun,mazeSess);

    % extract bursting spikes
    disp('Extract bursting spikes Hilbert for ctrl trials');
    methodTheta = 0; 
    BurstAll_CtrlOnly(path,fileName,methodTheta,onlyRun,mazeSess);
    
    disp('Extract bursting spikes Linear for ctrl trials');
    methodTheta = 1;
    BurstAll_CtrlOnly(path,fileName,methodTheta,onlyRun,mazeSess);
    
    disp('Single neuron correlation distance')
    neuronSpikeCorrDist_CtrlOnly(path,fileName,onlyRun,spaceBin,mazeSess);
    meanNeuronSpikeCorrDist_CtrlOnly(path,fileName,onlyRun,mazeSess);
    
    disp('Firing field detection based on correlation distance and spatial information')
    FieldDetection_CtrlOnly(path,fileName,onlyRun,spaceBin,0,mazeSess); % mazeSess is only used for plotting purpose
    
    clearvars;
end
