function ProcessingMice_smTr2P(path, fileName, onlyRun, mazeSess)
% process the data after getting the DataStructure file
%
% e.g.: ProcessingMice_smTr('./','xzvr_PR1-20170804-01_DataStructure_mazeSection1_TrialType1',1)
%
% by Yingxue: 08/24/2017
    
    if(nargin == 2)
        onlyRun = 1;
        mazeSess = 1;
    elseif(nargin ==3)
        mazeSess = 1;
    end

    GlobalConst2P;
    
    % get the basic information
    disp('Get basic information');
    BasicInfo_smTr2P(path,fileName);
    
    % get spikes during the running period
    disp('Get spikes during the running period');
    SpikeDuringRun2P(path,fileName);
    
    % get spikes during the immobile period
    disp('Get spikes during the immobile period');
    SpikeImmobile2P(path,fileName);
        
    % calculate running speed and extract licking information
    disp('Calculate running speed and extract licking information');
    RunSpeed2P(path,fileName,0);
     
    % calculate mean firing rate
    disp('Calculate mean firing rate');
    MeanFiringRate2P(path,fileName,onlyRun);
   
    % smooth the spike train over time
    disp('Smooth spike trains over time');
    ConvSpikeTrainTimePar_smTr2P(path,fileName,timeBin,onlyRun,1);
         
    % smooth the spike trains over dist
    disp('Smooth spike trains over dist');
    ConvSpikeTrainDistPar_smTr2P(path,fileName,spaceBin,onlyRun,1);
    
    % calculate peak firing rate
    disp('Calculate peak firing rate');
    PeakFiringRate_smTr2P(path,fileName,spaceBin,1,onlyRun,0);
    
    % calculate CCG for each subsession
%     disp('Calculate CCG');
%     CalCCG2P(path,fileName,onlyRun);

    %  calculate spatial information
    disp('Calculate spatial information');
    GetFRMapInfo2P(path,fileName,onlyRun);
   
    disp('Single neuron correlation distance')
    neuronSpikeCorrDist_smTr2P(path,fileName,onlyRun,mazeSess,spaceBin,corrIntervalD);
    meanNeuronSpikeCorrDist_smTr2P(path,fileName,onlyRun,mazeSess,corrIntervalD);
    FieldDetection_smTr2P(path,fileName,onlyRun,2,spaceBin,mazeSess,corrIntervalD);
    
    close all;
    disp('FIN')

end
