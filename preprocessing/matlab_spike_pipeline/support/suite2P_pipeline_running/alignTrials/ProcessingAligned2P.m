function ProcessingAligned2P(path,fileName,onlyRun,mazeSess)
% no cue passive --- condition 1
% middle cue passive --- condition 2
% 0.5s start cue passive --- condition 3
% 0.5s start cue active --- condition 4 
%
% ProcessingAligned2P('./','A582-20211111-02_DataStructure_mazeSection1_TrialType1',1,1)

    GlobalConst2P;
    
    %% align the spikes based on different run onset/reward/cue onset
    disp('Align trials based on run onset')
    alignToRunOnset2P1(path, fileName,mazeSess);
   
    disp('Align trials based on reward onset')
    alignToReward2P(path, fileName,mazeSess);
    
    disp('Align trials based on cue onset')
    alignToCue2P(path, fileName,mazeSess);

    disp('Align trials based on cue offset')
    alignToCueOff2P(path, fileName,mazeSess);
    
    %% get behavior parameters
    disp('Get behavior parameters')
    getBehParameters2P(path,fileName,mazeSess);
   
    % smooth spike trains
    disp('Convolve spike trains with Gaussian kernel (Run,Rew,CueOn)')
    ConvSpikeTrain_Aligned2P(path, fileName,onlyRun,mazeSess);
    
    disp(['Calculate peak firing rate with Gaussian kernel, only for spikes aligned' ...
        'to running onset'])
    PeakFiringRate_Aligned2P(path,fileName,1,onlyRun,mazeSess,0);

    disp('Calculate lick over distance')
    LickOverDistAligned2P(path, fileName, onlyRun, mazeSess);
    
    disp('Calculate running speed over distance')
    RunSpeedOverDistAligned2P(path, fileName, onlyRun, mazeSess);
     
    disp('Calculate mean firing rate')
    MeanFiringRateAligned2P(path,fileName,onlyRun,mazeSess,0);
     
    %% phase estimation for run onset
    disp(['Convolve spike trains with Gaussian kernel. only for spikes aligned' ...
        'to running onset, included spikes before run onset']);
    ConvSpikeTrain_AlignedRunOnset2P(path,fileName,onlyRun,mazeSess);
    
    disp(['Calculate peak firing rate with Gaussian kernel, only for spikes aligned' ...
        'to running onset, included spikes before run onset'])
    PeakFiringRate_AlignedRunOnset2P(path,fileName,1,onlyRun, mazeSess,0);

    disp(['Convolve spike trains with Gaussian kernel. only for spikes aligned' ...
        'to reward, included spikes before reward']);
    ConvSpikeTrain_AlignedReward2P(path,fileName,onlyRun,mazeSess); %% added on 9/7/2023

    disp('Get spikes per neuron per trial per distance bin')
    getTimePerDistBin2P(path,fileName,onlyRun,mazeSess);
    disp('Convolve spike trains with Gaussian kernel (distance)')
    ConvSpikeTrainDist_Aligned2P(path, fileName,spaceBin,onlyRun,mazeSess);
    
    delete(gcp('nocreate'));
    parpool('local',4);
    
    % single neuron corr between trials
    disp('Single neuron correlation T')
    neuronSpikeCorrT2P(path,fileName,onlyRun,mazeSess,corrIntervalT);
    meanNeuronSpikeCorrT2P(path,fileName,onlyRun,mazeSess,corrIntervalT);
    
    disp('Single neuron correlation distance')
    neuronSpikeCorrDist2P(path,fileName,onlyRun,mazeSess,corrIntervalD);
    meanNeuronSpikeCorrDist2P(path,fileName,onlyRun,mazeSess,corrIntervalD);
    
    %% single neuron trial similarity
    disp('Single neuron cosine similarity')
    spikeTrainSimilarityT2P(path,fileName,onlyRun,mazeSess,corrIntervalT);    
    meanSpikeTrainSimilarityT2P(path,fileName,onlyRun,mazeSess,corrIntervalT);

    %% population corrT between trials
    disp('Population vector correlation time')
    popSpikeCorrT2P(path,fileName,onlyRun,mazeSess,corrIntervalT);
    meanPopSpikeCorrT2P(path,fileName,onlyRun,mazeSess,corrIntervalT);
   
    %% population corrDist between trials
    disp('Population vector correlation distance')
    popSimilarityT2P(path,fileName,onlyRun,mazeSess,corrIntervalT);
    meanPopSpikeSimT2P(path,fileName,onlyRun,mazeSess,corrIntervalT);

    % calculate spatial information
    disp('Calculate spatial information');
    GetFRMapInfo_Aligned2P(path,fileName,mazeSess,onlyRun,intervalTSpInfo);

    disp('Calculate spatial information (aligned to cue)');
    GetFRMapInfo_AlignedCue2P(path,fileName,mazeSess,onlyRun,intervalTSpInfo);
    
    disp('Firing field detection based on correlation distance and spatial information')
    neuronSpikeCorrT2P(path,fileName,onlyRun,mazeSess,intervalTSpInfo);
    meanNeuronSpikeCorrT2P(path,fileName,onlyRun,mazeSess,intervalTSpInfo);
    FieldDetectionAligned2P(path,fileName,onlyRun,2,mazeSess,intervalTSpInfo); % mazeSess is only used for plotting purpose
    FieldDetectionAligned2PAllTrials(path,fileName,onlyRun,2,mazeSess,intervalTSpInfo); % detecting fields using all trials
    
    delete(gcp('nocreate'));
    
    begT = 3;
    for i = 1:2
        plotdFFAligned2P(path,fileName,mazeSess,i,intervalTSpInfo,begT,onlyRun); % plot dFF for each neuron
        close all;
    end
    
end
