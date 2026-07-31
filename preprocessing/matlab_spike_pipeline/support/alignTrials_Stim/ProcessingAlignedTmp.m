function ProcessingAligned(path,fileName,onlyRun,mazeSess,cond)
% no cue passive --- condition 1
% middle cue passive --- condition 2
% 0.5s start cue passive --- condition 3
% 0.5s start cue active --- condition 4 

    GlobalConst;
    
    spaceBin = 20; % mm
    intervalTSpInfo = 10; % sec
    corrIntervalT = 20; % sec
    corrIntervalTMin = -10;
    corrIntervalD = 1800; % mm
    tc = 500; % ms
    cost = 2/(tc/1000*sampleFq);
    methodTheta = 1;

%     disp('Calculate LFP theta power') 
%     ThetaPower_Hilbert(path,fileName);
%     
%     % align the spikes based on different run onset/reward/cue onset
%     disp('Align trials based on run onset')
%     alignToRunOnset(path, fileName,mazeSess);
%    
%     disp('Align trials based on reward onset')
%     alignToReward(path, fileName,mazeSess);
%     
%     disp('Align trials based on cue onset')
%     alignToCue(path, fileName,mazeSess);
% 
%     disp('Align trials based on cue offset')
%     alignToCueOff(path, fileName,mazeSess);
%     
%     % get behavior parameters
%     disp('Get behavior parameters')
%     getBehParameters(path,fileName,mazeSess);
%     
%     disp('Get spikes per neuron per trial')
%     getSpikesPerNeuPerTr(path,fileName,onlyRun,mazeSess);
% 
%     % lfp theta phase over time
%     disp('lfp theta phase over time')
%     thetaPhaseOverTimeAligned(path,fileName,mazeSess);
    
%     
%     if(cond == 2)
%         disp('Plot spike rasters comparing good and bad trials (for non-stimulated trials), (Run, CueOff)')
%         plotSpikeRaster_GoodVsBadBeh_RunVsCueOff(path, fileName, onlyRun, mazeSess);  
%     else
%         disp('Plot spike rasters comparing good and bad trials (for non-stimulated trials), (Run, Rew, CueOnset)')
%         plotSpikeRaster_GoodVsBadBeh(path, fileName, onlyRun, mazeSess);  
%     end
%      
%     %% smooth spike trains
%     disp('Convolve spike trains with Gaussian kernel (Run,Rew,CueOn)')
%     ConvSpikeTrain_Aligned(path, fileName,onlyRun,mazeSess);
%     
%     disp(['Calculate peak firing rate with Gaussian kernel, only for spikes aligned' ...
%         'to running onset'])
%     PeakFiringRate_Aligned(path,fileName,1,onlyRun,mazeSess,0);
% 
%     disp('Calculate lick over distance')
%     LickOverDistAligned(path, fileName, onlyRun, mazeSess);
%     
%     disp('Calculate running speed over distance')
%     RunSpeedOverDistAligned(path, fileName, onlyRun, mazeSess);
%     
%     disp('Calculate mean firing rate')
%     MeanFiringRateAligned(path,fileName,onlyRun,mazeSess,0);
%     
%     disp('Calculate spike theta phase change over a trial after aligned to run onset, linear')
%     methodTheta = 1; 
%     ThetaPhaseP2PAlignRun(path,fileName,methodTheta,onlyRun,mazeSess);
% 
%     disp('Calculate spike theta phase change over a trial after aligned to run onset, linear, in segments after running starts')
%     ThetaPhaseP2PAlignRunSeg(path,fileName,methodTheta,onlyRun,mazeSess);
% 
%     methodTheta = 0; 
%     disp('Calculate spike theta phase change over a trial after aligned to run onset, Hilbert')
%     ThetaPhaseP2PAlignRun(path,fileName,methodTheta,onlyRun,mazeSess);
% 
%     disp('Calculate spike theta phase change over a trial after aligned to run onset, Hilbert, in segments after running starts')
%     ThetaPhaseP2PAlignRunSeg(path,fileName,methodTheta,onlyRun,mazeSess);
% 
%     % phase estimation for run onset
%     disp(['Convolve spike trains with Gaussian kernel. only for spikes aligned' ...
%         'to running onset, included spikes before run onset']);
%     ConvSpikeTrain_AlignedRunOnset(path,fileName,onlyRun,mazeSess);
%     
%     disp(['Calculate peak firing rate with Gaussian kernel, only for spikes aligned' ...
%         'to running onset, included spikes before run onset'])
%     PeakFiringRate_AlignedRunOnset(path,fileName,1,onlyRun, mazeSess,0);
% 
% %     disp('Plot spike phase for non-stimulated trials aligned to run onset')
% %     plotSpikePhaseRunOnset(path, fileName, onlyRun, mazeSess);
%     
%     disp(['Calculate spike theta phase change over a trial after aligned to run onset'...
%          'included spikes before run onset'])
%     ThetaPhaseP2PRunOnset(path,fileName,methodTheta,onlyRun,mazeSess);
% %     
% % %     disp('Plot spike phase vs. theta cycle for non-stimulated trials aligned to run onset')
% % %     plotSpikePhaseVsCRunOnset(path, fileName, onlyRun, mazeSess, methodTheta);
% % 
%     disp('Plot theta frequency and amplitude for non-stimulated trials aligned to run onset')
%     plotSpikePhaseFreqRunOnset(path, fileName, onlyRun, mazeSess);
%     
%     disp('Plot mean and std of theta frequency and amplitude for non-stimulated good trials aligned to run onset')
%     plotSpikeThetaFreqAmpRunOnsetMeanStd(path, fileName, onlyRun, mazeSess);
% 
%     disp('Concatenating spikes from each trials during run');
%     ConcatenateSpikesAlignedRun(path,fileName,onlyRun,mazeSess);
% 
%     % calculate CCG for each subsession
%     disp('Calculate CCG');
%     CalCCGAlignedRun(path,fileName,onlyRun,mazeSess);
% 
%     % calculate theta modulation
%     disp('Calculate theta modulation');
%     ThetaModulationAlignedRun(path,fileName,onlyRun,mazeSess);
%         
%     disp('Extract bursting spikes Hilbert');
%     methodTheta = 0; 
%     BurstAllAlignedRun(path,fileName,methodTheta,0,onlyRun,mazeSess);
%     
%     disp('Extract bursting spikes Linear');
%     methodTheta = 1;
%     BurstAllAlignedRun(path,fileName,methodTheta,0,onlyRun,mazeSess);
    
%     % don't need to run AutoCorrSession because the same has been done for the whole
%     % recording session in DetectInt
%     disp('Extract autocorrelation parameters')
%     AutoCorrSession(path,fileName, onlyRun, mazeSess);

%     % calculate spatial information
%     disp('Calculate spatial information');
%     GetFRMapInfo_Aligned(path,fileName,onlyRun,mazeSess,intervalTSpInfo);
    
%     disp('Get spikes per neuron per trial per distance bin')
%     getTimePerDistBin(path,fileName,onlyRun,mazeSess);
%     disp('Convolve spike trains with Gaussian kernel (distance)')
%     ConvSpikeTrainDist_Aligned(path, fileName,spaceBin,onlyRun,mazeSess);
%     
%     %% single neuron corr between trials
%     disp('Single neuron correlation T')
%     neuronSpikeCorrT(path,fileName,onlyRun,mazeSess,corrIntervalT,corrIntervalTMin);
%     meanNeuronSpikeCorrT(path,fileName,onlyRun,mazeSess,corrIntervalT);
%     
%     disp('Single neuron correlation distance')
%     neuronSpikeCorrDist(path,fileName,onlyRun,mazeSess,corrIntervalD);
%     meanNeuronSpikeCorrDist(path,fileName,onlyRun,mazeSess,corrIntervalD);
%      
%     %% single neuron trial similarity
%     disp('Single neuron cosine similarity')
%     spikeTrainSimilarityT(path,fileName,onlyRun,mazeSess,corrIntervalT,corrIntervalTMin);    
%     meanSpikeTrainSimilarityT(path,fileName,onlyRun,mazeSess,corrIntervalT);
%     
%     %% single neuron Victor & Purpura spike time distance
%     disp('Single neuron similarity Victor & Purpura (spike time and interval)')
%     spikeTrainSimilarityVP(path,fileName,onlyRun,mazeSess,cost,corrIntervalT,corrIntervalTMin);
%     meanSpikeTrainSimilarityVP(path,fileName,onlyRun,mazeSess,cost,corrIntervalT);
%     %% single neuron Victor & Purpura spike time interval
%     meanSpikeTrainSimilarityVPI(path,fileName,onlyRun,mazeSess,cost,corrIntervalT);
%     
%     %% single neuron Victor & Purpura spike time interval
%     disp('Single neuron similarity Van Rossum')
%     spikeTrainSimilarityVanRossum(path,fileName,onlyRun,mazeSess,tc,corrIntervalT,corrIntervalTMin);
%     meanSpikeTrainSimilarityVanRossum(path,fileName,onlyRun,mazeSess,tc,corrIntervalT);
%     
%     % population corrT between trials
%     disp('Population vector correlation time')
%     popSpikeCorrT(path,fileName,onlyRun,mazeSess,corrIntervalT,corrIntervalTMin);
%     meanPopSpikeCorrT(path,fileName,onlyRun,mazeSess,corrIntervalT);
%     
%     %% population corrDist between trials
%     disp('Population vector correlation distance')
%     popSimilarityT(path,fileName,onlyRun,mazeSess,corrIntervalT,corrIntervalTMin);
%     meanPopSpikeSimT(path,fileName,onlyRun,mazeSess,corrIntervalT);

    disp('Firing field detection based on correlation distance and spatial information')
    neuronSpikeCorrT(path,fileName,onlyRun,mazeSess,intervalTSpInfo,corrIntervalTMin);
    meanNeuronSpikeCorrT(path,fileName,onlyRun,mazeSess,intervalTSpInfo);
%     FieldDetectionAligned(path,fileName,onlyRun,2,mazeSess,intervalTSpInfo); % mazeSess is only used for plotting purpose

end
