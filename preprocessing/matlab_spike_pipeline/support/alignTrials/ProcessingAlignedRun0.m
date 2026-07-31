function ProcessingAlignedRun0(path,fileName,mazeSess,cond)
%% process the aligned data including all the spikes, including the stopping period

% no cue passive --- condition 1
% middle cue passive --- condition 2
% 0.5s start cue passive --- condition 3
% 0.5s start cue active --- condition 4 

    GlobalConst;
    
    onlyRun = 0;
    
    spaceBin = 20; % mm
    intervalTSpInfo = 10; % sec
    corrIntervalT = 20; % sec
    corrIntervalTMin = -10;
    corrIntervalD = 1800; % mm
    tc = 500; % ms
    cost = 2/(tc/1000*sampleFq);
    methodTheta = 1;

    disp('Get spikes per neuron per trial')
    getSpikesPerNeuPerTr(path,fileName,onlyRun,mazeSess);
% 
    % smooth spike trains
    disp('Convolve spike trains with Gaussian kernel (Run,Rew,CueOn)')
    ConvSpikeTrain_Aligned1(path, fileName,onlyRun,mazeSess);
    
    disp(['Calculate peak firing rate with Gaussian kernel, only for spikes aligned' ...
        'to running onset'])
    PeakFiringRate_Aligned(path,fileName,1,onlyRun,mazeSess,0);

    %% lfp theta phase over time
    disp('lfp theta phase over time')
    thetaPhaseOverTimeAligned(path,fileName,mazeSess,onlyRun);
    
    disp('Calculate lick over distance')
    LickOverDistAligned(path, fileName, onlyRun, mazeSess);
    
    disp('Calculate running speed over distance')
    RunSpeedOverDistAligned(path, fileName, onlyRun, mazeSess);
    
    disp('Calculate mean firing rate')
    MeanFiringRateAligned(path,fileName,onlyRun,mazeSess,0);
    
    disp('Calculate spike theta phase change over a trial after aligned to run onset, linear')
    methodTheta = 1; 
    ThetaPhaseP2PAlignRun(path,fileName,methodTheta,onlyRun,mazeSess);

    methodTheta = 0; 
    disp('Calculate spike theta phase change over a trial after aligned to run onset, Hilbert')
    ThetaPhaseP2PAlignRun(path,fileName,methodTheta,onlyRun,mazeSess);

    %% changed by Yingxue on 3/4/2022
    %% reward onset aligned
    disp(['Convolve spike trains with Gaussian kernel. only for spikes aligned' ...
        'to reward onset, included spikes before reward']);
    ConvSpikeTrain_AlignedOnset(path,fileName,onlyRun,mazeSess,2);
    
    disp(['Calculate peak firing rate with Gaussian kernel, only for spikes aligned' ...
        'to reward onset, included spikes before reward'])
    PeakFiringRate_AlignedOnset(path,fileName,1,onlyRun, mazeSess,2,0);
    
    %% run onset aligned
    disp(['Convolve spike trains with Gaussian kernel. only for spikes aligned' ...
        'to running onset, included spikes before run onset']);
    ConvSpikeTrain_AlignedOnset(path,fileName,onlyRun,mazeSess,1);
    
    disp(['Calculate peak firing rate with Gaussian kernel, only for spikes aligned' ...
        'to running onset, included spikes before run onset'])
    PeakFiringRate_AlignedOnset(path,fileName,1,onlyRun, mazeSess,1,0);
    %%
    
    if(cond == 2)
        disp('Plot spike rasters comparing good and bad trials (for non-stimulated trials), (Run, CueOff)')
        plotSpikeRaster_GoodVsBadBeh_RunVsCueOff(path, fileName, onlyRun, mazeSess);  
    else
        disp('Plot spike rasters comparing good and bad trials (for non-stimulated trials), (Run, Rew, CueOnset)')
        plotSpikeRaster_GoodVsBadBeh(path, fileName, onlyRun, mazeSess, cond);  
        % cond = 1: only the align-to-run plots are extended to before run
        % cond = 0: both the align-to-run and align-to-rew plots are extended to
        % before the align onset
        % added "cond" on 3/4/2022
    end
    
    close all;
    
    disp('Detect neurons with significant peak')
    detectNeuWithPeak(path,fileName,onlyRun, mazeSess);

% %     disp('Plot spike phase for non-stimulated trials aligned to run onset')
% %     plotSpikePhaseRunOnset(path, fileName, onlyRun, mazeSess);
%     
    disp(['Calculate spike theta phase change over a trial after aligned to run onset'...
         'included spikes before run onset'])
    ThetaPhaseP2PRunOnset(path,fileName,methodTheta,onlyRun,mazeSess);
%     
% %     disp('Plot spike phase vs. theta cycle for non-stimulated trials aligned to run onset')
% %     plotSpikePhaseVsCRunOnset(path, fileName, onlyRun, mazeSess, methodTheta);
% 
% %     disp('Plot theta frequency and amplitude for non-stimulated trials aligned to run onset')
% %     plotSpikePhaseFreqRunOnset(path, fileName, onlyRun, mazeSess);
%     
    disp('Plot mean and std of theta frequency and amplitude for non-stimulated good trials aligned to run onset')
    plotSpikeThetaFreqAmpRunOnsetMeanStd(path, fileName, onlyRun, mazeSess);

end
