function ProcessingAlignedInclNonRun(path,fileName,mazeSess)
% no cue passive --- condition 1
% middle cue passive --- condition 2
% 0.5s start cue passive --- condition 3
% 0.5s start cue active --- condition 4 

    GlobalConst;
    
    onlyRun = 0;
    spaceBin = 20; % mm
    corrIntervalT = 20; % sec
    corrIntervalD = 1800; % mm
    tc = 500; % ms
    cost = 2/(tc/1000*sampleFq);
    methodTheta = 1;

    disp('Get spikes per neuron per trial')
    getSpikesPerNeuPerTr(path,fileName,onlyRun,mazeSess);
    
    disp('Calculate running speed over time, aligned to run onset, include non-run')
    %% plotSpikeThetaFreqAmpRunOnsetMeanStd() also calculated the speed over time for good non-stim trials
    RunSpeedAlignedNonRun(path, fileName, onlyRun, mazeSess);

%     if(cond == 2)
%         disp('Plot spike rasters comparing good and bad trials (for non-stimulated trials), (Run, CueOff)')
%         plotSpikeRaster_GoodVsBadBeh_RunVsCueOff(path, fileName, onlyRun, mazeSess);  
%     else
%         disp('Plot spike rasters comparing good and bad trials (for non-stimulated trials), (Run, Rew, CueOnset)')
%         plotSpikeRaster_GoodVsBadBeh(path, fileName, onlyRun, mazeSess);  
%     end
%     
    % phase estimation for run onset
    disp(['Convolve spike trains with Gaussian kernel. only for spikes aligned' ...
        'to running onset, included spikes before run onset']);
    ConvSpikeTrain_AlignedRunOnset(path,fileName,onlyRun,mazeSess);
    
    disp(['Calculate peak firing rate with Gaussian kernel, only for spikes aligned' ...
        'to running onset, included spikes before run onset'])
    PeakFiringRate_AlignedRunOnset(path,fileName,1,onlyRun, mazeSess,0);
    
%     disp('Plot spike phase for non-stimulated trials aligned to run onset')
%     plotSpikePhaseRunOnset(path, fileName, onlyRun, mazeSess);
%     
    disp(['Calculate theta phase change over a trial after aligned to run onset'...
         'included spikes before run onset'])
    ThetaPhaseP2PRunOnset(path,fileName,methodTheta,onlyRun,mazeSess);

    disp('Calculate theta phase change over a trial after aligned to run onset, linear, in segments before running starts,linear')
    ThetaPhaseP2PAlignRunSegBefRun(path,fileName,methodTheta,onlyRun,mazeSess);
%     
% %     disp('Plot spike phase vs. theta cycle for non-stimulated trials aligned to run onset')
% %     plotSpikePhaseVsCRunOnset(path, fileName, onlyRun, mazeSess, methodTheta);
% 
%     disp('Plot theta frequency and amplitude for non-stimulated trials aligned to run onset')
%     plotSpikePhaseFreqRunOnset(path, fileName, onlyRun, mazeSess);
%     
%     disp('Plot mean and std of theta frequency and amplitude for non-stimulated good trials aligned to run onset')
%     plotSpikeThetaFreqAmpRunOnsetMeanStd(path, fileName, onlyRun, mazeSess);
% 

end
