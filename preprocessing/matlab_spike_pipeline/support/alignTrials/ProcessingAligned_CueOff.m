function ProcessingAligned_CueOff(path,fileName,onlyRun,mazeSess)

    GlobalConst;
    
    spaceBin = 20; % mm
    corrIntervalT = 20; % sec
    corrIntervalD = 1800; % mm
    tc = 500; % ms
    cost = 2/(tc/1000*sampleFq);
    methodTheta = 1;
    
%     % phase estimation for cue offset 
%     disp('Convolve spike trains with Gaussian kernel (CueOff)')
%     ConvSpikeTrain_AlignedCueOff(path, fileName,onlyRun,mazeSess);
% 
%     disp(['Convolve spike trains with Gaussian kernel. only for spikes aligned' ...
%         'to cue offset, included spikes before cue offset']);
%     ConvSpikeTrain_AlignedCueOffWithBef(path,fileName,onlyRun,mazeSess);
%     
%     disp(['Calculate peak firing rate with Gaussian kernel, only for spikes aligned' ...
%         'to cue offset, included spikes before cue offset']);
%     PeakFiringRate_AlignedCueOff(path,fileName,1,onlyRun,mazeSess,0);
%     
% %     disp('Plot spike phase for non-stimulated trials aligned to cue offset')
% %     plotSpikePhaseCueOffset(path, fileName, onlyRun, mazeSess);
%     
    disp('Calculate theta phase change over a trial after aligned to cue offset')
    ThetaPhaseP2PCueOffset(path,fileName,methodTheta,onlyRun,mazeSess);
    
%     disp('Plot spike phase vs. theta cycle for non-stimulated trials aligned to cue offset')
%     plotSpikePhaseVsCCueOffset(path, fileName, onlyRun, mazeSess, methodTheta);
%     
%     disp('Plot spike phase and theta frequency for non-stimulated trials aligned to cue offset')
%     plotSpikePhaseFreqCueOffset(path, fileName, onlyRun, mazeSess);
    
end
