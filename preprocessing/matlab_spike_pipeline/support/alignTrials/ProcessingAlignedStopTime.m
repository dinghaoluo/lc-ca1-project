function ProcessingAlignedStopTime(path,fileName,mazeSess,onlyRun)
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
 
    ThetaA
    
    disp(['Calculate peak firing rate with Gaussian kernel, only for spikes aligned' ...
        'to running onset'])
    PeakFiringRate_AlignedStopTime(path,fileName,1,onlyRun,mazeSess,0);

    disp(['Calculate peak firing rate with Gaussian kernel, only for spikes aligned' ...
        'to running onset, included spikes before run onset'])
    PeakFiringRate_AlignedRunOnset(path,fileName,1,onlyRun, mazeSess,0);
    
end
