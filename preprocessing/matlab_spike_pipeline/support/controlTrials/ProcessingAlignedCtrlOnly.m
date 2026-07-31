function ProcessingAlignedCtrlOnly(path,fileName,onlyRun,mazeSess)

    GlobalConst;
    
    spaceBin = 20; % mm
    intervalTSpInfo = 10; % sec
    corrIntervalT = 20; % sec
    corrIntervalTMin = -10;
    corrIntervalD = 1800; % mm
    tc = 500; % ms
    cost = 2/(tc/1000*sampleFq);
    methodTheta = 1;

    disp(['Calculate peak firing rate with Gaussian kernel, only for spikes aligned' ...
        'to running onset (Ctrl)'])
    PeakFiringRate_AlignedCtrlOnly(path,fileName,1,onlyRun,mazeSess,0);

    %% lfp theta phase over time
    disp('lfp theta phase over time (Ctrl)')
    thetaPhaseOverTimeAlignedCtrlOnly(path,fileName,mazeSess,onlyRun);
    
    disp('Calculate lick over distance (Ctrl)')
    LickOverDistAlignedCtrlOnly(path, fileName, onlyRun, mazeSess);
    
    disp('Calculate running speed over distance (Ctrl)')
    RunSpeedOverDistAlignedCtrlOnly(path, fileName, onlyRun, mazeSess);
    
    disp('Calculate mean firing rate (Ctrl)')
    MeanFiringRateAlignedCtrlOnly(path,fileName,onlyRun,mazeSess,0);
    
    disp('Calculate spike theta phase change over a trial after aligned to run onset, linear (Ctrl)')
    methodTheta = 1; 
    ThetaPhaseP2PAlignRunCtrlOnly(path,fileName,methodTheta,onlyRun,mazeSess);

    methodTheta = 0; 
    disp('Calculate spike theta phase change over a trial after aligned to run onset, Hilbert (Ctrl)')
    ThetaPhaseP2PAlignRunCtrlOnly(path,fileName,methodTheta,onlyRun,mazeSess);

    disp(['Calculate peak firing rate with Gaussian kernel, only for spikes aligned' ...
        'to running onset, included spikes before run onset (Ctrl)'])
    PeakFiringRate_AlignedRunOnsetCtrlOnly(path,fileName,1,onlyRun, mazeSess,0);

    disp(['Calculate spike theta phase change over a trial after aligned to run onset'...
         'included spikes before run onset (Ctrl)'])
    ThetaPhaseP2PRunOnsetCtrlOnly(path,fileName,methodTheta,onlyRun,mazeSess);
   
    disp('Concatenating spikes from each trials during run (Ctrl)');
    ConcatenateSpikesAlignedRunCtrlOnly(path,fileName,onlyRun,mazeSess);

    %% calculate CCG for each subsession
    disp('Calculate CCG (Ctrl)');
    CalCCGAlignedRunCtrlOnly(path,fileName,onlyRun,mazeSess);

    %% calculate theta modulation
    disp('Calculate theta modulation (Ctrl)');
    ThetaModulationAlignedRunCtrlOnly(path,fileName,onlyRun,mazeSess);
        
    disp('Extract bursting spikes Hilbert (Ctrl)');
    methodTheta = 0; 
    BurstAllAlignedRunCtrlOnly(path,fileName,methodTheta,0,onlyRun,mazeSess);
    
    disp('Extract bursting spikes Linear (Ctrl)');
    methodTheta = 1;
    BurstAllAlignedRunCtrlOnly(path,fileName,methodTheta,0,onlyRun,mazeSess);
    
    close all;

    % calculate spatial information
    disp('Calculate spatial information (Ctrl)');
    GetFRMapInfo_AlignedCtrlOnly(path,fileName,mazeSess,onlyRun,intervalTSpInfo);

    disp('Calculate spatial information (aligned to cue) (Ctrl)');
    GetFRMapInfo_AlignedCueCtrlOnly(path,fileName,mazeSess,onlyRun,intervalTSpInfo);
    
    disp('Firing field detection based on correlation distance and spatial information (Ctrl)')
    meanNeuronSpikeCorrTCtrlOnly(path,fileName,onlyRun,mazeSess,intervalTSpInfo);
    FieldDetectionAlignedCtrlOnly(path,fileName,onlyRun,0,mazeSess,intervalTSpInfo); % mazeSess is only used for plotting purpose
        
    clearvars;
end
