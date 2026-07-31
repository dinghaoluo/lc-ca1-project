function ProcessingAlignedBeh(path,fileName,onlyRun,mazeSess)

    GlobalConst;
    
    %% generate behavior data structure
    GenerateBehavDataStructures_smallTMPFI(path,fileName);
    
    %% align the spikes based on different run onset/reward/cue onset
    disp('Align trials based on run onset')
    alignToRunOnsetBeh(path, fileName,mazeSess);
    
    disp('Align trials based on reward onset')
    alignToRewardBeh(path, fileName,mazeSess);
    
    disp('Align trials based on cue onset')
    alignToCueBeh(path, fileName,mazeSess);
    
    %% get behavior parameters
    disp('Get behavior parameters')
    getBehParametersBeh(path,fileName,mazeSess);
        
    %% calculate lick over distance
    disp('Calculate lick over distance')
    LickOverDistBeh(path, fileName, mazeSess)
    
    %% calculate running speed over distance
    disp('Calculate running speed over distance')
    RunSpeedOverDistBeh(path, fileName, onlyRun, mazeSess);
end
