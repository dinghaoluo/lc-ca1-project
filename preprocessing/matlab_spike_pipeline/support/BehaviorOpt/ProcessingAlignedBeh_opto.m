function ProcessingAlignedBeh_opto(path,fileName,onlyRun,mazeSess)

    GlobalConst;
    
    %% generate behavior data structure
    GenerateBehavDataStructures_smallTMPFI_opto(path,fileName);
    
    for mz = mazeSess
        %% align the spikes based on different run onset/reward/cue onset
        disp('Align trials based on run onset')
        alignToRunOnsetBeh1_opto(path, fileName,mz);

        disp('Align trials based on reward onset')
        alignToRewardBeh(path, fileName,mz);

        disp('Align trials based on cue onset')
        alignToCueBeh(path, fileName,mz);

        %% get behavior parameters
        disp('Get behavior parameters')
        getBehParametersBeh_opto(path,fileName,mz);

        %% calculate lick over distance
        disp('Calculate lick over distance')
        LickOverDistBeh_opto(path, fileName, mz)

        %% calculate running speed over distance
        disp('Calculate running speed over distance')
        RunSpeedOverDistBeh_opto(path, fileName, onlyRun, mz);
        
        %% calculate running speed over time
        disp('Calculate running speed over time')
        RunSpeedOverTimeBeh_opto(path, fileName, onlyRun, mz);
        
        %% calculate trial by trial licking similarity
        LickOverDistSim_opto(path,fileName,mz);
        
        %% calculate trial by trial speed similarity
        RunSpeedOverDistSim_opto(path,fileName,onlyRun,mz);
    end
end
