function ProcessingMice_smTrSpeedLick(path, fileName, onlyRun, mazeSess)
% calculate licking and speed over distance    

    disp('lick over distance')
    LickOverDist2P(path, fileName, mazeSess);
    
    disp('speed over distance')
    RunSpeedOverDist2P(path, fileName, onlyRun, mazeSess);
    
    clearvars;
end
