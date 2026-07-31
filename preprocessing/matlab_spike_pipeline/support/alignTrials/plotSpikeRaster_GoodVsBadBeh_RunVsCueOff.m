function plotSpikeRaster_GoodVsBadBeh_RunVsCueOff(path, fileName, onlyRun, mazeSess)
% plot spike rasters and separate trials based on the animal behavior

    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    indBadBeh = behPar.indTrBadBeh;
    
    goodTrials = find(indBadBeh == 0);
    badTrials = find(indBadBeh == 1);
    trialNo = [goodTrials badTrials];
    disp(['Number of good trials = ' num2str(length(goodTrials))]);
    
    plotSpikeRaster_alignedRunVsCueOff(path, fileName, onlyRun, mazeSess, trialNo);
    
    %pause; 
    %close all;
end
