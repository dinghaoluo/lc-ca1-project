function plotSpikeRaster_GoodVsBadBehNoStim(path, fileName, onlyRun, mazeSess)
% plot spike rasters and separate trials based on the animal behavior

    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    indnoStim = find(behPar.stimOn == 1);
    if(isempty(indnoStim))
        indBadBeh = behPar.indTrBadBeh;
    else
        indnoStim = 1:indnoStim(1)-1;
        indBadBeh = behPar.indTrBadBeh(indnoStim);
    end
    
    goodTrials = find(indBadBeh == 0);
    badTrials = find(indBadBeh == 1);
    trialNo = [goodTrials badTrials];
    disp(['Number of good trials = ' num2str(length(goodTrials))]);
    
    plotSpikeRaster_aligned(path, fileName, onlyRun, mazeSess, trialNo);
    
    pause; 
    close all;
end
