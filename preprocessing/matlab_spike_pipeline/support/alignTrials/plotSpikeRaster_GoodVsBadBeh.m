function plotSpikeRaster_GoodVsBadBeh(path, fileName, onlyRun, mazeSess, cond)
% plot spike rasters and separate trials based on the animal behavior

    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    fullPath = [path fileName '_Info.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to run file does not exist');
        return;
    end
    load(fullPath,'beh');
        
    GlobalConst;
    ind = beh.mazeSess == mazeSess;
    
    if(isfield(beh,'indStimLap'))
        trialNo = find(beh.indStimLap(ind) == 0);
    else
        trialNo = 1:length(behPar.indTrBadBeh);
    end
    
    indBadBeh = behPar.indTrBadBeh;
    
    goodTrials = intersect(find(indBadBeh == 0),trialNo);
    badTrials = intersect(find(indBadBeh == 1),trialNo);
    trialNo = [goodTrials badTrials];
    disp(['Number of good trials = ' num2str(length(goodTrials))]);
    
    plotSpikeRaster_aligned(path, fileName, onlyRun, mazeSess, trialNo, cond);
    
    %pause; 
    %close all;
end
