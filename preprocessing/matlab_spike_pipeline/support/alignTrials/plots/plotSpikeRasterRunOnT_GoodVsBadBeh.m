function plotSpikeRasterRunOnT_GoodVsBadBeh(path, fileName, onlyRun, mazeSess, NeuronNo)
% plot spike rasters and separate trials based on the animal behavior

    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    indBadBeh = behPar.indTrBadBeh;
    
    goodTrials = find(indBadBeh == 0);
    goodTrials = goodTrials(goodTrials ~= 1);
    badTrials = find(indBadBeh == 1);
    badTrials = badTrials(badTrials ~= 1);
    trialNo = [goodTrials badTrials];
    disp(['Number of good trials = ' num2str(length(goodTrials))]);
    nGoodTrials = length(goodTrials);
    startCueToRunBad = behPar.startCueToRun(badTrials);
    rewardToRunBad = behPar.rewardToRun(badTrials);
    
    %% plot spike raster in the order of trials (for bad trials)
%     plotSpikeRasterAndRunOnT_aligned_notOrdered(path,fileName,onlyRun,...
%         mazeSess,badTrials,NeuronNo,length(badTrials),'badTr');
    
    plotSpikeRasterAndRunOnT_aligned(path,fileName,onlyRun,...
        mazeSess,badTrials,NeuronNo,length(badTrials),startCueToRunBad,rewardToRunBad,...
        'badTr');
    
    %% plot spike raster in the order of trials (for good trials)
    plotSpikeRasterAndRunOnT_aligned_notOrdered(path,fileName,onlyRun,...
        mazeSess,goodTrials,NeuronNo,nGoodTrials,'goodTr');
    
    %% plot spike raster in the order of trials
    plotSpikeRasterAndRunOnT_aligned_notOrdered(path,fileName,onlyRun,...
        mazeSess,trialNo,NeuronNo,nGoodTrials,'');
    
    %% plot spike raster in the order of startCue to Runonset
    startCueToRunGood = behPar.startCueToRun(goodTrials);
    rewardToRunGood = behPar.rewardToRun(goodTrials);
    [startCueToRunGood,indGood] = sort(startCueToRunGood,'descend');
    goodTrials = goodTrials(indGood);
       
    startCueToRunBad = behPar.startCueToRun(badTrials);
    rewardToRunBad = behPar.rewardToRun(badTrials);
    [startCueToRunBad,indBad] = sort(startCueToRunBad,'descend');
    badTrials = badTrials(indBad);
%     trialNo = [goodTrials badTrials];  
%     startCueToRun = [startCueToRunGood startCueToRunBad];
%     rewardToRun = [rewardToRunGood(indGood) rewardToRunBad(indBad)];
    trialNo = [goodTrials];
    startCueToRun = [startCueToRunGood];
    rewardToRun = [rewardToRunGood(indGood)];
    nGoodTrials = length(indGood);
    
    plotSpikeRasterAndRunOnT_aligned(path,fileName,onlyRun,mazeSess,trialNo,NeuronNo,...
        nGoodTrials,startCueToRun,rewardToRun,'');
    
end
