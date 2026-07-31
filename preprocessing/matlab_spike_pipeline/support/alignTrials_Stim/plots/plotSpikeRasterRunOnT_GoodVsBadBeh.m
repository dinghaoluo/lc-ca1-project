function plotSpikeRasterRunOnT_GoodVsBadBeh(path, fileName, onlyRun, NeuronNo)
% plot spike rasters and separate trials based on the animal behavior

    fullPath = [path fileName '_behPar.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    indBadBeh = behPar.indTrBadBeh;
    
    goodTrials = find(indBadBeh == 0);
    goodTrials = goodTrials(goodTrials ~= 1);
    badTrials = find(indBadBeh == 1);
    badTrials = goodTrials(badTrials ~= 1);
    trialNo = [goodTrials badTrials];
    disp(['Number of good trials = ' num2str(length(goodTrials))]);
    
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
    plotSpikeRasterAndRunOnT_aligned(path,fileName,onlyRun,trialNo,NeuronNo,...
        nGoodTrials,startCueToRun,rewardToRun);
end
