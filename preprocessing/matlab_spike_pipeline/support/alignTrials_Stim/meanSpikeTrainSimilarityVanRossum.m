function meanSpikeTrainSimilarityVanRossum(path,fileName,onlyRun,mazeSess,tc,intervalT)
% single neuron level mean spike correlation across trials
% van Rossum distance
% e.g.: meanSpikeTrainSimilarityVanRossum('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,100,8)

    fullPath = [path fileName '_spikeTrainSimiVanR_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_tc' ...
            num2str(tc) '_intT' num2str(intervalT) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _spikeTrainSimiVanR_Run file does not exist');
        return;
    end
    load(fullPath);
    
    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    fileNameCorr = [fileName '_mean_spikeTrainSimiVanR_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_tc' ...
            num2str(tc) '_intT' num2str(intervalT) '.mat'];
    
    neuronNo = length(spikeTrainSimVRRun);
    trialNo = size(spikeTrainSimVRRun{1},1);
    nElem = trialNo*trialNo-trialNo;
    
    indBadBeh = behPar.indTrBadBeh;
    indGoodTr = indBadBeh == 0;
    indBadTr = indBadBeh == 1;
    nGoodTr = sum(indGoodTr);
    nBadTr = sum(indBadTr);
    nElemGood = nGoodTr*nGoodTr-nGoodTr;
    nElemBad = nBadTr*nBadTr-nBadTr;
    
    for n = 1:neuronNo
        %% non zero elements
        simArr = triu(spikeTrainSimVRRun{n},1);
        simArr = abs(simArr(:));
        meanSimVRRun.mean(n) = sum(simArr(isnan(simArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrRun{n} == 1);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVRRun.meanNZ(n) = sum(simArr(isnan(simArr) == 0))/nElemNonZeroRun;
        meanSimVRRun.nNonZeroTr(n) = nNonZeroTr;
                
        simArr = triu(spikeTrainSimVRRew{n},1);
        simArr = abs(simArr(:));
        meanSimVRRew.mean(n) = sum(simArr(isnan(simArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrRew{n} == 1);
        nElemNonZeroRew = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVRRew.meanNZ(n) = sum(simArr(isnan(simArr) == 0))/nElemNonZeroRew;
        meanSimVRRew.nNonZeroTr(n) = nNonZeroTr;
        
        simArr = triu(spikeTrainSimVRCue{n},1);
        simArr = abs(simArr(:));
        meanSimVRCue.mean(n) = sum(simArr(isnan(simArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrCue{n} == 1);
        nElemNonZeroCue = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVRCue.meanNZ(n) = sum(simArr(isnan(simArr) == 0))/nElemNonZeroCue;
        meanSimVRCue.nNonZeroTr(n) = nNonZeroTr;
        
        %% non zero elements in good trials
        goodTrSim = spikeTrainSimVRRun{n}(indGoodTr,indGoodTr);
        goodTrSim = triu(goodTrSim,1);
        goodTrSim = abs(goodTrSim(:));
        meanSimVRRun.meanGood(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGood;
        indGoodNonZeroTrRun = indBadBeh == 0 & nonZeroTrRun{n} == 1;
        nGoodNonZeroTrRun = sum(indGoodNonZeroTrRun);
        nElemGoodNonZeroTrRun = (nGoodNonZeroTrRun*nGoodNonZeroTrRun-nGoodNonZeroTrRun)/2;
        meanSimVRRun.meanGoodNZ(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGoodNonZeroTrRun;
        meanSimVRRun.nGoodNonZeroTr(n) = nGoodNonZeroTrRun;
        
        goodTrSim = spikeTrainSimVRRew{n}(indGoodTr,indGoodTr);
        goodTrSim = triu(goodTrSim,1);
        goodTrSim = abs(goodTrSim(:));
        meanSimVRRew.meanGood(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGood;
        indGoodNonZeroTrRew = indBadBeh == 0 & nonZeroTrRew{n} == 1;
        nGoodNonZeroTrRew = sum(indGoodNonZeroTrRew);
        nElemGoodNonZeroTrRew = (nGoodNonZeroTrRew*nGoodNonZeroTrRew-nGoodNonZeroTrRew)/2;
        meanSimVRRew.meanGoodNZ(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGoodNonZeroTrRew;
        meanSimVRRew.nGoodNonZeroTr(n) = nGoodNonZeroTrRew;
        
        goodTrSim = spikeTrainSimVRCue{n}(indGoodTr,indGoodTr);
        goodTrSim = triu(goodTrSim,1);
        goodTrSim = abs(goodTrSim(:));
        meanSimVRCue.meanGood(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGood;
        indGoodNonZeroTrCue = indBadBeh == 0 & nonZeroTrCue{n} == 1;
        nGoodNonZeroTrCue = sum(indGoodNonZeroTrCue);
        nElemGoodNonZeroTrCue = (nGoodNonZeroTrCue*nGoodNonZeroTrCue-nGoodNonZeroTrCue)/2;
        meanSimVRCue.meanGoodNZ(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGoodNonZeroTrCue;
        meanSimVRCue.nGoodNonZeroTr(n) = nGoodNonZeroTrCue;
        
        %% non zero elements in bad trials
        badTrSim = spikeTrainSimVRRun{n}(indBadTr,indBadTr);
        badTrSim = triu(badTrSim,1);
        badTrSim = abs(badTrSim(:));
        meanSimVRRun.meanBad(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBad;
        indBadNonZeroTrRun = indBadBeh == 1 & nonZeroTrRun{n} == 1;
        nBadNonZeroTrRun = sum(indBadNonZeroTrRun);
        nElemBadNonZeroTrRun = (nBadNonZeroTrRun*nBadNonZeroTrRun-nBadNonZeroTrRun)/2;
        meanSimVRRun.meanBadNZ(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBadNonZeroTrRun;
        meanSimVRRun.nBadNonZeroTr(n) = nBadNonZeroTrRun;
        
        badTrSim = spikeTrainSimVRRew{n}(indBadTr,indBadTr);
        badTrSim = triu(badTrSim,1);
        badTrSim = abs(badTrSim(:));
        meanSimVRRew.meanBad(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBad;
        indBadNonZeroTrRew = indBadBeh == 1 & nonZeroTrRew{n} == 1;
        nBadNonZeroTrRew = sum(indBadNonZeroTrRew);
        nElemBadNonZeroTrRew = (nBadNonZeroTrRew*nBadNonZeroTrRew-nBadNonZeroTrRew)/2;
        meanSimVRRew.meanBadNZ(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBadNonZeroTrRew;
        meanSimVRRew.nBadNonZeroTr(n) = nBadNonZeroTrRew;
        
        badTrSim = spikeTrainSimVRCue{n}(indBadTr,indBadTr);
        badTrSim = triu(badTrSim,1);
        badTrSim = abs(badTrSim(:));
        meanSimVRCue.meanBad(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBad;
        indBadNonZeroTrCue = indBadBeh == 1 & nonZeroTrCue{n} == 1;
        nBadNonZeroTrCue = sum(indBadNonZeroTrCue);
        nElemBadNonZeroTrCue = (nBadNonZeroTrCue*nBadNonZeroTrCue-nBadNonZeroTrCue)/2;
        meanSimVRCue.meanBadNZ(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBadNonZeroTrCue;
        meanSimVRCue.nBadNonZeroTr(n) = nBadNonZeroTrCue;
    end
    
    % compare single neuron correlation between good and bad trials --
    % aligned to run
    plotCompCorr(meanSimVRRun.meanGood,meanSimVRRun.meanBad,...
        'Neu Tr SimVR Run - good trials','Neu Tr SimVR Run - bad trials');
    
%     % compare single neuron correlation between good and bad trials --
%     % aligned to reward
%     plotCompCorr(meanSimVRRew.meanGood,meanSimVRRew.meanBad,...
%         'Neu Tr SimVR Rew - good trials','Neu Tr SimVR Rew - bad trials');
%         
%     % compare single neuron correlation between good and bad trials --
%     % aligned to cue
%     plotCompCorr(meanSimVRCue.meanGood,meanSimVRCue.meanBad,...
%         'Neu Tr SimVR Cue - good trials','Neu Tr SimVR Cue - bad trials');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward
%     plotCompCorr(meanSimVRRun.mean,meanSimVRRew.mean,...
%         'Neu Tr SimVR Run','Neu Tr SimVR Rew');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue
%     plotCompCorr(meanSimVRRun.mean,meanSimVRCue.mean,...
%         'Neu Tr SimVR Run','Neu Tr SimVR Cue');
%     
    % compare single neuron correlation between aligned to run and aligned
    % to reward for good trials
    plotCompCorr(meanSimVRRun.meanGood,meanSimVRRew.meanGood,...
        'Neu Tr SimVR Run - good trials','Neu Tr SimVR Rew - good trials');
   
    % compare single neuron correlation between aligned to run and aligned
    % to cue for good trials
    plotCompCorr(meanSimVRRun.meanGood,meanSimVRCue.meanGood,...
        'Neu Tr SimVR Run - good trials','Neu Tr SimVR Cue - good trials');
    
    %% only consider nonzero trials
%     % compare single neuron correlation between good and bad trials --
%     % aligned to run
%     indNeu = meanSimVRRun.nGoodNonZeroTr > 10 & meanSimVRRun.nBadNonZeroTr > 10;
%     plotCompCorr(meanSimVRRun.meanGoodNZ(indNeu),meanSimVRRun.meanBadNZ(indNeu),...
%         'Neu None-zero Tr SimVR Run - good trials','Neu None-zero Tr SimVR Run - bad trials');
%     
%     % compare single neuron correlation between good and bad trials --
%     % aligned to reward
%     indNeu = meanSimVRRew.nGoodNonZeroTr > 10 & meanSimVRRew.nBadNonZeroTr > 10;
%     plotCompCorr(meanSimVRRew.meanGoodNZ(indNeu),meanSimVRRew.meanBadNZ(indNeu),...
%         'Neu None-zero Tr SimVR Rew - good trials','Neu None-zero Tr SimVR Rew - bad trials');
%         
%     % compare single neuron correlation between good and bad trials --
%     % aligned to cue
%     indNeu = meanSimVRCue.nGoodNonZeroTr > 10 & meanSimVRCue.nBadNonZeroTr > 10;
%     plotCompCorr(meanSimVRCue.meanGoodNZ(indNeu),meanSimVRCue.meanBadNZ(indNeu),...
%         'Neu None-zero Tr SimVR Cue - good trials','Neu None-zero Tr SimVR Cue - bad trials');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward
%     indNeu = meanSimVRRun.nNonZeroTr > 10 & meanSimVRRew.nNonZeroTr > 10;
%     plotCompCorr(meanSimVRRun.meanNZ(indNeu),meanSimVRRew.meanNZ(indNeu),...
%         'Neu None-zero Tr SimVR Run ','Neu None-zero Tr SimVR Rew');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue
%     indNeu = meanSimVRRun.nNonZeroTr > 10 & meanSimVRCue.nNonZeroTr > 10;
%     plotCompCorr(meanSimVRRun.meanNZ(indNeu),meanSimVRCue.meanNZ(indNeu),...
%         'Neu None-zero Tr SimVR Run','Neu None-zero Tr SimVR Cue');
%     
    % compare single neuron correlation between aligned to run and aligned
    % to reward for good trials
    indNeu = meanSimVRRun.nGoodNonZeroTr > 5 & meanSimVRRew.nGoodNonZeroTr > 5;
    plotCompCorr(meanSimVRRun.meanGoodNZ(indNeu),meanSimVRRew.meanGoodNZ(indNeu),...
        'Neu None-zero Tr SimVR Run - good trials','Neu None-zero Tr SimVR Rew - good trials');
   
    % compare single neuron correlation between aligned to run and aligned
    % to cue for good trials
    indNeu = meanSimVRRun.nGoodNonZeroTr > 5 & meanSimVRCue.nGoodNonZeroTr > 5;
    plotCompCorr(meanSimVRRun.meanGoodNZ(indNeu),meanSimVRCue.meanGoodNZ(indNeu),...
        'Neu None-zero Tr SimVR Run - good trials','Neu None-zero Tr SimVR Cue - good trials');
    
    save([path fileNameCorr],'meanSimVRRun','meanSimVRRew','meanSimVRCue');
end

function plotCompCorr(x,y,xlab,ylab)
    figure
    plot(x,y,'ro');
    hold on;
    corrAll = [x y];
    maxCorr = max(corrAll);
    minCorr = min(corrAll);
    plot([minCorr maxCorr],[minCorr maxCorr],'k:');
    xlim(gca,[minCorr maxCorr]);
    ylim(gca,[minCorr maxCorr]);
    xlabel(xlab);
    ylabel(ylab);
end
