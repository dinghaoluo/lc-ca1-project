function meanNeuronSpikeCorrT(path,fileName,onlyRun,mazeSess,intervalT)
% single neuron level mean spike correlation across trials
% meanNeuronSpikeCorrT('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,10)

    if(nargin == 4)
        intervalT = 0;
    end
    fullPath = [path fileName '_spikesCorrTAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    fileNameCorr = [fileName '_meanSpikesCorrTAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _spikesCorrTAligned_Run file does not exist');
        return;
    end
    load(fullPath);
    
    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    neuronNo = length(spikeCorrTRun);
    trialNo = size(spikeCorrTRun{1},1);
    nElem = (trialNo*trialNo-trialNo)/2;
    
    indBadBeh = behPar.indTrBadBeh;
    indGoodTr = indBadBeh == 0;
    indBadTr = indBadBeh == 1;
    nGoodTr = sum(indGoodTr);
    nBadTr = sum(indBadTr);
    nElemGood = (nGoodTr*nGoodTr-nGoodTr)/2;
    nElemBad = (nBadTr*nBadTr-nBadTr)/2;
       
    for n = 1:neuronNo
        corrArr = triu(spikeCorrTRun{n},1);
        corrArr = corrArr(:);        
        meanCorrTRun.mean(n) = sum(corrArr(isnan(corrArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrRun{n} == 1);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanCorrTRun.meanNZ(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanCorrTRun.nNonZeroTr(n) = nNonZeroTr;
        
        corrArr = triu(spikeCorrTRew{n},1);
        corrArr = corrArr(:);  
        meanCorrTRew.mean(n) = sum(corrArr(isnan(corrArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrRew{n} == 1);
        nElemNonZeroRew = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanCorrTRew.meanNZ(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRew;
        meanCorrTRew.nNonZeroTr(n) = nNonZeroTr;
        
        corrArr = triu(spikeCorrTCue{n},1);
        corrArr = corrArr(:);
        meanCorrTCue.mean(n) = sum(corrArr(isnan(corrArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrCue{n} == 1);
        nElemNonZeroCue = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanCorrTCue.meanNZ(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroCue;
        meanCorrTCue.nNonZeroTr(n) = nNonZeroTr;
        
        goodTrCorr = spikeCorrTRun{n}(indGoodTr,indGoodTr);
        goodTrCorr = triu(goodTrCorr,1);
        goodTrCorr = goodTrCorr(:);
        meanCorrTRun.meanGood(n) = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGood;
        indGoodNonZeroTrRun = indBadBeh == 0 & nonZeroTrRun{n} == 1;
        nGoodNonZeroTrRun = sum(indGoodNonZeroTrRun);
        nElemGoodNonZeroTrRun = (nGoodNonZeroTrRun*nGoodNonZeroTrRun-nGoodNonZeroTrRun)/2;
        meanCorrTRun.meanGoodNZ(n) = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGoodNonZeroTrRun;
        meanCorrTRun.nGoodNonZeroTr(n) = nGoodNonZeroTrRun;
        
        goodTrCorr = triu(spikeCorrTRew{n}(indGoodTr,indGoodTr),1);
        goodTrCorr = goodTrCorr(:);
        meanCorrTRew.meanGood(n) = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGood;
        indGoodNonZeroTrRew = indBadBeh == 0 & nonZeroTrRew{n} == 1;
        nGoodNonZeroTrRew = sum(indGoodNonZeroTrRew);
        nElemGoodNonZeroTrRew = (nGoodNonZeroTrRew*nGoodNonZeroTrRew-nGoodNonZeroTrRew)/2;
        meanCorrTRew.meanGoodNZ(n) = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGoodNonZeroTrRew;
        meanCorrTRew.nGoodNonZeroTr(n) = nGoodNonZeroTrRew;
        
        goodTrCorr = triu(spikeCorrTCue{n}(indGoodTr,indGoodTr),1);
        goodTrCorr = goodTrCorr(:);
        meanCorrTCue.meanGood(n) = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGood;
        indGoodNonZeroTrCue = indBadBeh == 0 & nonZeroTrCue{n} == 1;
        nGoodNonZeroTrCue = sum(indGoodNonZeroTrCue);
        nElemGoodNonZeroTrCue = (nGoodNonZeroTrCue*nGoodNonZeroTrCue-nGoodNonZeroTrCue)/2;
        meanCorrTCue.meanGoodNZ(n) = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGoodNonZeroTrCue;
        meanCorrTCue.nGoodNonZeroTr(n) = nGoodNonZeroTrCue;
        
        badTrCorr = triu(spikeCorrTRun{n}(indBadTr,indBadTr),1);
        badTrCorr = badTrCorr(:);
        meanCorrTRun.meanBad(n) = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBad;
        indBadNonZeroTrRun = indBadBeh == 1 & nonZeroTrRun{n} == 1;
        nBadNonZeroTrRun = sum(indBadNonZeroTrRun);
        nElemBadNonZeroTrRun = (nBadNonZeroTrRun*nBadNonZeroTrRun-nBadNonZeroTrRun)/2;
        meanCorrTRun.meanBadNZ(n) = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBadNonZeroTrRun;
        meanCorrTRun.nBadNonZeroTr(n) = nBadNonZeroTrRun;
        
        badTrCorr = triu(spikeCorrTRew{n}(indBadTr,indBadTr),1);
        badTrCorr = badTrCorr(:);
        meanCorrTRew.meanBad(n) = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBad;
        indBadNonZeroTrRew = indBadBeh == 1 & nonZeroTrRew{n} == 1;
        nBadNonZeroTrRew = sum(indBadNonZeroTrRew);
        nElemBadNonZeroTrRew = (nBadNonZeroTrRew*nBadNonZeroTrRew-nBadNonZeroTrRew)/2;
        meanCorrTRew.meanBadNZ(n) = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBadNonZeroTrRew;
        meanCorrTRew.nBadNonZeroTr(n) = nBadNonZeroTrRew;
        
        badTrCorr = triu(spikeCorrTCue{n}(indBadTr,indBadTr),1);
        badTrCorr = badTrCorr(:);
        meanCorrTCue.meanBad(n) = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBad;
        indBadNonZeroTrCue = indBadBeh == 1 & nonZeroTrCue{n} == 1;
        nBadNonZeroTrCue = sum(indBadNonZeroTrCue);
        nElemBadNonZeroTrCue = (nBadNonZeroTrCue*nBadNonZeroTrCue-nBadNonZeroTrCue)/2;
        meanCorrTCue.meanBadNZ(n) = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBadNonZeroTrCue;
        meanCorrTCue.nBadNonZeroTr(n) = nBadNonZeroTrCue;
    end
    
    % compare single neuron correlation between good and bad trials --
    % aligned to run
    plotCompCorr(meanCorrTRun.meanGood,meanCorrTRun.meanBad,...
        'Neu Tr CorrT Run - good trials','Neu Tr CorrT Run - bad trials');
    
%     % compare single neuron correlation between good and bad trials --
%     % aligned to reward
%     plotCompCorr(meanCorrTRew.meanGood,meanCorrTRew.meanBad,...
%         'Neu Tr CorrT Rew - good trials','Neu Tr CorrT Rew - bad trials');
%         
%     % compare single neuron correlation between good and bad trials --
%     % aligned to cue
%     plotCompCorr(meanCorrTCue.meanGood,meanCorrTCue.meanBad,...
%         'Neu Tr CorrT Cue - good trials','Neu Tr CorrT Cue - bad trials');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward
%     plotCompCorr(meanCorrTRun.mean,meanCorrTRew.mean,...
%         'Neu Tr CorrT Run','Neu Tr CorrT Rew');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue
%     plotCompCorr(meanCorrTRun.mean,meanCorrTCue.mean,...
%         'Neu Tr CorrT Run','Neu Tr CorrT Cue');
%     
    % compare single neuron correlation between aligned to run and aligned
    % to reward for good trials
    plotCompCorr(meanCorrTRun.meanGood,meanCorrTRew.meanGood,...
        'Neu Tr CorrT Run - good trials','Neu Tr CorrT Rew - good trials');
   
    % compare single neuron correlation between aligned to run and aligned
    % to cue for good trials
    plotCompCorr(meanCorrTRun.meanGood,meanCorrTCue.meanGood,...
        'Neu Tr CorrT Run - good trials','Neu Tr CorrT Cue - good trials');
%     
%     %% only consider nonzero trials
%     % compare single neuron correlation between good and bad trials --
%     % aligned to run
%     indNeu = meanCorrTRun.nGoodNonZeroTr > 10 & meanCorrTRun.nBadNonZeroTr > 5;
%     plotCompCorr(meanCorrTRun.meanGoodNZ(indNeu),meanCorrTRun.meanBadNZ(indNeu),...
%         'Neu None-zero Tr CorrT Run - good trials','Neu None-zero Tr CorrT Run - bad trials');
%     
%     % compare single neuron correlation between good and bad trials --
%     % aligned to reward
%     indNeu = meanCorrTRew.nGoodNonZeroTr > 10 & meanCorrTRew.nBadNonZeroTr > 5;
%     plotCompCorr(meanCorrTRew.meanGoodNZ(indNeu),meanCorrTRew.meanBadNZ(indNeu),...
%         'Neu None-zero Tr CorrT Rew - good trials','Neu None-zero Tr CorrT Rew - bad trials');
%         
%     % compare single neuron correlation between good and bad trials --
%     % aligned to cue
%     indNeu = meanCorrTCue.nGoodNonZeroTr > 10 & meanCorrTCue.nBadNonZeroTr > 5;
%     plotCompCorr(meanCorrTCue.meanGoodNZ(indNeu),meanCorrTCue.meanBadNZ(indNeu),...
%         'Neu None-zero Tr CorrT Cue - good trials','Neu None-zero Tr CorrT Cue - bad trials');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward
%     indNeu = meanCorrTRun.nNonZeroTr > 10 & meanCorrTRew.nNonZeroTr > 10;
%     plotCompCorr(meanCorrTRun.meanNZ(indNeu),meanCorrTRew.meanNZ(indNeu),...
%         'Neu None-zero Tr CorrT Run ','Neu None-zero Tr CorrT Rew');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue
%     indNeu = meanCorrTRun.nNonZeroTr > 10 & meanCorrTCue.nNonZeroTr > 10;
%     plotCompCorr(meanCorrTRun.meanNZ(indNeu),meanCorrTCue.meanNZ(indNeu),...
%         'Neu None-zero Tr CorrT Run','Neu None-zero Tr CorrT Cue');
%     
    % compare single neuron correlation between aligned to run and aligned
    % to reward for good trials
    indNeu = meanCorrTRun.nGoodNonZeroTr > 5 & meanCorrTRew.nGoodNonZeroTr > 5;
    if(sum(indNeu) > 0)
        plotCompCorr(meanCorrTRun.meanGoodNZ(indNeu),meanCorrTRew.meanGoodNZ(indNeu),...
            'Neu None-zero Tr CorrT Run - good trials','Neu None-zero Tr CorrT Rew - good trials');
    end
   
    % compare single neuron correlation between aligned to run and aligned
    % to cue for good trials
    indNeu = meanCorrTRun.nGoodNonZeroTr > 5 & meanCorrTCue.nGoodNonZeroTr > 5;
    if(sum(indNeu) > 0)
        plotCompCorr(meanCorrTRun.meanGoodNZ(indNeu),meanCorrTCue.meanGoodNZ(indNeu),...
            'Neu None-zero Tr CorrT Run - good trials','Neu None-zero Tr CorrT Cue - good trials');
    end
    
    save([path fileNameCorr],'meanCorrTRun','meanCorrTRew','meanCorrTCue');
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
