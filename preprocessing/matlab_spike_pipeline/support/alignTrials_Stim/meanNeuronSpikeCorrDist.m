function meanNeuronSpikeCorrDist(path,fileName,onlyRun,mazeSess,intervalD)
% single neuron level mean spike correlation across trials
% meanNeuronSpikeCorrDist('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,500)

    if(nargin == 4)
        intervalD = 0;
    end
    fullPath = [path fileName '_spikesCorrDistAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intD' ...
            num2str(intervalD) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _spikesCorrDistAligned_Run file does not exist');
        return;
    end
    load(fullPath);
    
    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    fileNameCorr = [fileName '_meanSpikesCorrDistAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intD' ...
            num2str(intervalD) '.mat'];
    
    neuronNo = length(spikeCorrDistRun);
    trialNo = size(spikeCorrDistRun{1},1);
    nElem = (trialNo*trialNo-trialNo)/2;
    
    indBadBeh = behPar.indTrBadBeh;
    indGoodTr = indBadBeh == 0;
    indBadTr = indBadBeh == 1;
    nGoodTr = sum(indGoodTr);
    nBadTr = sum(indBadTr);
    nElemGood = (nGoodTr*nGoodTr-nGoodTr)/2;
    nElemBad = (nBadTr*nBadTr-nBadTr)/2;
    for n = 1:neuronNo
        corrArr = triu(spikeCorrDistRun{n},1);
        corrArr = corrArr(:); 
        meanCorrDistRun.mean(n) = sum(corrArr(isnan(corrArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrRun{n} == 1);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanCorrDistRun.meanNZ(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanCorrDistRun.nNonZeroTr(n) = nNonZeroTr;
        
        corrArr = triu(spikeCorrDistRew{n},1);
        corrArr = corrArr(:);  
        meanCorrDistRew.mean(n) = sum(corrArr(isnan(corrArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrRew{n} == 1);
        nElemNonZeroRew = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanCorrDistRew.meanNZ(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRew;
        meanCorrDistRew.nNonZeroTr(n) = nNonZeroTr;
        
        corrArr = triu(spikeCorrDistCue{n},1);
        corrArr = corrArr(:);
        meanCorrDistCue.mean(n) = sum(corrArr(isnan(corrArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrCue{n} == 1);
        nElemNonZeroCue = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanCorrDistCue.meanNZ(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroCue;
        meanCorrDistCue.nNonZeroTr(n) = nNonZeroTr;
        
        goodTrCorr = spikeCorrDistRun{n}(indGoodTr,indGoodTr);
        goodTrCorr = triu(goodTrCorr,1);
        goodTrCorr = goodTrCorr(:);
        meanCorrDistRun.meanGood(n) = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGood;
        indGoodNonZeroTrRun = indBadBeh == 0 & nonZeroTrRun{n} == 1;
        nGoodNonZeroTrRun = sum(indGoodNonZeroTrRun);
        nElemGoodNonZeroTrRun = (nGoodNonZeroTrRun*nGoodNonZeroTrRun-nGoodNonZeroTrRun)/2;
        meanCorrDistRun.meanGoodNZ(n) = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGoodNonZeroTrRun;
        meanCorrDistRun.nGoodNonZeroTr(n) = nGoodNonZeroTrRun;
        
        goodTrCorr = triu(spikeCorrDistRew{n}(indGoodTr,indGoodTr),1);
        goodTrCorr = goodTrCorr(:);
        meanCorrDistRew.meanGood(n) = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGood;
        indGoodNonZeroTrRew = indBadBeh == 0 & nonZeroTrRew{n} == 1;
        nGoodNonZeroTrRew = sum(indGoodNonZeroTrRew);
        nElemGoodNonZeroTrRew = (nGoodNonZeroTrRew*nGoodNonZeroTrRew-nGoodNonZeroTrRew)/2;
        meanCorrDistRew.meanGoodNZ(n) = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGoodNonZeroTrRew;
        meanCorrDistRew.nGoodNonZeroTr(n) = nGoodNonZeroTrRew;
        
        goodTrCorr = triu(spikeCorrDistCue{n}(indGoodTr,indGoodTr),1);
        goodTrCorr = goodTrCorr(:);
        meanCorrDistCue.meanGood(n) = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGood;
        indGoodNonZeroTrCue = indBadBeh == 0 & nonZeroTrCue{n} == 1;
        nGoodNonZeroTrCue = sum(indGoodNonZeroTrCue);
        nElemGoodNonZeroTrCue = (nGoodNonZeroTrCue*nGoodNonZeroTrCue-nGoodNonZeroTrCue)/2;
        meanCorrDistCue.meanGoodNZ(n) = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGoodNonZeroTrCue;
        meanCorrDistCue.nGoodNonZeroTr(n) = nGoodNonZeroTrCue;
        
        badTrCorr = triu(spikeCorrDistRun{n}(indBadTr,indBadTr),1);
        badTrCorr = badTrCorr(:);
        meanCorrDistRun.meanBad(n) = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBad;
        indBadNonZeroTrRun = indBadBeh == 1 & nonZeroTrRun{n} == 1;
        nBadNonZeroTrRun = sum(indBadNonZeroTrRun);
        nElemBadNonZeroTrRun = (nBadNonZeroTrRun*nBadNonZeroTrRun-nBadNonZeroTrRun)/2;
        meanCorrDistRun.meanBadNZ(n) = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBadNonZeroTrRun;
        meanCorrDistRun.nBadNonZeroTr(n) = nBadNonZeroTrRun;
        
        badTrCorr = triu(spikeCorrDistRew{n}(indBadTr,indBadTr),1);
        badTrCorr = badTrCorr(:);
        meanCorrDistRew.meanBad(n) = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBad;
        indBadNonZeroTrRew = indBadBeh == 1 & nonZeroTrRew{n} == 1;
        nBadNonZeroTrRew = sum(indBadNonZeroTrRew);
        nElemBadNonZeroTrRew = (nBadNonZeroTrRew*nBadNonZeroTrRew-nBadNonZeroTrRew)/2;
        meanCorrDistRew.meanBadNZ(n) = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBadNonZeroTrRew;
        meanCorrDistRew.nBadNonZeroTr(n) = nBadNonZeroTrRew;
        
        badTrCorr = spikeCorrDistCue{n}(indBadTr,indBadTr);
        badTrCorr = badTrCorr(:);
        meanCorrDistCue.meanBad(n) = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBad;
        indBadNonZeroTrCue = indBadBeh == 1 & nonZeroTrCue{n} == 1;
        nBadNonZeroTrCue = sum(indBadNonZeroTrCue);
        nElemBadNonZeroTrCue = (nBadNonZeroTrCue*nBadNonZeroTrCue-nBadNonZeroTrCue)/2;
        meanCorrDistCue.meanBadNZ(n) = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBadNonZeroTrCue;
        meanCorrDistCue.nBadNonZeroTr(n) = nBadNonZeroTrCue;
    end
    
    % compare single neuron correlation between good and bad trials --
    % aligned to run
    plotCompCorr(meanCorrDistRun.meanGood,meanCorrDistRun.meanBad,...
        'Neu Tr CorrDist Run - good trials','Neu Tr CorrDist Run - bad trials');
    
%     % compare single neuron correlation between good and bad trials --
%     % aligned to reward
%     plotCompCorr(meanCorrDistRew.meanGood,meanCorrDistRew.meanBad,...
%         'Neu Tr CorrDist Rew - good trials','Neu Tr CorrDist Rew - bad trials');
%         
%     % compare single neuron correlation between good and bad trials --
%     % aligned to cue
%     plotCompCorr(meanCorrDistCue.meanGood,meanCorrDistCue.meanBad,...
%         'Neu Tr CorrDist Cue - good trials','Neu Tr CorrDist Cue - bad trials');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward
%     plotCompCorr(meanCorrDistRun.mean,meanCorrDistRew.mean,...
%         'Neu Tr CorrDist Run','Neu Tr CorrDist Rew');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue
%     plotCompCorr(meanCorrDistRun.mean,meanCorrDistCue.mean,...
%         'Neu Tr CorrDist Run','Neu Tr CorrDist Cue');
%     
    % compare single neuron correlation between aligned to run and aligned
    % to reward for good trials
    plotCompCorr(meanCorrDistRun.meanGood,meanCorrDistRew.meanGood,...
        'Neu Tr CorrDist Run - good trials','Neu Tr CorrDist Rew - good trials');
   
    % compare single neuron correlation between aligned to run and aligned
    % to cue for good trials
    plotCompCorr(meanCorrDistRun.meanGood,meanCorrDistCue.meanGood,...
        'Neu Tr CorrDist Run - good trials','Neu Tr CorrDist Cue - good trials');
    
%     %% only consider nonzero trials
%     % compare single neuron correlation between good and bad trials --
%     % aligned to run
%     indNeu = meanCorrDistRun.nGoodNonZeroTr > 10 & meanCorrDistRun.nBadNonZeroTr > 5;
%     plotCompCorr(meanCorrDistRun.meanGoodNZ(indNeu),meanCorrDistRun.meanBadNZ(indNeu),...
%         'Neu None-zero Tr CorrDist Run - good trials','Neu None-zero Tr CorrDist Run - bad trials');
%     
%     % compare single neuron correlation between good and bad trials --
%     % aligned to reward
%     indNeu = meanCorrDistRew.nGoodNonZeroTr > 10 & meanCorrDistRew.nBadNonZeroTr > 5;
%     plotCompCorr(meanCorrDistRew.meanGoodNZ(indNeu),meanCorrDistRew.meanBadNZ(indNeu),...
%         'Neu None-zero Tr CorrDist Rew - good trials','Neu None-zero Tr CorrDist Rew - bad trials');
%         
%     % compare single neuron correlation between good and bad trials --
%     % aligned to cue
%     indNeu = meanCorrDistCue.nGoodNonZeroTr > 10 & meanCorrDistCue.nBadNonZeroTr > 5;
%     plotCompCorr(meanCorrDistCue.meanGoodNZ(indNeu),meanCorrDistCue.meanBadNZ(indNeu),...
%         'Neu None-zero Tr CorrDist Cue - good trials','Neu None-zero Tr CorrDist Cue - bad trials');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward
%     indNeu = meanCorrDistRun.nNonZeroTr > 10 & meanCorrDistRew.nNonZeroTr > 10;
%     plotCompCorr(meanCorrDistRun.meanNZ(indNeu),meanCorrDistRew.meanNZ(indNeu),...
%         'Neu None-zero Tr CorrDist Run ','Neu None-zero Tr CorrDist Rew');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue
%     indNeu = meanCorrDistRun.nNonZeroTr > 10 & meanCorrDistCue.nNonZeroTr > 10;
%     plotCompCorr(meanCorrDistRun.meanNZ(indNeu),meanCorrDistCue.meanNZ(indNeu),...
%         'Neu None-zero Tr CorrDist Run','Neu None-zero Tr CorrDist Cue');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward for good trials
%     indNeu = meanCorrDistRun.nGoodNonZeroTr > 10 & meanCorrDistRew.nGoodNonZeroTr > 10;
%     plotCompCorr(meanCorrDistRun.meanGoodNZ(indNeu),meanCorrDistRew.meanGoodNZ(indNeu),...
%         'Neu None-zero Tr CorrDist Run - good trials','Neu None-zero Tr CorrDist Rew - good trials');
%    
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue for good trials
%     indNeu = meanCorrDistRun.nGoodNonZeroTr > 10 & meanCorrDistCue.nGoodNonZeroTr > 10;
%     plotCompCorr(meanCorrDistRun.meanGoodNZ(indNeu),meanCorrDistCue.meanGoodNZ(indNeu),...
%         'Neu None-zero Tr CorrDist Run - good trials','Neu None-zero Tr CorrDist Cue - good trials');    
    
    save([path fileNameCorr],'meanCorrDistRun','meanCorrDistRew','meanCorrDistCue');
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
