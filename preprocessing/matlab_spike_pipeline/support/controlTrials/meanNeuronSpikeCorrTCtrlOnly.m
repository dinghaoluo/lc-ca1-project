function meanNeuronSpikeCorrTCtrlOnly(path,fileName,onlyRun,mazeSess,intervalT)
% single neuron level mean spike correlation across trials
% meanNeuronSpikeCorrT('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,10)

    if(nargin == 4)
        intervalT = 0;
    end
    fullPath = [path fileName '_spikesCorrTAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    fileNameCorr = [fileName '_meanSpikesCorrTAlignedCtrl_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
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
    
    fullPath = [path fileName '_PeakFRAlignedCtrl_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];    
    if(exist(fullPath) == 0)
        disp('The _PeakFRAlignedCtrl file does not exist');
        return;
    end
    load(fullPath,'trialNoNonStim');
    
    neuronNo = length(spikeCorrTRun);
    trialNo = length(trialNoNonStim);
    nElem = (trialNo*trialNo-trialNo)/2;
    
    for n = 1:neuronNo
        corrArr = triu(spikeCorrTRun{n}(trialNoNonStim,trialNoNonStim),1);
        corrArr = corrArr(:);        
        meanCorrTRunCtrl.mean(n) = sum(corrArr(isnan(corrArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrRun{n}(trialNoNonStim) == 1);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanCorrTRunCtrl.meanNZ(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanCorrTRunCtrl.nNonZeroTr(n) = nNonZeroTr;
        
        corrArr = triu(spikeCorrTRew{n}(trialNoNonStim,trialNoNonStim),1);
        corrArr = corrArr(:);  
        meanCorrTRewCtrl.mean(n) = sum(corrArr(isnan(corrArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrRew{n}(trialNoNonStim) == 1);
        nElemNonZeroRew = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanCorrTRewCtrl.meanNZ(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRew;
        meanCorrTRewCtrl.nNonZeroTr(n) = nNonZeroTr;
        
        corrArr = triu(spikeCorrTCue{n}(trialNoNonStim,trialNoNonStim),1);
        corrArr = corrArr(:);
        meanCorrTCueCtrl.mean(n) = sum(corrArr(isnan(corrArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrCue{n}(trialNoNonStim) == 1);
        nElemNonZeroCue = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanCorrTCueCtrl.meanNZ(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroCue;
        meanCorrTCueCtrl.nNonZeroTr(n) = nNonZeroTr;
        
    end
    
    % compare single neuron correlation between aligned to run and aligned
    % to reward for good trials
    plotCompCorr(meanCorrTRunCtrl.mean,meanCorrTRewCtrl.mean,...
        'Neu Tr CorrT Run - ctrl trials','Neu Tr CorrT Rew - ctrl trials');
   
    % compare single neuron correlation between aligned to run and aligned
    % to cue for good trials
    plotCompCorr(meanCorrTRunCtrl.mean,meanCorrTCueCtrl.mean,...
        'Neu Tr CorrT Run - ctrl trials','Neu Tr CorrT Cue - ctrl trials');

    % compare single neuron correlation between aligned to run and aligned
    % to reward for good trials
    indNeu = meanCorrTRunCtrl.nNonZeroTr > 5 & meanCorrTRewCtrl.nNonZeroTr > 5;
    if(sum(indNeu) > 0)
        plotCompCorr(meanCorrTRunCtrl.meanNZ(indNeu),meanCorrTRewCtrl.meanNZ(indNeu),...
            'Neu None-zero Tr CorrT Run - ctrl trials','Neu None-zero Tr CorrT Rew - ctrl trials');
    end
   
    % compare single neuron correlation between aligned to run and aligned
    % to cue for good trials
    indNeu = meanCorrTRunCtrl.nNonZeroTr > 5 & meanCorrTCueCtrl.nNonZeroTr > 5;
    if(sum(indNeu) > 0)
        plotCompCorr(meanCorrTRunCtrl.meanNZ(indNeu),meanCorrTCueCtrl.meanNZ(indNeu),...
            'Neu None-zero Tr CorrT Run - ctrl trials','Neu None-zero Tr CorrT Cue - ctrl trials');
    end
    
    save([path fileNameCorr],'meanCorrTRunCtrl','meanCorrTRewCtrl','meanCorrTCueCtrl');
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
