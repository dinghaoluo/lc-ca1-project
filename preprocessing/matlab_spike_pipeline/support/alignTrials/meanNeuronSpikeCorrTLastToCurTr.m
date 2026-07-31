function meanNeuronSpikeCorrTLastToCurTr(path,fileName,onlyRun,mazeSess,intervalT,intervalTMin)
% single neuron level mean spike correlation across trials
% meanNeuronSpikeCorrT('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,10)

    if(nargin == 4)
        intervalT = 0;
    end
    fullPath = [path fileName '_spikesCorrTAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    fileNameCorr = [fileName '_meanSpikesCorrTAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalTMin) '-' num2str(intervalT) '.mat'];
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
    
    neuronNo = length(spikeCorrTRun_LasttoCurTr);
    trialNo = size(spikeCorrTRun_LasttoCurTr{1},1);
    nElem = (trialNo*trialNo-trialNo)/2;
    
    indBadBeh = behPar.indTrBadBeh;
    indGoodTr = indBadBeh == 0;
    indBadTr = indBadBeh == 1;
    nGoodTr = sum(indGoodTr);
    nBadTr = sum(indBadTr);
    nElemGood = (nGoodTr*nGoodTr-nGoodTr)/2;
    nElemBad = (nBadTr*nBadTr-nBadTr)/2;
       
    for n = 1:neuronNo
        corrArr = triu(spikeCorrTRun_LasttoCurTr{n},1);
        corrArr = corrArr(:);        
        meanCorrTRun_LasttoCurTr.mean(n) = sum(corrArr(isnan(corrArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrRun_LasttoCurTr{n} == 1);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanCorrTRun_LasttoCurTr.meanNZ(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanCorrTRun_LasttoCurTr.nNonZeroTr(n) = nNonZeroTr;
        
        corrArr = triu(spikeCorrTRew_LasttoCurTr{n},1);
        corrArr = corrArr(:);  
        meanCorrTRew_LasttoCurTr.mean(n) = sum(corrArr(isnan(corrArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrRew_LasttoCurTr{n} == 1);
        nElemNonZeroRew = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanCorrTRew_LasttoCurTr.meanNZ(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRew;
        meanCorrTRew_LasttoCurTr.nNonZeroTr(n) = nNonZeroTr;
        
        corrArr = triu(spikeCorrTCue_LasttoCurTr{n},1);
        corrArr = corrArr(:);
        meanCorrTCue_LasttoCurTr.mean(n) = sum(corrArr(isnan(corrArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrCue_LasttoCurTr{n} == 1);
        nElemNonZeroCue = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanCorrTCue_LasttoCurTr.meanNZ(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroCue;
        meanCorrTCue_LasttoCurTr.nNonZeroTr(n) = nNonZeroTr;
        
        goodTrCorr = spikeCorrTRun_LasttoCurTr{n}(indGoodTr,indGoodTr);
        goodTrCorr = triu(goodTrCorr,1);
        goodTrCorr = goodTrCorr(:);
        meanCorrTRun_LasttoCurTr.meanGood(n) = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGood;
        indGoodNonZeroTrRun = indBadBeh == 0 & nonZeroTrRun_LasttoCurTr{n} == 1;
        nGoodNonZeroTrRun = sum(indGoodNonZeroTrRun);
        nElemGoodNonZeroTrRun = (nGoodNonZeroTrRun*nGoodNonZeroTrRun-nGoodNonZeroTrRun)/2;
        meanCorrTRun_LasttoCurTr.meanGoodNZ(n) = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGoodNonZeroTrRun;
        meanCorrTRun_LasttoCurTr.nGoodNonZeroTr(n) = nGoodNonZeroTrRun;
        
        goodTrCorr = triu(spikeCorrTRew_LasttoCurTr{n}(indGoodTr,indGoodTr),1);
        goodTrCorr = goodTrCorr(:);
        meanCorrTRew_LasttoCurTr.meanGood(n) = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGood;
        indGoodNonZeroTrRew = indBadBeh == 0 & nonZeroTrRew_LasttoCurTr{n} == 1;
        nGoodNonZeroTrRew = sum(indGoodNonZeroTrRew);
        nElemGoodNonZeroTrRew = (nGoodNonZeroTrRew*nGoodNonZeroTrRew-nGoodNonZeroTrRew)/2;
        meanCorrTRew_LasttoCurTr.meanGoodNZ(n) = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGoodNonZeroTrRew;
        meanCorrTRew_LasttoCurTr.nGoodNonZeroTr(n) = nGoodNonZeroTrRew;
        
        goodTrCorr = triu(spikeCorrTCue_LasttoCurTr{n}(indGoodTr,indGoodTr),1);
        goodTrCorr = goodTrCorr(:);
        meanCorrTCue_LasttoCurTr.meanGood(n) = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGood;
        indGoodNonZeroTrCue = indBadBeh == 0 & nonZeroTrCue_LasttoCurTr{n} == 1;
        nGoodNonZeroTrCue = sum(indGoodNonZeroTrCue);
        nElemGoodNonZeroTrCue = (nGoodNonZeroTrCue*nGoodNonZeroTrCue-nGoodNonZeroTrCue)/2;
        meanCorrTCue_LasttoCurTr.meanGoodNZ(n) = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGoodNonZeroTrCue;
        meanCorrTCue_LasttoCurTr.nGoodNonZeroTr(n) = nGoodNonZeroTrCue;
        
        badTrCorr = triu(spikeCorrTRun_LasttoCurTr{n}(indBadTr,indBadTr),1);
        badTrCorr = badTrCorr(:);
        meanCorrTRun_LasttoCurTr.meanBad(n) = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBad;
        indBadNonZeroTrRun = indBadBeh == 1 & nonZeroTrRun_LasttoCurTr{n} == 1;
        nBadNonZeroTrRun = sum(indBadNonZeroTrRun);
        nElemBadNonZeroTrRun = (nBadNonZeroTrRun*nBadNonZeroTrRun-nBadNonZeroTrRun)/2;
        meanCorrTRun_LasttoCurTr.meanBadNZ(n) = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBadNonZeroTrRun;
        meanCorrTRun_LasttoCurTr.nBadNonZeroTr(n) = nBadNonZeroTrRun;
        
        badTrCorr = triu(spikeCorrTRew_LasttoCurTr{n}(indBadTr,indBadTr),1);
        badTrCorr = badTrCorr(:);
        meanCorrTRew_LasttoCurTr.meanBad(n) = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBad;
        indBadNonZeroTrRew = indBadBeh == 1 & nonZeroTrRew_LasttoCurTr{n} == 1;
        nBadNonZeroTrRew = sum(indBadNonZeroTrRew);
        nElemBadNonZeroTrRew = (nBadNonZeroTrRew*nBadNonZeroTrRew-nBadNonZeroTrRew)/2;
        meanCorrTRew_LasttoCurTr.meanBadNZ(n) = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBadNonZeroTrRew;
        meanCorrTRew_LasttoCurTr.nBadNonZeroTr(n) = nBadNonZeroTrRew;
        
        badTrCorr = triu(spikeCorrTCue_LasttoCurTr{n}(indBadTr,indBadTr),1);
        badTrCorr = badTrCorr(:);
        meanCorrTCue_LasttoCurTr.meanBad(n) = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBad;
        indBadNonZeroTrCue = indBadBeh == 1 & nonZeroTrCue_LasttoCurTr{n} == 1;
        nBadNonZeroTrCue = sum(indBadNonZeroTrCue);
        nElemBadNonZeroTrCue = (nBadNonZeroTrCue*nBadNonZeroTrCue-nBadNonZeroTrCue)/2;
        meanCorrTCue_LasttoCurTr.meanBadNZ(n) = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBadNonZeroTrCue;
        meanCorrTCue_LasttoCurTr.nBadNonZeroTr(n) = nBadNonZeroTrCue;
    end
    
    % compare single neuron correlation between good and bad trials --
    % aligned to run
    plotCompCorr(meanCorrTRun_LasttoCurTr.meanGood,meanCorrTRun_LasttoCurTr.meanBad,...
        'Neu Tr CorrT Run - good trials','Neu Tr CorrT Run - bad trials');
    
%     % compare single neuron correlation between good and bad trials --
%     % aligned to reward
%     plotCompCorr(meanCorrTRew_LasttoCurTr.meanGood,meanCorrTRew_LasttoCurTr.meanBad,...
%         'Neu Tr CorrT Rew - good trials','Neu Tr CorrT Rew - bad trials');
%         
%     % compare single neuron correlation between good and bad trials --
%     % aligned to cue
%     plotCompCorr(meanCorrTCue_LasttoCurTr.meanGood,meanCorrTCue_LasttoCurTr.meanBad,...
%         'Neu Tr CorrT Cue - good trials','Neu Tr CorrT Cue - bad trials');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward
%     plotCompCorr(meanCorrTRun_LasttoCurTr.mean,meanCorrTRew_LasttoCurTr.mean,...
%         'Neu Tr CorrT Run','Neu Tr CorrT Rew');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue
%     plotCompCorr(meanCorrTRun_LasttoCurTr.mean,meanCorrTCue_LasttoCurTr.mean,...
%         'Neu Tr CorrT Run','Neu Tr CorrT Cue');
%     
    % compare single neuron correlation between aligned to run and aligned
    % to reward for good trials
    plotCompCorr(meanCorrTRun_LasttoCurTr.meanGood,meanCorrTRew_LasttoCurTr.meanGood,...
        'Neu Tr CorrT Run - good trials','Neu Tr CorrT Rew - good trials');
   
    % compare single neuron correlation between aligned to run and aligned
    % to cue for good trials
    plotCompCorr(meanCorrTRun_LasttoCurTr.meanGood,meanCorrTCue_LasttoCurTr.meanGood,...
        'Neu Tr CorrT Run - good trials','Neu Tr CorrT Cue - good trials');
%     
%     %% only consider nonzero trials
%     % compare single neuron correlation between good and bad trials --
%     % aligned to run
%     indNeu = meanCorrTRun_LasttoCurTr.nGoodNonZeroTr > 10 & meanCorrTRun_LasttoCurTr.nBadNonZeroTr > 5;
%     plotCompCorr(meanCorrTRun_LasttoCurTr.meanGoodNZ(indNeu),meanCorrTRun_LasttoCurTr.meanBadNZ(indNeu),...
%         'Neu None-zero Tr CorrT Run - good trials','Neu None-zero Tr CorrT Run - bad trials');
%     
%     % compare single neuron correlation between good and bad trials --
%     % aligned to reward
%     indNeu = meanCorrTRew_LasttoCurTr.nGoodNonZeroTr > 10 & meanCorrTRew_LasttoCurTr.nBadNonZeroTr > 5;
%     plotCompCorr(meanCorrTRew_LasttoCurTr.meanGoodNZ(indNeu),meanCorrTRew_LasttoCurTr.meanBadNZ(indNeu),...
%         'Neu None-zero Tr CorrT Rew - good trials','Neu None-zero Tr CorrT Rew - bad trials');
%         
%     % compare single neuron correlation between good and bad trials --
%     % aligned to cue
%     indNeu = meanCorrTCue_LasttoCurTr.nGoodNonZeroTr > 10 & meanCorrTCue_LasttoCurTr.nBadNonZeroTr > 5;
%     plotCompCorr(meanCorrTCue_LasttoCurTr.meanGoodNZ(indNeu),meanCorrTCue_LasttoCurTr.meanBadNZ(indNeu),...
%         'Neu None-zero Tr CorrT Cue - good trials','Neu None-zero Tr CorrT Cue - bad trials');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward
%     indNeu = meanCorrTRun_LasttoCurTr.nNonZeroTr > 10 & meanCorrTRew_LasttoCurTr.nNonZeroTr > 10;
%     plotCompCorr(meanCorrTRun_LasttoCurTr.meanNZ(indNeu),meanCorrTRew_LasttoCurTr.meanNZ(indNeu),...
%         'Neu None-zero Tr CorrT Run ','Neu None-zero Tr CorrT Rew');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue
%     indNeu = meanCorrTRun_LasttoCurTr.nNonZeroTr > 10 & meanCorrTCue_LasttoCurTr.nNonZeroTr > 10;
%     plotCompCorr(meanCorrTRun_LasttoCurTr.meanNZ(indNeu),meanCorrTCue_LasttoCurTr.meanNZ(indNeu),...
%         'Neu None-zero Tr CorrT Run','Neu None-zero Tr CorrT Cue');
%     
    % compare single neuron correlation between aligned to run and aligned
    % to reward for good trials
    indNeu = meanCorrTRun_LasttoCurTr.nGoodNonZeroTr > 5 & meanCorrTRew_LasttoCurTr.nGoodNonZeroTr > 5;
    if(sum(indNeu) > 0)
        plotCompCorr(meanCorrTRun_LasttoCurTr.meanGoodNZ(indNeu),meanCorrTRew_LasttoCurTr.meanGoodNZ(indNeu),...
            'Neu None-zero Tr CorrT Run - good trials','Neu None-zero Tr CorrT Rew - good trials');
    end
   
    % compare single neuron correlation between aligned to run and aligned
    % to cue for good trials
    indNeu = meanCorrTRun_LasttoCurTr.nGoodNonZeroTr > 5 & meanCorrTCue_LasttoCurTr.nGoodNonZeroTr > 5;
    if(sum(indNeu) > 0)
        plotCompCorr(meanCorrTRun_LasttoCurTr.meanGoodNZ(indNeu),meanCorrTCue_LasttoCurTr.meanGoodNZ(indNeu),...
            'Neu None-zero Tr CorrT Run - good trials','Neu None-zero Tr CorrT Cue - good trials');
    end
    
    save([path fileNameCorr],'meanCorrTRun_LasttoCurTr','meanCorrTRew_LasttoCurTr','meanCorrTCue_LasttoCurTr');
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
