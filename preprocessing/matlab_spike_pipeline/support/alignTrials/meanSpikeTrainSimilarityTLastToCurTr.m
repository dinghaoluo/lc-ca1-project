function meanSpikeTrainSimilarityTLastToCurTr(path,fileName,onlyRun,mazeSess,intervalT,intervalTMin)
% single neuron level mean spike correlation across trials
% meanSpikeTrainSimilarityT('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,8)

    fullPath = [path fileName '_spikeTrainSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _spikeTrainSimT_Run file does not exist');
        return;
    end
    load(fullPath);
    
    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    fileNameCorr = [fileName '_meanSpikeTrainSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalTMin) '-' num2str(intervalT) '.mat'];
    
    neuronNo = length(spikeTrainSimTRun_LasttoCurTr);
    trialNo = size(spikeTrainSimTRun_LasttoCurTr{1},1);
    nElem = trialNo*trialNo-trialNo;
    
    indBadBeh = behPar.indTrBadBeh;
    indGoodTr = indBadBeh == 0;
    indBadTr = indBadBeh == 1;
    nGoodTr = sum(indGoodTr);
    nBadTr = sum(indBadTr);
    nElemGood = nGoodTr*nGoodTr-nGoodTr;
    nElemBad = nBadTr*nBadTr-nBadTr;
    
    for n = 1:neuronNo
        simArr = triu(spikeTrainSimTRun_LasttoCurTr{n},1);
        simArr = simArr(:);
        meanSimTRun_LasttoCurTr.mean(n) = sum(simArr(isnan(simArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrRun_LasttoCurTr{n} == 1);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimTRun_LasttoCurTr.meanNZ(n) = sum(simArr(isnan(simArr) == 0))/nElemNonZeroRun;
        meanSimTRun_LasttoCurTr.nNonZeroTr(n) = nNonZeroTr;
                
        simArr = triu(spikeTrainSimTRew_LasttoCurTr{n},1);
        simArr = simArr(:);
        meanSimTRew_LasttoCurTr.mean(n) = sum(simArr(isnan(simArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrRew_LasttoCurTr{n} == 1);
        nElemNonZeroRew = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimTRew_LasttoCurTr.meanNZ(n) = sum(simArr(isnan(simArr) == 0))/nElemNonZeroRew;
        meanSimTRew_LasttoCurTr.nNonZeroTr(n) = nNonZeroTr;
        
        simArr = triu(spikeTrainSimTCue_LasttoCurTr{n},1);
        simArr = simArr(:);
        meanSimTCue_LasttoCurTr.mean(n) = sum(simArr(isnan(simArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrCue_LasttoCurTr{n} == 1);
        nElemNonZeroCue = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimTCue_LasttoCurTr.meanNZ(n) = sum(simArr(isnan(simArr) == 0))/nElemNonZeroCue;
        meanSimTCue_LasttoCurTr.nNonZeroTr(n) = nNonZeroTr;
        
        goodTrSim = spikeTrainSimTRun_LasttoCurTr{n}(indGoodTr,indGoodTr);
        goodTrSim = triu(goodTrSim,1);
        goodTrSim = goodTrSim(:);
        meanSimTRun_LasttoCurTr.meanGood(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGood;
        indGoodNonZeroTrRun = indBadBeh == 0 & nonZeroTrRun_LasttoCurTr{n} == 1;
        nGoodNonZeroTrRun = sum(indGoodNonZeroTrRun);
        nElemGoodNonZeroTrRun = (nGoodNonZeroTrRun*nGoodNonZeroTrRun-nGoodNonZeroTrRun)/2;
        meanSimTRun_LasttoCurTr.meanGoodNZ(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGoodNonZeroTrRun;
        meanSimTRun_LasttoCurTr.nGoodNonZeroTr(n) = nGoodNonZeroTrRun;
        
        goodTrSim = spikeTrainSimTRew_LasttoCurTr{n}(indGoodTr,indGoodTr);
        goodTrSim = triu(goodTrSim,1);
        goodTrSim = goodTrSim(:);
        meanSimTRew_LasttoCurTr.meanGood(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGood;
        indGoodNonZeroTrRew = indBadBeh == 0 & nonZeroTrRew_LasttoCurTr{n} == 1;
        nGoodNonZeroTrRew = sum(indGoodNonZeroTrRew);
        nElemGoodNonZeroTrRew = (nGoodNonZeroTrRew*nGoodNonZeroTrRew-nGoodNonZeroTrRew)/2;
        meanSimTRew_LasttoCurTr.meanGoodNZ(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGoodNonZeroTrRew;
        meanSimTRew_LasttoCurTr.nGoodNonZeroTr(n) = nGoodNonZeroTrRew;
        
        goodTrSim = spikeTrainSimTCue_LasttoCurTr{n}(indGoodTr,indGoodTr);
        goodTrSim = triu(goodTrSim,1);
        goodTrSim = goodTrSim(:);
        meanSimTCue_LasttoCurTr.meanGood(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGood;
        indGoodNonZeroTrCue = indBadBeh == 0 & nonZeroTrCue_LasttoCurTr{n} == 1;
        nGoodNonZeroTrCue = sum(indGoodNonZeroTrCue);
        nElemGoodNonZeroTrCue = (nGoodNonZeroTrCue*nGoodNonZeroTrCue-nGoodNonZeroTrCue)/2;
        meanSimTCue_LasttoCurTr.meanGoodNZ(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGoodNonZeroTrCue;
        meanSimTCue_LasttoCurTr.nGoodNonZeroTr(n) = nGoodNonZeroTrCue;
        
        badTrSim = spikeTrainSimTRun_LasttoCurTr{n}(indBadTr,indBadTr);
        badTrSim = triu(badTrSim,1);
        badTrSim = badTrSim(:);
        meanSimTRun_LasttoCurTr.meanBad(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBad;
        indBadNonZeroTrRun = indBadBeh == 1 & nonZeroTrRun_LasttoCurTr{n} == 1;
        nBadNonZeroTrRun = sum(indBadNonZeroTrRun);
        nElemBadNonZeroTrRun = (nBadNonZeroTrRun*nBadNonZeroTrRun-nBadNonZeroTrRun)/2;
        meanSimTRun_LasttoCurTr.meanBadNZ(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBadNonZeroTrRun;
        meanSimTRun_LasttoCurTr.nBadNonZeroTr(n) = nBadNonZeroTrRun;
        
        badTrSim = spikeTrainSimTRew_LasttoCurTr{n}(indBadTr,indBadTr);
        badTrSim = triu(badTrSim,1);
        badTrSim = badTrSim(:);
        meanSimTRew_LasttoCurTr.meanBad(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBad;
        indBadNonZeroTrRew = indBadBeh == 1 & nonZeroTrRew_LasttoCurTr{n} == 1;
        nBadNonZeroTrRew = sum(indBadNonZeroTrRew);
        nElemBadNonZeroTrRew = (nBadNonZeroTrRew*nBadNonZeroTrRew-nBadNonZeroTrRew)/2;
        meanSimTRew_LasttoCurTr.meanBadNZ(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBadNonZeroTrRew;
        meanSimTRew_LasttoCurTr.nBadNonZeroTr(n) = nBadNonZeroTrRew;
        
        badTrSim = spikeTrainSimTCue_LasttoCurTr{n}(indBadTr,indBadTr);
        badTrSim = triu(badTrSim,1);
        badTrSim = badTrSim(:);
        meanSimTCue_LasttoCurTr.meanBad(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBad;
        indBadNonZeroTrCue = indBadBeh == 1 & nonZeroTrCue_LasttoCurTr{n} == 1;
        nBadNonZeroTrCue = sum(indBadNonZeroTrCue);
        nElemBadNonZeroTrCue = (nBadNonZeroTrCue*nBadNonZeroTrCue-nBadNonZeroTrCue)/2;
        meanSimTCue_LasttoCurTr.meanBadNZ(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBadNonZeroTrCue;
        meanSimTCue_LasttoCurTr.nBadNonZeroTr(n) = nBadNonZeroTrCue;
    end
    
%     % compare single neuron correlation between good and bad trials --
%     % aligned to run
    plotCompCorr(meanSimTRun_LasttoCurTr.meanGood,meanSimTRun_LasttoCurTr.meanBad,...
        'Neu Tr SimT Run - good trials','Neu Tr SimT Run - bad trials');
%     
%     % compare single neuron correlation between good and bad trials --
%     % aligned to reward
%     plotCompCorr(meanSimTRew_LasttoCurTr.meanGood,meanSimTRew_LasttoCurTr.meanBad,...
%         'Neu Tr SimT Rew - good trials','Neu Tr SimT Rew - bad trials');
%         
%     % compare single neuron correlation between good and bad trials --
%     % aligned to cue
%     plotCompCorr(meanSimTCue_LasttoCurTr.meanGood,meanSimTCue_LasttoCurTr.meanBad,...
%         'Neu Tr SimT Cue - good trials','Neu Tr SimT Cue - bad trials');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward
%     plotCompCorr(meanSimTRun_LasttoCurTr.mean,meanSimTRew_LasttoCurTr.mean,...
%         'Neu Tr SimT Run','Neu Tr SimT Rew');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue
%     plotCompCorr(meanSimTRun_LasttoCurTr.mean,meanSimTCue_LasttoCurTr.mean,...
%         'Neu Tr SimT Run','Neu Tr SimT Cue');
%     
    % compare single neuron correlation between aligned to run and aligned
    % to reward for good trials
    plotCompCorr(meanSimTRun_LasttoCurTr.meanGood,meanSimTRew_LasttoCurTr.meanGood,...
        'Neu Tr SimT Run - good trials','Neu Tr SimT Rew - good trials');
   
    % compare single neuron correlation between aligned to run and aligned
    % to cue for good trials
    plotCompCorr(meanSimTRun_LasttoCurTr.meanGood,meanSimTCue_LasttoCurTr.meanGood,...
        'Neu Tr SimT Run - good trials','Neu Tr SimT Cue - good trials');
    
    %% only consider nonzero trials
    % compare single neuron correlation between good and bad trials --
    % aligned to run
%     indNeu = meanSimTRun_LasttoCurTr.nGoodNonZeroTr > 10 & meanSimTRun_LasttoCurTr.nBadNonZeroTr > 5;
%     plotCompCorr(meanSimTRun_LasttoCurTr.meanGoodNZ(indNeu),meanSimTRun_LasttoCurTr.meanBadNZ(indNeu),...
%         'Neu None-zero Tr SimT Run - good trials','Neu None-zero Tr SimT Run - bad trials');
    
%     % compare single neuron correlation between good and bad trials --
%     % aligned to reward
%     indNeu = meanSimTRew_LasttoCurTr.nGoodNonZeroTr > 10 & meanSimTRew_LasttoCurTr.nBadNonZeroTr > 10;
%     plotCompCorr(meanSimTRew_LasttoCurTr.meanGoodNZ(indNeu),meanSimTRew_LasttoCurTr.meanBadNZ(indNeu),...
%         'Neu None-zero Tr SimT Rew - good trials','Neu None-zero Tr SimT Rew - bad trials');
%         
%     % compare single neuron correlation between good and bad trials --
%     % aligned to cue
%     indNeu = meanSimTCue_LasttoCurTr.nGoodNonZeroTr > 10 & meanSimTCue_LasttoCurTr.nBadNonZeroTr > 10;
%     plotCompCorr(meanSimTCue_LasttoCurTr.meanGoodNZ(indNeu),meanSimTCue_LasttoCurTr.meanBadNZ(indNeu),...
%         'Neu None-zero Tr SimT Cue - good trials','Neu None-zero Tr SimT Cue - bad trials');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward
%     indNeu = meanSimTRun_LasttoCurTr.nNonZeroTr > 10 & meanSimTRew_LasttoCurTr.nNonZeroTr > 10;
%     plotCompCorr(meanSimTRun_LasttoCurTr.meanNZ(indNeu),meanSimTRew_LasttoCurTr.meanNZ(indNeu),...
%         'Neu None-zero Tr SimT Run ','Neu None-zero Tr SimT Rew');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue
%     indNeu = meanSimTRun_LasttoCurTr.nNonZeroTr > 10 & meanSimTCue_LasttoCurTr.nNonZeroTr > 10;
%     plotCompCorr(meanSimTRun_LasttoCurTr.meanNZ(indNeu),meanSimTCue_LasttoCurTr.meanNZ(indNeu),...
%         'Neu None-zero Tr SimT Run','Neu None-zero Tr SimT Cue');
    
    % compare single neuron correlation between aligned to run and aligned
    % to reward for good trials
    indNeu = meanSimTRun_LasttoCurTr.nGoodNonZeroTr > 5 & meanSimTRew_LasttoCurTr.nGoodNonZeroTr > 5;
    plotCompCorr(meanSimTRun_LasttoCurTr.meanGoodNZ(indNeu),meanSimTRew_LasttoCurTr.meanGoodNZ(indNeu),...
        'Neu None-zero Tr SimT Run - good trials','Neu None-zero Tr SimT Rew - good trials');
   
    % compare single neuron correlation between aligned to run and aligned
    % to cue for good trials
    indNeu = meanSimTRun_LasttoCurTr.nGoodNonZeroTr > 5 & meanSimTCue_LasttoCurTr.nGoodNonZeroTr > 5;
    plotCompCorr(meanSimTRun_LasttoCurTr.meanGoodNZ(indNeu),meanSimTCue_LasttoCurTr.meanGoodNZ(indNeu),...
        'Neu None-zero Tr SimT Run - good trials','Neu None-zero Tr SimT Cue - good trials');
    
    save([path fileNameCorr],'meanSimTRun_LasttoCurTr','meanSimTRew_LasttoCurTr','meanSimTCue_LasttoCurTr');
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
