function meanSpikeTrainSimilarityVPI(path,fileName,onlyRun,mazeSess,cost,intervalT)
% single neuron level mean spike correlation across trials
% Victor & Purpura spike time interval
% e.g.: meanSpikeTrainSimilarityVPI('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,100,8)

    s = num2str(cost);
    ind = strfind(s,'.');
    s(ind) = 'p';
    fullPath = [path fileName '_spikeTrainSimVP_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_q' ...
            s '_intT' num2str(intervalT) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _spikeTrainSimVP_Run file does not exist');
        return;
    end
    load(fullPath,'spikeTrainSimVPIRun','spikeTrainSimVPIRew','spikeTrainSimVPICue',...
        'nonZeroTrRun','nonZeroTrRew','nonZeroTrCue');
    
    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    fileNameCorr = [fileName '_meanSpikeTrainSimVP_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_q' ...
            s '_intT' num2str(intervalT) '.mat'];
    
    neuronNo = length(spikeTrainSimVPIRun);
    trialNo = size(spikeTrainSimVPIRun{1},1);
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
        simArr = triu(spikeTrainSimVPIRun{n},1);
        simArr = abs(simArr(:));
        meanSimVPIRun.mean(n) = sum(simArr(isnan(simArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrRun{n} == 1);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPIRun.meanNZ(n) = sum(simArr(isnan(simArr) == 0))/nElemNonZeroRun;
        meanSimVPIRun.nNonZeroTr(n) = nNonZeroTr;
                
        simArr = triu(spikeTrainSimVPIRew{n},1);
        simArr = abs(simArr(:));
        meanSimVPIRew.mean(n) = sum(simArr(isnan(simArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrRew{n} == 1);
        nElemNonZeroRew = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPIRew.meanNZ(n) = sum(simArr(isnan(simArr) == 0))/nElemNonZeroRew;
        meanSimVPIRew.nNonZeroTr(n) = nNonZeroTr;
        
        simArr = triu(spikeTrainSimVPICue{n},1);
        simArr = abs(simArr(:));
        meanSimVPICue.mean(n) = sum(simArr(isnan(simArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrCue{n} == 1);
        nElemNonZeroCue = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPICue.meanNZ(n) = sum(simArr(isnan(simArr) == 0))/nElemNonZeroCue;
        meanSimVPICue.nNonZeroTr(n) = nNonZeroTr;
        
        %% non zero elements in good trials
        goodTrSim = spikeTrainSimVPIRun{n}(indGoodTr,indGoodTr);
        goodTrSim = triu(goodTrSim,1);
        goodTrSim = abs(goodTrSim(:));
        meanSimVPIRun.meanGood(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGood;
        indGoodNonZeroTrRun = indBadBeh == 0 & nonZeroTrRun{n} == 1;
        nGoodNonZeroTrRun = sum(indGoodNonZeroTrRun);
        nElemGoodNonZeroTrRun = (nGoodNonZeroTrRun*nGoodNonZeroTrRun-nGoodNonZeroTrRun)/2;
        meanSimVPIRun.meanGoodNZ(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGoodNonZeroTrRun;
        meanSimVPIRun.nGoodNonZeroTr(n) = nGoodNonZeroTrRun;
        
        goodTrSim = spikeTrainSimVPIRew{n}(indGoodTr,indGoodTr);
        goodTrSim = triu(goodTrSim,1);
        goodTrSim = abs(goodTrSim(:));
        meanSimVPIRew.meanGood(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGood;
        indGoodNonZeroTrRew = indBadBeh == 0 & nonZeroTrRew{n} == 1;
        nGoodNonZeroTrRew = sum(indGoodNonZeroTrRew);
        nElemGoodNonZeroTrRew = (nGoodNonZeroTrRew*nGoodNonZeroTrRew-nGoodNonZeroTrRew)/2;
        meanSimVPIRew.meanGoodNZ(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGoodNonZeroTrRew;
        meanSimVPIRew.nGoodNonZeroTr(n) = nGoodNonZeroTrRew;
        
        goodTrSim = spikeTrainSimVPICue{n}(indGoodTr,indGoodTr);
        goodTrSim = triu(goodTrSim,1);
        goodTrSim = abs(goodTrSim(:));
        meanSimVPICue.meanGood(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGood;
        indGoodNonZeroTrCue = indBadBeh == 0 & nonZeroTrCue{n} == 1;
        nGoodNonZeroTrCue = sum(indGoodNonZeroTrCue);
        nElemGoodNonZeroTrCue = (nGoodNonZeroTrCue*nGoodNonZeroTrCue-nGoodNonZeroTrCue)/2;
        meanSimVPICue.meanGoodNZ(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGoodNonZeroTrCue;
        meanSimVPICue.nGoodNonZeroTr(n) = nGoodNonZeroTrCue;
        
        %% non zero elements in bad trials
        badTrSim = spikeTrainSimVPIRun{n}(indBadTr,indBadTr);
        badTrSim = triu(badTrSim,1);
        badTrSim = abs(badTrSim(:));
        meanSimVPIRun.meanBad(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBad;
        indBadNonZeroTrRun = indBadBeh == 1 & nonZeroTrRun{n} == 1;
        nBadNonZeroTrRun = sum(indBadNonZeroTrRun);
        nElemBadNonZeroTrRun = (nBadNonZeroTrRun*nBadNonZeroTrRun-nBadNonZeroTrRun)/2;
        meanSimVPIRun.meanBadNZ(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBadNonZeroTrRun;
        meanSimVPIRun.nBadNonZeroTr(n) = nBadNonZeroTrRun;
        
        badTrSim = spikeTrainSimVPIRew{n}(indBadTr,indBadTr);
        badTrSim = triu(badTrSim,1);
        badTrSim = abs(badTrSim(:));
        meanSimVPIRew.meanBad(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBad;
        indBadNonZeroTrRew = indBadBeh == 1 & nonZeroTrRew{n} == 1;
        nBadNonZeroTrRew = sum(indBadNonZeroTrRew);
        nElemBadNonZeroTrRew = (nBadNonZeroTrRew*nBadNonZeroTrRew-nBadNonZeroTrRew)/2;
        meanSimVPIRew.meanBadNZ(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBadNonZeroTrRew;
        meanSimVPIRew.nBadNonZeroTr(n) = nBadNonZeroTrRew;
        
        badTrSim = spikeTrainSimVPICue{n}(indBadTr,indBadTr);
        badTrSim = triu(badTrSim,1);
        badTrSim = abs(badTrSim(:));
        meanSimVPICue.meanBad(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBad;
        indBadNonZeroTrCue = indBadBeh == 1 & nonZeroTrCue{n} == 1;
        nBadNonZeroTrCue = sum(indBadNonZeroTrCue);
        nElemBadNonZeroTrCue = (nBadNonZeroTrCue*nBadNonZeroTrCue-nBadNonZeroTrCue)/2;
        meanSimVPICue.meanBadNZ(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBadNonZeroTrCue;
        meanSimVPICue.nBadNonZeroTr(n) = nBadNonZeroTrCue;
    end
    
    % compare single neuron correlation between good and bad trials --
    % aligned to run
    plotCompCorr(meanSimVPIRun.meanGood,meanSimVPIRun.meanBad,...
        'Neu Tr SimVPI Run - good trials','Neu Tr SimVPI Run - bad trials');
    
%     % compare single neuron correlation between good and bad trials --
%     % aligned to reward
%     plotCompCorr(meanSimVPIRew.meanGood,meanSimVPIRew.meanBad,...
%         'Neu Tr SimVPI Rew - good trials','Neu Tr SimVPI Rew - bad trials');
%         
%     % compare single neuron correlation between good and bad trials --
%     % aligned to cue
%     plotCompCorr(meanSimVPICue.meanGood,meanSimVPICue.meanBad,...
%         'Neu Tr SimVPI Cue - good trials','Neu Tr SimVPI Cue - bad trials');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward
%     plotCompCorr(meanSimVPIRun.mean,meanSimVPIRew.mean,...
%         'Neu Tr SimVPI Run','Neu Tr SimVPI Rew');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue
%     plotCompCorr(meanSimVPIRun.mean,meanSimVPICue.mean,...
%         'Neu Tr SimVPI Run','Neu Tr SimVPI Cue');
%     
%     % compare single neuron correlation between aligned to run and aligned
    % to reward for good trials
    plotCompCorr(meanSimVPIRun.meanGood,meanSimVPIRew.meanGood,...
        'Neu Tr SimVPI Run - good trials','Neu Tr SimVPI Rew - good trials');
   
    % compare single neuron correlation between aligned to run and aligned
    % to cue for good trials
    plotCompCorr(meanSimVPIRun.meanGood,meanSimVPICue.meanGood,...
        'Neu Tr SimVPI Run - good trials','Neu Tr SimVPI Cue - good trials');
    
%     %% only consider nonzero trials
%     % compare single neuron correlation between good and bad trials --
%     % aligned to run
%     indNeu = meanSimVPIRun.nGoodNonZeroTr > 10 & meanSimVPIRun.nBadNonZeroTr > 10;
%     plotCompCorr(meanSimVPIRun.meanGoodNZ(indNeu),meanSimVPIRun.meanBadNZ(indNeu),...
%         'Neu None-zero Tr SimVPI Run - good trials','Neu None-zero Tr SimVPI Run - bad trials');
%     
%     % compare single neuron correlation between good and bad trials --
%     % aligned to reward
%     indNeu = meanSimVPIRew.nGoodNonZeroTr > 10 & meanSimVPIRew.nBadNonZeroTr > 10;
%     plotCompCorr(meanSimVPIRew.meanGoodNZ(indNeu),meanSimVPIRew.meanBadNZ(indNeu),...
%         'Neu None-zero Tr SimVPI Rew - good trials','Neu None-zero Tr SimVPI Rew - bad trials');
%         
%     % compare single neuron correlation between good and bad trials --
%     % aligned to cue
%     indNeu = meanSimVPICue.nGoodNonZeroTr > 10 & meanSimVPICue.nBadNonZeroTr > 10;
%     plotCompCorr(meanSimVPICue.meanGoodNZ(indNeu),meanSimVPICue.meanBadNZ(indNeu),...
%         'Neu None-zero Tr SimVPI Cue - good trials','Neu None-zero Tr SimVPI Cue - bad trials');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward
%     indNeu = meanSimVPIRun.nNonZeroTr > 10 & meanSimVPIRew.nNonZeroTr > 10;
%     plotCompCorr(meanSimVPIRun.meanNZ(indNeu),meanSimVPIRew.meanNZ(indNeu),...
%         'Neu None-zero Tr SimVPI Run ','Neu None-zero Tr SimVPI Rew');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue
%     indNeu = meanSimVPIRun.nNonZeroTr > 10 & meanSimVPICue.nNonZeroTr > 10;
%     plotCompCorr(meanSimVPIRun.meanNZ(indNeu),meanSimVPICue.meanNZ(indNeu),...
%         'Neu None-zero Tr SimVPI Run','Neu None-zero Tr SimVPI Cue');
%     
    % compare single neuron correlation between aligned to run and aligned
    % to reward for good trials
    indNeu = meanSimVPIRun.nGoodNonZeroTr > 5 & meanSimVPIRew.nGoodNonZeroTr > 5;
    plotCompCorr(meanSimVPIRun.meanGoodNZ(indNeu),meanSimVPIRew.meanGoodNZ(indNeu),...
        'Neu None-zero Tr SimVPI Run - good trials','Neu None-zero Tr SimVPI Rew - good trials');
   
    % compare single neuron correlation between aligned to run and aligned
    % to cue for good trials
    indNeu = meanSimVPIRun.nGoodNonZeroTr > 5 & meanSimVPICue.nGoodNonZeroTr > 5;
    plotCompCorr(meanSimVPIRun.meanGoodNZ(indNeu),meanSimVPICue.meanGoodNZ(indNeu),...
        'Neu None-zero Tr SimVPI Run - good trials','Neu None-zero Tr SimVPI Cue - good trials');
    
    save([path fileNameCorr],'meanSimVPIRun','meanSimVPIRew','meanSimVPICue','-append');
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
