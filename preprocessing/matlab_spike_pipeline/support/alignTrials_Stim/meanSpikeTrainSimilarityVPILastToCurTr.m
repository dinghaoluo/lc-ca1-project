function meanSpikeTrainSimilarityVPILastToCurTr(path,fileName,onlyRun,mazeSess,cost,intervalT,intervalTMin)
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
    load(fullPath,'spikeTrainSimVPIRun_LasttoCurTr','spikeTrainSimVPIRew_LasttoCurTr','spikeTrainSimVPICue_LasttoCurTr',...
        'nonZeroTrRun_LasttoCurTr','nonZeroTrRew_LasttoCurTr','nonZeroTrCue_LasttoCurTr');
    
    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    fileNameCorr = [fileName '_meanSpikeTrainSimVP_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_q' ...
            s '_intT'  num2str(intervalTMin) '-' num2str(intervalT) '.mat'];
    
    neuronNo = length(spikeTrainSimVPIRun_LasttoCurTr);
    trialNo = size(spikeTrainSimVPIRun_LasttoCurTr{1},1);
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
        simArr = triu(spikeTrainSimVPIRun_LasttoCurTr{n},1);
        simArr = abs(simArr(:));
        meanSimVPIRun_LasttoCurTr.mean(n) = sum(simArr(isnan(simArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrRun_LasttoCurTr{n} == 1);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPIRun_LasttoCurTr.meanNZ(n) = sum(simArr(isnan(simArr) == 0))/nElemNonZeroRun;
        meanSimVPIRun_LasttoCurTr.nNonZeroTr(n) = nNonZeroTr;
                
        simArr = triu(spikeTrainSimVPIRew_LasttoCurTr{n},1);
        simArr = abs(simArr(:));
        meanSimVPIRew_LasttoCurTr.mean(n) = sum(simArr(isnan(simArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrRew_LasttoCurTr{n} == 1);
        nElemNonZeroRew = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPIRew_LasttoCurTr.meanNZ(n) = sum(simArr(isnan(simArr) == 0))/nElemNonZeroRew;
        meanSimVPIRew_LasttoCurTr.nNonZeroTr(n) = nNonZeroTr;
        
        simArr = triu(spikeTrainSimVPICue_LasttoCurTr{n},1);
        simArr = abs(simArr(:));
        meanSimVPICue_LasttoCurTr.mean(n) = sum(simArr(isnan(simArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrCue_LasttoCurTr{n} == 1);
        nElemNonZeroCue = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPICue_LasttoCurTr.meanNZ(n) = sum(simArr(isnan(simArr) == 0))/nElemNonZeroCue;
        meanSimVPICue_LasttoCurTr.nNonZeroTr(n) = nNonZeroTr;
        
        %% non zero elements in good trials
        goodTrSim = spikeTrainSimVPIRun_LasttoCurTr{n}(indGoodTr,indGoodTr);
        goodTrSim = triu(goodTrSim,1);
        goodTrSim = abs(goodTrSim(:));
        meanSimVPIRun_LasttoCurTr.meanGood(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGood;
        indGoodNonZeroTrRun = indBadBeh == 0 & nonZeroTrRun_LasttoCurTr{n} == 1;
        nGoodNonZeroTrRun = sum(indGoodNonZeroTrRun);
        nElemGoodNonZeroTrRun = (nGoodNonZeroTrRun*nGoodNonZeroTrRun-nGoodNonZeroTrRun)/2;
        meanSimVPIRun_LasttoCurTr.meanGoodNZ(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGoodNonZeroTrRun;
        meanSimVPIRun_LasttoCurTr.nGoodNonZeroTr(n) = nGoodNonZeroTrRun;
        
        goodTrSim = spikeTrainSimVPIRew_LasttoCurTr{n}(indGoodTr,indGoodTr);
        goodTrSim = triu(goodTrSim,1);
        goodTrSim = abs(goodTrSim(:));
        meanSimVPIRew_LasttoCurTr.meanGood(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGood;
        indGoodNonZeroTrRew = indBadBeh == 0 & nonZeroTrRew_LasttoCurTr{n} == 1;
        nGoodNonZeroTrRew = sum(indGoodNonZeroTrRew);
        nElemGoodNonZeroTrRew = (nGoodNonZeroTrRew*nGoodNonZeroTrRew-nGoodNonZeroTrRew)/2;
        meanSimVPIRew_LasttoCurTr.meanGoodNZ(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGoodNonZeroTrRew;
        meanSimVPIRew_LasttoCurTr.nGoodNonZeroTr(n) = nGoodNonZeroTrRew;
        
        goodTrSim = spikeTrainSimVPICue_LasttoCurTr{n}(indGoodTr,indGoodTr);
        goodTrSim = triu(goodTrSim,1);
        goodTrSim = abs(goodTrSim(:));
        meanSimVPICue_LasttoCurTr.meanGood(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGood;
        indGoodNonZeroTrCue = indBadBeh == 0 & nonZeroTrCue_LasttoCurTr{n} == 1;
        nGoodNonZeroTrCue = sum(indGoodNonZeroTrCue);
        nElemGoodNonZeroTrCue = (nGoodNonZeroTrCue*nGoodNonZeroTrCue-nGoodNonZeroTrCue)/2;
        meanSimVPICue_LasttoCurTr.meanGoodNZ(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGoodNonZeroTrCue;
        meanSimVPICue_LasttoCurTr.nGoodNonZeroTr(n) = nGoodNonZeroTrCue;
        
        %% non zero elements in bad trials
        badTrSim = spikeTrainSimVPIRun_LasttoCurTr{n}(indBadTr,indBadTr);
        badTrSim = triu(badTrSim,1);
        badTrSim = abs(badTrSim(:));
        meanSimVPIRun_LasttoCurTr.meanBad(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBad;
        indBadNonZeroTrRun = indBadBeh == 1 & nonZeroTrRun_LasttoCurTr{n} == 1;
        nBadNonZeroTrRun = sum(indBadNonZeroTrRun);
        nElemBadNonZeroTrRun = (nBadNonZeroTrRun*nBadNonZeroTrRun-nBadNonZeroTrRun)/2;
        meanSimVPIRun_LasttoCurTr.meanBadNZ(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBadNonZeroTrRun;
        meanSimVPIRun_LasttoCurTr.nBadNonZeroTr(n) = nBadNonZeroTrRun;
        
        badTrSim = spikeTrainSimVPIRew_LasttoCurTr{n}(indBadTr,indBadTr);
        badTrSim = triu(badTrSim,1);
        badTrSim = abs(badTrSim(:));
        meanSimVPIRew_LasttoCurTr.meanBad(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBad;
        indBadNonZeroTrRew = indBadBeh == 1 & nonZeroTrRew_LasttoCurTr{n} == 1;
        nBadNonZeroTrRew = sum(indBadNonZeroTrRew);
        nElemBadNonZeroTrRew = (nBadNonZeroTrRew*nBadNonZeroTrRew-nBadNonZeroTrRew)/2;
        meanSimVPIRew_LasttoCurTr.meanBadNZ(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBadNonZeroTrRew;
        meanSimVPIRew_LasttoCurTr.nBadNonZeroTr(n) = nBadNonZeroTrRew;
        
        badTrSim = spikeTrainSimVPICue_LasttoCurTr{n}(indBadTr,indBadTr);
        badTrSim = triu(badTrSim,1);
        badTrSim = abs(badTrSim(:));
        meanSimVPICue_LasttoCurTr.meanBad(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBad;
        indBadNonZeroTrCue = indBadBeh == 1 & nonZeroTrCue_LasttoCurTr{n} == 1;
        nBadNonZeroTrCue = sum(indBadNonZeroTrCue);
        nElemBadNonZeroTrCue = (nBadNonZeroTrCue*nBadNonZeroTrCue-nBadNonZeroTrCue)/2;
        meanSimVPICue_LasttoCurTr.meanBadNZ(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBadNonZeroTrCue;
        meanSimVPICue_LasttoCurTr.nBadNonZeroTr(n) = nBadNonZeroTrCue;
    end
    
    % compare single neuron correlation between good and bad trials --
    % aligned to run
    plotCompCorr(meanSimVPIRun_LasttoCurTr.meanGood,meanSimVPIRun_LasttoCurTr.meanBad,...
        'Neu Tr SimVPI Run - good trials','Neu Tr SimVPI Run - bad trials');
    
%     % compare single neuron correlation between good and bad trials --
%     % aligned to reward
%     plotCompCorr(meanSimVPIRew_LasttoCurTr.meanGood,meanSimVPIRew_LasttoCurTr.meanBad,...
%         'Neu Tr SimVPI Rew - good trials','Neu Tr SimVPI Rew - bad trials');
%         
%     % compare single neuron correlation between good and bad trials --
%     % aligned to cue
%     plotCompCorr(meanSimVPICue_LasttoCurTr.meanGood,meanSimVPICue_LasttoCurTr.meanBad,...
%         'Neu Tr SimVPI Cue - good trials','Neu Tr SimVPI Cue - bad trials');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward
%     plotCompCorr(meanSimVPIRun_LasttoCurTr.mean,meanSimVPIRew_LasttoCurTr.mean,...
%         'Neu Tr SimVPI Run','Neu Tr SimVPI Rew');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue
%     plotCompCorr(meanSimVPIRun_LasttoCurTr.mean,meanSimVPICue_LasttoCurTr.mean,...
%         'Neu Tr SimVPI Run','Neu Tr SimVPI Cue');
%     
%     % compare single neuron correlation between aligned to run and aligned
    % to reward for good trials
    plotCompCorr(meanSimVPIRun_LasttoCurTr.meanGood,meanSimVPIRew_LasttoCurTr.meanGood,...
        'Neu Tr SimVPI Run - good trials','Neu Tr SimVPI Rew - good trials');
   
    % compare single neuron correlation between aligned to run and aligned
    % to cue for good trials
    plotCompCorr(meanSimVPIRun_LasttoCurTr.meanGood,meanSimVPICue_LasttoCurTr.meanGood,...
        'Neu Tr SimVPI Run - good trials','Neu Tr SimVPI Cue - good trials');
    
%     %% only consider nonzero trials
%     % compare single neuron correlation between good and bad trials --
%     % aligned to run
%     indNeu = meanSimVPIRun_LasttoCurTr.nGoodNonZeroTr > 10 & meanSimVPIRun_LasttoCurTr.nBadNonZeroTr > 10;
%     plotCompCorr(meanSimVPIRun_LasttoCurTr.meanGoodNZ(indNeu),meanSimVPIRun_LasttoCurTr.meanBadNZ(indNeu),...
%         'Neu None-zero Tr SimVPI Run - good trials','Neu None-zero Tr SimVPI Run - bad trials');
%     
%     % compare single neuron correlation between good and bad trials --
%     % aligned to reward
%     indNeu = meanSimVPIRew_LasttoCurTr.nGoodNonZeroTr > 10 & meanSimVPIRew_LasttoCurTr.nBadNonZeroTr > 10;
%     plotCompCorr(meanSimVPIRew_LasttoCurTr.meanGoodNZ(indNeu),meanSimVPIRew_LasttoCurTr.meanBadNZ(indNeu),...
%         'Neu None-zero Tr SimVPI Rew - good trials','Neu None-zero Tr SimVPI Rew - bad trials');
%         
%     % compare single neuron correlation between good and bad trials --
%     % aligned to cue
%     indNeu = meanSimVPICue_LasttoCurTr.nGoodNonZeroTr > 10 & meanSimVPICue_LasttoCurTr.nBadNonZeroTr > 10;
%     plotCompCorr(meanSimVPICue_LasttoCurTr.meanGoodNZ(indNeu),meanSimVPICue_LasttoCurTr.meanBadNZ(indNeu),...
%         'Neu None-zero Tr SimVPI Cue - good trials','Neu None-zero Tr SimVPI Cue - bad trials');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward
%     indNeu = meanSimVPIRun_LasttoCurTr.nNonZeroTr > 10 & meanSimVPIRew_LasttoCurTr.nNonZeroTr > 10;
%     plotCompCorr(meanSimVPIRun_LasttoCurTr.meanNZ(indNeu),meanSimVPIRew_LasttoCurTr.meanNZ(indNeu),...
%         'Neu None-zero Tr SimVPI Run ','Neu None-zero Tr SimVPI Rew');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue
%     indNeu = meanSimVPIRun_LasttoCurTr.nNonZeroTr > 10 & meanSimVPICue_LasttoCurTr.nNonZeroTr > 10;
%     plotCompCorr(meanSimVPIRun_LasttoCurTr.meanNZ(indNeu),meanSimVPICue_LasttoCurTr.meanNZ(indNeu),...
%         'Neu None-zero Tr SimVPI Run','Neu None-zero Tr SimVPI Cue');
%     
    % compare single neuron correlation between aligned to run and aligned
    % to reward for good trials
    indNeu = meanSimVPIRun_LasttoCurTr.nGoodNonZeroTr > 5 & meanSimVPIRew_LasttoCurTr.nGoodNonZeroTr > 5;
    plotCompCorr(meanSimVPIRun_LasttoCurTr.meanGoodNZ(indNeu),meanSimVPIRew_LasttoCurTr.meanGoodNZ(indNeu),...
        'Neu None-zero Tr SimVPI Run - good trials','Neu None-zero Tr SimVPI Rew - good trials');
   
    % compare single neuron correlation between aligned to run and aligned
    % to cue for good trials
    indNeu = meanSimVPIRun_LasttoCurTr.nGoodNonZeroTr > 5 & meanSimVPICue_LasttoCurTr.nGoodNonZeroTr > 5;
    plotCompCorr(meanSimVPIRun_LasttoCurTr.meanGoodNZ(indNeu),meanSimVPICue_LasttoCurTr.meanGoodNZ(indNeu),...
        'Neu None-zero Tr SimVPI Run - good trials','Neu None-zero Tr SimVPI Cue - good trials');
    
    save([path fileNameCorr],'meanSimVPIRun_LasttoCurTr','meanSimVPIRew_LasttoCurTr','meanSimVPICue_LasttoCurTr','-append');
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
