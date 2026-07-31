function meanSpikeTrainSimilarityVPLastToCurTr(path,fileName,onlyRun,mazeSess,cost,intervalT,intervalTMin)
% single neuron level mean spike correlation across trials
% Victor & Purpura spike time distance
% e.g.: meanSpikeTrainSimilarityVP('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,0.016,8)

    s = num2str(cost);
    ind = strfind(s,'.');
    s(ind) = 'p';
    fullPath = [path fileName '_spikeTrainSimVP_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_q' ...
            s '_intT' num2str(intervalT) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _spikeTrainSimiVP_Run file does not exist');
        return;
    end
    load(fullPath,'spikeTrainSimVPRun_LasttoCurTr','spikeTrainSimVPRew_LasttoCurTr','spikeTrainSimVPCue_LasttoCurTr',...
        'nonZeroTrRun_LasttoCurTr','nonZeroTrRew_LasttoCurTr','nonZeroTrCue_LasttoCurTr');
    
    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    fileNameCorr = [fileName '_meanSpikeTrainSimVP_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_q' ...
            s '_intT' num2str(intervalTMin) '-' num2str(intervalT) '.mat'];
    
    neuronNo = length(spikeTrainSimVPRun_LasttoCurTr);
    trialNo = size(spikeTrainSimVPRun_LasttoCurTr{1},1);
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
        simArr = triu(spikeTrainSimVPRun_LasttoCurTr{n},1);
        simArr = abs(simArr(:));
        meanSimVPRun_LasttoCurTr.mean(n) = sum(simArr(isnan(simArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrRun_LasttoCurTr{n} == 1);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPRun_LasttoCurTr.meanNZ(n) = sum(simArr(isnan(simArr) == 0))/nElemNonZeroRun;
        meanSimVPRun_LasttoCurTr.nNonZeroTr(n) = nNonZeroTr;
                
        simArr = triu(spikeTrainSimVPRew_LasttoCurTr{n},1);
        simArr = abs(simArr(:));
        meanSimVPRew_LasttoCurTr.mean(n) = sum(simArr(isnan(simArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrRew_LasttoCurTr{n} == 1);
        nElemNonZeroRew = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPRew_LasttoCurTr.meanNZ(n) = sum(simArr(isnan(simArr) == 0))/nElemNonZeroRew;
        meanSimVPRew_LasttoCurTr.nNonZeroTr(n) = nNonZeroTr;
        
        simArr = triu(spikeTrainSimVPCue_LasttoCurTr{n},1);
        simArr = abs(simArr(:));
        meanSimVPCue_LasttoCurTr.mean(n) = sum(simArr(isnan(simArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrCue_LasttoCurTr{n} == 1);
        nElemNonZeroCue = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPCue_LasttoCurTr.meanNZ(n) = sum(simArr(isnan(simArr) == 0))/nElemNonZeroCue;
        meanSimVPCue_LasttoCurTr.nNonZeroTr(n) = nNonZeroTr;
        
        %% non zero elements in good trials
        goodTrSim = spikeTrainSimVPRun_LasttoCurTr{n}(indGoodTr,indGoodTr);
        goodTrSim = triu(goodTrSim,1);
        goodTrSim = abs(goodTrSim(:));
        meanSimVPRun_LasttoCurTr.meanGood(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGood;
        indGoodNonZeroTrRun = indBadBeh == 0 & nonZeroTrRun_LasttoCurTr{n} == 1;
        nGoodNonZeroTrRun = sum(indGoodNonZeroTrRun);
        nElemGoodNonZeroTrRun = (nGoodNonZeroTrRun*nGoodNonZeroTrRun-nGoodNonZeroTrRun)/2;
        meanSimVPRun_LasttoCurTr.meanGoodNZ(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGoodNonZeroTrRun;
        meanSimVPRun_LasttoCurTr.nGoodNonZeroTr(n) = nGoodNonZeroTrRun;
        
        goodTrSim = spikeTrainSimVPRew_LasttoCurTr{n}(indGoodTr,indGoodTr);
        goodTrSim = triu(goodTrSim,1);
        goodTrSim = abs(goodTrSim(:));
        meanSimVPRew_LasttoCurTr.meanGood(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGood;
        indGoodNonZeroTrRew = indBadBeh == 0 & nonZeroTrRew_LasttoCurTr{n} == 1;
        nGoodNonZeroTrRew = sum(indGoodNonZeroTrRew);
        nElemGoodNonZeroTrRew = (nGoodNonZeroTrRew*nGoodNonZeroTrRew-nGoodNonZeroTrRew)/2;
        meanSimVPRew_LasttoCurTr.meanGoodNZ(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGoodNonZeroTrRew;
        meanSimVPRew_LasttoCurTr.nGoodNonZeroTr(n) = nGoodNonZeroTrRew;
        
        goodTrSim = spikeTrainSimVPCue_LasttoCurTr{n}(indGoodTr,indGoodTr);
        goodTrSim = triu(goodTrSim,1);
        goodTrSim = abs(goodTrSim(:));
        meanSimVPCue_LasttoCurTr.meanGood(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGood;
        indGoodNonZeroTrCue = indBadBeh == 0 & nonZeroTrCue_LasttoCurTr{n} == 1;
        nGoodNonZeroTrCue = sum(indGoodNonZeroTrCue);
        nElemGoodNonZeroTrCue = (nGoodNonZeroTrCue*nGoodNonZeroTrCue-nGoodNonZeroTrCue)/2;
        meanSimVPCue_LasttoCurTr.meanGoodNZ(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGoodNonZeroTrCue;
        meanSimVPCue_LasttoCurTr.nGoodNonZeroTr(n) = nGoodNonZeroTrCue;
        
        %% non zero elements in bad trials
        badTrSim = spikeTrainSimVPRun_LasttoCurTr{n}(indBadTr,indBadTr);
        badTrSim = triu(badTrSim,1);
        badTrSim = abs(badTrSim(:));
        meanSimVPRun_LasttoCurTr.meanBad(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBad;
        indBadNonZeroTrRun = indBadBeh == 1 & nonZeroTrRun_LasttoCurTr{n} == 1;
        nBadNonZeroTrRun = sum(indBadNonZeroTrRun);
        nElemBadNonZeroTrRun = (nBadNonZeroTrRun*nBadNonZeroTrRun-nBadNonZeroTrRun)/2;
        meanSimVPRun_LasttoCurTr.meanBadNZ(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBadNonZeroTrRun;
        meanSimVPRun_LasttoCurTr.nBadNonZeroTr(n) = nBadNonZeroTrRun;
        
        badTrSim = spikeTrainSimVPRew_LasttoCurTr{n}(indBadTr,indBadTr);
        badTrSim = triu(badTrSim,1);
        badTrSim = abs(badTrSim(:));
        meanSimVPRew_LasttoCurTr.meanBad(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBad;
        indBadNonZeroTrRew = indBadBeh == 1 & nonZeroTrRew_LasttoCurTr{n} == 1;
        nBadNonZeroTrRew = sum(indBadNonZeroTrRew);
        nElemBadNonZeroTrRew = (nBadNonZeroTrRew*nBadNonZeroTrRew-nBadNonZeroTrRew)/2;
        meanSimVPRew_LasttoCurTr.meanBadNZ(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBadNonZeroTrRew;
        meanSimVPRew_LasttoCurTr.nBadNonZeroTr(n) = nBadNonZeroTrRew;
        
        badTrSim = spikeTrainSimVPCue_LasttoCurTr{n}(indBadTr,indBadTr);
        badTrSim = triu(badTrSim,1);
        badTrSim = abs(badTrSim(:));
        meanSimVPCue_LasttoCurTr.meanBad(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBad;
        indBadNonZeroTrCue = indBadBeh == 1 & nonZeroTrCue_LasttoCurTr{n} == 1;
        nBadNonZeroTrCue = sum(indBadNonZeroTrCue);
        nElemBadNonZeroTrCue = (nBadNonZeroTrCue*nBadNonZeroTrCue-nBadNonZeroTrCue)/2;
        meanSimVPCue_LasttoCurTr.meanBadNZ(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBadNonZeroTrCue;
        meanSimVPCue_LasttoCurTr.nBadNonZeroTr(n) = nBadNonZeroTrCue;
    end
    
    % compare single neuron correlation between good and bad trials --
    % aligned to run
    plotCompCorr(meanSimVPRun_LasttoCurTr.meanGood,meanSimVPRun_LasttoCurTr.meanBad,...
        'Neu Tr SimVP Run - good trials','Neu Tr SimVP Run - bad trials');
    
%     % compare single neuron correlation between good and bad trials --
%     % aligned to reward
%     plotCompCorr(meanSimVPRew_LasttoCurTr.meanGood,meanSimVPRew_LasttoCurTr.meanBad,...
%         'Neu Tr SimVP Rew - good trials','Neu Tr SimVP Rew - bad trials');
%         
%     % compare single neuron correlation between good and bad trials --
%     % aligned to cue
%     plotCompCorr(meanSimVPCue_LasttoCurTr.meanGood,meanSimVPCue_LasttoCurTr.meanBad,...
%         'Neu Tr SimVP Cue - good trials','Neu Tr SimVP Cue - bad trials');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward
%     plotCompCorr(meanSimVPRun_LasttoCurTr.mean,meanSimVPRew_LasttoCurTr.mean,...
%         'Neu Tr SimVP Run','Neu Tr SimVP Rew');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue
%     plotCompCorr(meanSimVPRun_LasttoCurTr.mean,meanSimVPCue_LasttoCurTr.mean,...
%         'Neu Tr SimVP Run','Neu Tr SimVP Cue');
    
    % compare single neuron correlation between aligned to run and aligned
    % to reward for good trials
    plotCompCorr(meanSimVPRun_LasttoCurTr.meanGood,meanSimVPRew_LasttoCurTr.meanGood,...
        'Neu Tr SimVP Run - good trials','Neu Tr SimVP Rew - good trials');
   
    % compare single neuron correlation between aligned to run and aligned
    % to cue for good trials
    plotCompCorr(meanSimVPRun_LasttoCurTr.meanGood,meanSimVPCue_LasttoCurTr.meanGood,...
        'Neu Tr SimVP Run - good trials','Neu Tr SimVP Cue - good trials');
    
%     %% only consider nonzero trials
%     % compare single neuron correlation between good and bad trials --
%     % aligned to run
%     indNeu = meanSimVPRun_LasttoCurTr.nGoodNonZeroTr > 10 & meanSimVPRun_LasttoCurTr.nBadNonZeroTr > 10;
%     plotCompCorr(meanSimVPRun_LasttoCurTr.meanGoodNZ(indNeu),meanSimVPRun_LasttoCurTr.meanBadNZ(indNeu),...
%         'Neu None-zero Tr SimVP Run - good trials','Neu None-zero Tr SimVP Run - bad trials');
%     
%     % compare single neuron correlation between good and bad trials --
%     % aligned to reward
%     indNeu = meanSimVPRew_LasttoCurTr.nGoodNonZeroTr > 10 & meanSimVPRew_LasttoCurTr.nBadNonZeroTr > 10;
%     plotCompCorr(meanSimVPRew_LasttoCurTr.meanGoodNZ(indNeu),meanSimVPRew_LasttoCurTr.meanBadNZ(indNeu),...
%         'Neu None-zero Tr SimVP Rew - good trials','Neu None-zero Tr SimVP Rew - bad trials');
%         
%     % compare single neuron correlation between good and bad trials --
%     % aligned to cue
%     indNeu = meanSimVPCue_LasttoCurTr.nGoodNonZeroTr > 10 & meanSimVPCue_LasttoCurTr.nBadNonZeroTr > 10;
%     plotCompCorr(meanSimVPCue_LasttoCurTr.meanGoodNZ(indNeu),meanSimVPCue_LasttoCurTr.meanBadNZ(indNeu),...
%         'Neu None-zero Tr SimVP Cue - good trials','Neu None-zero Tr SimVP Cue - bad trials');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward
%     indNeu = meanSimVPRun_LasttoCurTr.nNonZeroTr > 10 & meanSimVPRew_LasttoCurTr.nNonZeroTr > 10;
%     plotCompCorr(meanSimVPRun_LasttoCurTr.meanNZ(indNeu),meanSimVPRew_LasttoCurTr.meanNZ(indNeu),...
%         'Neu None-zero Tr SimVP Run ','Neu None-zero Tr SimVP Rew');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue
%     indNeu = meanSimVPRun_LasttoCurTr.nNonZeroTr > 10 & meanSimVPCue_LasttoCurTr.nNonZeroTr > 10;
%     plotCompCorr(meanSimVPRun_LasttoCurTr.meanNZ(indNeu),meanSimVPCue_LasttoCurTr.meanNZ(indNeu),...
%         'Neu None-zero Tr SimVP Run','Neu None-zero Tr SimVP Cue');
%     
    % compare single neuron correlation between aligned to run and aligned
    % to reward for good trials
    indNeu = meanSimVPRun_LasttoCurTr.nGoodNonZeroTr > 5 & meanSimVPRew_LasttoCurTr.nGoodNonZeroTr > 5;
    plotCompCorr(meanSimVPRun_LasttoCurTr.meanGoodNZ(indNeu),meanSimVPRew_LasttoCurTr.meanGoodNZ(indNeu),...
        'Neu None-zero Tr SimVP Run - good trials','Neu None-zero Tr SimVP Rew - good trials');
   
    % compare single neuron correlation between aligned to run and aligned
    % to cue for good trials
    indNeu = meanSimVPRun_LasttoCurTr.nGoodNonZeroTr > 5 & meanSimVPCue_LasttoCurTr.nGoodNonZeroTr > 5;
    plotCompCorr(meanSimVPRun_LasttoCurTr.meanGoodNZ(indNeu),meanSimVPCue_LasttoCurTr.meanGoodNZ(indNeu),...
        'Neu None-zero Tr SimVP Run - good trials','Neu None-zero Tr SimVP Cue - good trials');
    
    save([path fileNameCorr],'meanSimVPRun_LasttoCurTr','meanSimVPRew_LasttoCurTr','meanSimVPCue_LasttoCurTr');
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
