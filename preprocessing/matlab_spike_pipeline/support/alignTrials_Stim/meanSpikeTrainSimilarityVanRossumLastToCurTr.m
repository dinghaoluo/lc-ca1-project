function meanSpikeTrainSimilarityVanRossumLastToCurTr(path,fileName,onlyRun,mazeSess,tc,intervalT,intervalTMin)
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
            num2str(tc) '_intT' num2str(intervalTMin) '-' num2str(intervalT) '.mat'];
    
    neuronNo = length(spikeTrainSimVRRun_LasttoCurTr);
    trialNo = size(spikeTrainSimVRRun_LasttoCurTr{1},1);
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
        simArr = triu(spikeTrainSimVRRun_LasttoCurTr{n},1);
        simArr = abs(simArr(:));
        meanSimVRRun_LasttoCurTr.mean(n) = sum(simArr(isnan(simArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrRun_LasttoCurTr{n} == 1);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVRRun_LasttoCurTr.meanNZ(n) = sum(simArr(isnan(simArr) == 0))/nElemNonZeroRun;
        meanSimVRRun_LasttoCurTr.nNonZeroTr(n) = nNonZeroTr;
                
        simArr = triu(spikeTrainSimVRRew_LasttoCurTr{n},1);
        simArr = abs(simArr(:));
        meanSimVRRew_LasttoCurTr.mean(n) = sum(simArr(isnan(simArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrRew_LasttoCurTr{n} == 1);
        nElemNonZeroRew = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVRRew_LasttoCurTr.meanNZ(n) = sum(simArr(isnan(simArr) == 0))/nElemNonZeroRew;
        meanSimVRRew_LasttoCurTr.nNonZeroTr(n) = nNonZeroTr;
        
        simArr = triu(spikeTrainSimVRCue_LasttoCurTr{n},1);
        simArr = abs(simArr(:));
        meanSimVRCue_LasttoCurTr.mean(n) = sum(simArr(isnan(simArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrCue_LasttoCurTr{n} == 1);
        nElemNonZeroCue = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVRCue_LasttoCurTr.meanNZ(n) = sum(simArr(isnan(simArr) == 0))/nElemNonZeroCue;
        meanSimVRCue_LasttoCurTr.nNonZeroTr(n) = nNonZeroTr;
        
        %% non zero elements in good trials
        goodTrSim = spikeTrainSimVRRun_LasttoCurTr{n}(indGoodTr,indGoodTr);
        goodTrSim = triu(goodTrSim,1);
        goodTrSim = abs(goodTrSim(:));
        meanSimVRRun_LasttoCurTr.meanGood(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGood;
        indGoodNonZeroTrRun = indBadBeh == 0 & nonZeroTrRun_LasttoCurTr{n} == 1;
        nGoodNonZeroTrRun = sum(indGoodNonZeroTrRun);
        nElemGoodNonZeroTrRun = (nGoodNonZeroTrRun*nGoodNonZeroTrRun-nGoodNonZeroTrRun)/2;
        meanSimVRRun_LasttoCurTr.meanGoodNZ(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGoodNonZeroTrRun;
        meanSimVRRun_LasttoCurTr.nGoodNonZeroTr(n) = nGoodNonZeroTrRun;
        
        goodTrSim = spikeTrainSimVRRew_LasttoCurTr{n}(indGoodTr,indGoodTr);
        goodTrSim = triu(goodTrSim,1);
        goodTrSim = abs(goodTrSim(:));
        meanSimVRRew_LasttoCurTr.meanGood(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGood;
        indGoodNonZeroTrRew = indBadBeh == 0 & nonZeroTrRew_LasttoCurTr{n} == 1;
        nGoodNonZeroTrRew = sum(indGoodNonZeroTrRew);
        nElemGoodNonZeroTrRew = (nGoodNonZeroTrRew*nGoodNonZeroTrRew-nGoodNonZeroTrRew)/2;
        meanSimVRRew_LasttoCurTr.meanGoodNZ(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGoodNonZeroTrRew;
        meanSimVRRew_LasttoCurTr.nGoodNonZeroTr(n) = nGoodNonZeroTrRew;
        
        goodTrSim = spikeTrainSimVRCue_LasttoCurTr{n}(indGoodTr,indGoodTr);
        goodTrSim = triu(goodTrSim,1);
        goodTrSim = abs(goodTrSim(:));
        meanSimVRCue_LasttoCurTr.meanGood(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGood;
        indGoodNonZeroTrCue = indBadBeh == 0 & nonZeroTrCue_LasttoCurTr{n} == 1;
        nGoodNonZeroTrCue = sum(indGoodNonZeroTrCue);
        nElemGoodNonZeroTrCue = (nGoodNonZeroTrCue*nGoodNonZeroTrCue-nGoodNonZeroTrCue)/2;
        meanSimVRCue_LasttoCurTr.meanGoodNZ(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGoodNonZeroTrCue;
        meanSimVRCue_LasttoCurTr.nGoodNonZeroTr(n) = nGoodNonZeroTrCue;
        
        %% non zero elements in bad trials
        badTrSim = spikeTrainSimVRRun_LasttoCurTr{n}(indBadTr,indBadTr);
        badTrSim = triu(badTrSim,1);
        badTrSim = abs(badTrSim(:));
        meanSimVRRun_LasttoCurTr.meanBad(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBad;
        indBadNonZeroTrRun = indBadBeh == 1 & nonZeroTrRun_LasttoCurTr{n} == 1;
        nBadNonZeroTrRun = sum(indBadNonZeroTrRun);
        nElemBadNonZeroTrRun = (nBadNonZeroTrRun*nBadNonZeroTrRun-nBadNonZeroTrRun)/2;
        meanSimVRRun_LasttoCurTr.meanBadNZ(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBadNonZeroTrRun;
        meanSimVRRun_LasttoCurTr.nBadNonZeroTr(n) = nBadNonZeroTrRun;
        
        badTrSim = spikeTrainSimVRRew_LasttoCurTr{n}(indBadTr,indBadTr);
        badTrSim = triu(badTrSim,1);
        badTrSim = abs(badTrSim(:));
        meanSimVRRew_LasttoCurTr.meanBad(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBad;
        indBadNonZeroTrRew = indBadBeh == 1 & nonZeroTrRew_LasttoCurTr{n} == 1;
        nBadNonZeroTrRew = sum(indBadNonZeroTrRew);
        nElemBadNonZeroTrRew = (nBadNonZeroTrRew*nBadNonZeroTrRew-nBadNonZeroTrRew)/2;
        meanSimVRRew_LasttoCurTr.meanBadNZ(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBadNonZeroTrRew;
        meanSimVRRew_LasttoCurTr.nBadNonZeroTr(n) = nBadNonZeroTrRew;
        
        badTrSim = spikeTrainSimVRCue_LasttoCurTr{n}(indBadTr,indBadTr);
        badTrSim = triu(badTrSim,1);
        badTrSim = abs(badTrSim(:));
        meanSimVRCue_LasttoCurTr.meanBad(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBad;
        indBadNonZeroTrCue = indBadBeh == 1 & nonZeroTrCue_LasttoCurTr{n} == 1;
        nBadNonZeroTrCue = sum(indBadNonZeroTrCue);
        nElemBadNonZeroTrCue = (nBadNonZeroTrCue*nBadNonZeroTrCue-nBadNonZeroTrCue)/2;
        meanSimVRCue_LasttoCurTr.meanBadNZ(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBadNonZeroTrCue;
        meanSimVRCue_LasttoCurTr.nBadNonZeroTr(n) = nBadNonZeroTrCue;
    end
    
    % compare single neuron correlation between good and bad trials --
    % aligned to run
    plotCompCorr(meanSimVRRun_LasttoCurTr.meanGood,meanSimVRRun_LasttoCurTr.meanBad,...
        'Neu Tr SimVR Run - good trials','Neu Tr SimVR Run - bad trials');
    
%     % compare single neuron correlation between good and bad trials --
%     % aligned to reward
%     plotCompCorr(meanSimVRRew_LasttoCurTr.meanGood,meanSimVRRew_LasttoCurTr.meanBad,...
%         'Neu Tr SimVR Rew - good trials','Neu Tr SimVR Rew - bad trials');
%         
%     % compare single neuron correlation between good and bad trials --
%     % aligned to cue
%     plotCompCorr(meanSimVRCue_LasttoCurTr.meanGood,meanSimVRCue_LasttoCurTr.meanBad,...
%         'Neu Tr SimVR Cue - good trials','Neu Tr SimVR Cue - bad trials');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward
%     plotCompCorr(meanSimVRRun_LasttoCurTr.mean,meanSimVRRew_LasttoCurTr.mean,...
%         'Neu Tr SimVR Run','Neu Tr SimVR Rew');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue
%     plotCompCorr(meanSimVRRun_LasttoCurTr.mean,meanSimVRCue_LasttoCurTr.mean,...
%         'Neu Tr SimVR Run','Neu Tr SimVR Cue');
%     
    % compare single neuron correlation between aligned to run and aligned
    % to reward for good trials
    plotCompCorr(meanSimVRRun_LasttoCurTr.meanGood,meanSimVRRew_LasttoCurTr.meanGood,...
        'Neu Tr SimVR Run - good trials','Neu Tr SimVR Rew - good trials');
   
    % compare single neuron correlation between aligned to run and aligned
    % to cue for good trials
    plotCompCorr(meanSimVRRun_LasttoCurTr.meanGood,meanSimVRCue_LasttoCurTr.meanGood,...
        'Neu Tr SimVR Run - good trials','Neu Tr SimVR Cue - good trials');
    
    %% only consider nonzero trials
%     % compare single neuron correlation between good and bad trials --
%     % aligned to run
%     indNeu = meanSimVRRun_LasttoCurTr.nGoodNonZeroTr > 10 & meanSimVRRun_LasttoCurTr.nBadNonZeroTr > 10;
%     plotCompCorr(meanSimVRRun_LasttoCurTr.meanGoodNZ(indNeu),meanSimVRRun_LasttoCurTr.meanBadNZ(indNeu),...
%         'Neu None-zero Tr SimVR Run - good trials','Neu None-zero Tr SimVR Run - bad trials');
%     
%     % compare single neuron correlation between good and bad trials --
%     % aligned to reward
%     indNeu = meanSimVRRew_LasttoCurTr.nGoodNonZeroTr > 10 & meanSimVRRew_LasttoCurTr.nBadNonZeroTr > 10;
%     plotCompCorr(meanSimVRRew_LasttoCurTr.meanGoodNZ(indNeu),meanSimVRRew_LasttoCurTr.meanBadNZ(indNeu),...
%         'Neu None-zero Tr SimVR Rew - good trials','Neu None-zero Tr SimVR Rew - bad trials');
%         
%     % compare single neuron correlation between good and bad trials --
%     % aligned to cue
%     indNeu = meanSimVRCue_LasttoCurTr.nGoodNonZeroTr > 10 & meanSimVRCue_LasttoCurTr.nBadNonZeroTr > 10;
%     plotCompCorr(meanSimVRCue_LasttoCurTr.meanGoodNZ(indNeu),meanSimVRCue_LasttoCurTr.meanBadNZ(indNeu),...
%         'Neu None-zero Tr SimVR Cue - good trials','Neu None-zero Tr SimVR Cue - bad trials');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward
%     indNeu = meanSimVRRun_LasttoCurTr.nNonZeroTr > 10 & meanSimVRRew_LasttoCurTr.nNonZeroTr > 10;
%     plotCompCorr(meanSimVRRun_LasttoCurTr.meanNZ(indNeu),meanSimVRRew_LasttoCurTr.meanNZ(indNeu),...
%         'Neu None-zero Tr SimVR Run ','Neu None-zero Tr SimVR Rew');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue
%     indNeu = meanSimVRRun_LasttoCurTr.nNonZeroTr > 10 & meanSimVRCue_LasttoCurTr.nNonZeroTr > 10;
%     plotCompCorr(meanSimVRRun_LasttoCurTr.meanNZ(indNeu),meanSimVRCue_LasttoCurTr.meanNZ(indNeu),...
%         'Neu None-zero Tr SimVR Run','Neu None-zero Tr SimVR Cue');
%     
    % compare single neuron correlation between aligned to run and aligned
    % to reward for good trials
    indNeu = meanSimVRRun_LasttoCurTr.nGoodNonZeroTr > 5 & meanSimVRRew_LasttoCurTr.nGoodNonZeroTr > 5;
    plotCompCorr(meanSimVRRun_LasttoCurTr.meanGoodNZ(indNeu),meanSimVRRew_LasttoCurTr.meanGoodNZ(indNeu),...
        'Neu None-zero Tr SimVR Run - good trials','Neu None-zero Tr SimVR Rew - good trials');
   
    % compare single neuron correlation between aligned to run and aligned
    % to cue for good trials
    indNeu = meanSimVRRun_LasttoCurTr.nGoodNonZeroTr > 5 & meanSimVRCue_LasttoCurTr.nGoodNonZeroTr > 5;
    plotCompCorr(meanSimVRRun_LasttoCurTr.meanGoodNZ(indNeu),meanSimVRCue_LasttoCurTr.meanGoodNZ(indNeu),...
        'Neu None-zero Tr SimVR Run - good trials','Neu None-zero Tr SimVR Cue - good trials');
    
    save([path fileNameCorr],'meanSimVRRun_LasttoCurTr','meanSimVRRew_LasttoCurTr','meanSimVRCue_LasttoCurTr');
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
