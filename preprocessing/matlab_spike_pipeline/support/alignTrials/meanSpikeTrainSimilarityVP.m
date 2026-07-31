function meanSpikeTrainSimilarityVP(path,fileName,onlyRun,mazeSess,cost,intervalT)
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
    load(fullPath,'spikeTrainSimVPRun','spikeTrainSimVPRew','spikeTrainSimVPCue',...
        'nonZeroTrRun','nonZeroTrRew','nonZeroTrCue');
    
    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);    
    
    fullPath = [path fileName '_PeakFRAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];    
    if(exist(fullPath) == 0)
        disp('The _PeakFRAligned file does not exist');
        return;
    end
    load(fullPath,'trialNoNonStimGood','trialNoNonStimBad','trialNoStim','trialNoStimCtrl','pulseMeth');
    
    fileNameCorr = [fileName '_meanSpikeTrainSimVP_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_q' ...
            s '_intT' num2str(intervalT) '.mat'];
    
    neuronNo = length(spikeTrainSimVPRun);
    trialNo = size(spikeTrainSimVPRun{1},1);
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
        simArr = triu(spikeTrainSimVPRun{n},1);
        simArr = abs(simArr(:));
        meanSimVPRun.mean(n) = sum(simArr(isnan(simArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrRun{n} == 1);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPRun.meanNZ(n) = sum(simArr(isnan(simArr) == 0))/nElemNonZeroRun;
        meanSimVPRun.nNonZeroTr(n) = nNonZeroTr;
                
        simArr = triu(spikeTrainSimVPRew{n},1);
        simArr = abs(simArr(:));
        meanSimVPRew.mean(n) = sum(simArr(isnan(simArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrRew{n} == 1);
        nElemNonZeroRew = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPRew.meanNZ(n) = sum(simArr(isnan(simArr) == 0))/nElemNonZeroRew;
        meanSimVPRew.nNonZeroTr(n) = nNonZeroTr;
        
        simArr = triu(spikeTrainSimVPCue{n},1);
        simArr = abs(simArr(:));
        meanSimVPCue.mean(n) = sum(simArr(isnan(simArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrCue{n} == 1);
        nElemNonZeroCue = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPCue.meanNZ(n) = sum(simArr(isnan(simArr) == 0))/nElemNonZeroCue;
        meanSimVPCue.nNonZeroTr(n) = nNonZeroTr;
        
        %% non zero elements in good trials
        goodTrSim = spikeTrainSimVPRun{n}(indGoodTr,indGoodTr);
        goodTrSim = triu(goodTrSim,1);
        goodTrSim = abs(goodTrSim(:));
        meanSimVPRun.meanGood(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGood;
        indGoodNonZeroTrRun = indBadBeh == 0 & nonZeroTrRun{n} == 1;
        nGoodNonZeroTrRun = sum(indGoodNonZeroTrRun);
        nElemGoodNonZeroTrRun = (nGoodNonZeroTrRun*nGoodNonZeroTrRun-nGoodNonZeroTrRun)/2;
        meanSimVPRun.meanGoodNZ(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGoodNonZeroTrRun;
        meanSimVPRun.nGoodNonZeroTr(n) = nGoodNonZeroTrRun;
        
        goodTrSim = spikeTrainSimVPRew{n}(indGoodTr,indGoodTr);
        goodTrSim = triu(goodTrSim,1);
        goodTrSim = abs(goodTrSim(:));
        meanSimVPRew.meanGood(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGood;
        indGoodNonZeroTrRew = indBadBeh == 0 & nonZeroTrRew{n} == 1;
        nGoodNonZeroTrRew = sum(indGoodNonZeroTrRew);
        nElemGoodNonZeroTrRew = (nGoodNonZeroTrRew*nGoodNonZeroTrRew-nGoodNonZeroTrRew)/2;
        meanSimVPRew.meanGoodNZ(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGoodNonZeroTrRew;
        meanSimVPRew.nGoodNonZeroTr(n) = nGoodNonZeroTrRew;
        
        goodTrSim = spikeTrainSimVPCue{n}(indGoodTr,indGoodTr);
        goodTrSim = triu(goodTrSim,1);
        goodTrSim = abs(goodTrSim(:));
        meanSimVPCue.meanGood(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGood;
        indGoodNonZeroTrCue = indBadBeh == 0 & nonZeroTrCue{n} == 1;
        nGoodNonZeroTrCue = sum(indGoodNonZeroTrCue);
        nElemGoodNonZeroTrCue = (nGoodNonZeroTrCue*nGoodNonZeroTrCue-nGoodNonZeroTrCue)/2;
        meanSimVPCue.meanGoodNZ(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGoodNonZeroTrCue;
        meanSimVPCue.nGoodNonZeroTr(n) = nGoodNonZeroTrCue;
        
        %% non zero elements in bad trials
        badTrSim = spikeTrainSimVPRun{n}(indBadTr,indBadTr);
        badTrSim = triu(badTrSim,1);
        badTrSim = abs(badTrSim(:));
        meanSimVPRun.meanBad(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBad;
        indBadNonZeroTrRun = indBadBeh == 1 & nonZeroTrRun{n} == 1;
        nBadNonZeroTrRun = sum(indBadNonZeroTrRun);
        nElemBadNonZeroTrRun = (nBadNonZeroTrRun*nBadNonZeroTrRun-nBadNonZeroTrRun)/2;
        meanSimVPRun.meanBadNZ(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBadNonZeroTrRun;
        meanSimVPRun.nBadNonZeroTr(n) = nBadNonZeroTrRun;
        
        badTrSim = spikeTrainSimVPRew{n}(indBadTr,indBadTr);
        badTrSim = triu(badTrSim,1);
        badTrSim = abs(badTrSim(:));
        meanSimVPRew.meanBad(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBad;
        indBadNonZeroTrRew = indBadBeh == 1 & nonZeroTrRew{n} == 1;
        nBadNonZeroTrRew = sum(indBadNonZeroTrRew);
        nElemBadNonZeroTrRew = (nBadNonZeroTrRew*nBadNonZeroTrRew-nBadNonZeroTrRew)/2;
        meanSimVPRew.meanBadNZ(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBadNonZeroTrRew;
        meanSimVPRew.nBadNonZeroTr(n) = nBadNonZeroTrRew;
        
        badTrSim = spikeTrainSimVPCue{n}(indBadTr,indBadTr);
        badTrSim = triu(badTrSim,1);
        badTrSim = abs(badTrSim(:));
        meanSimVPCue.meanBad(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBad;
        indBadNonZeroTrCue = indBadBeh == 1 & nonZeroTrCue{n} == 1;
        nBadNonZeroTrCue = sum(indBadNonZeroTrCue);
        nElemBadNonZeroTrCue = (nBadNonZeroTrCue*nBadNonZeroTrCue-nBadNonZeroTrCue)/2;
        meanSimVPCue.meanBadNZ(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBadNonZeroTrCue;
        meanSimVPCue.nBadNonZeroTr(n) = nBadNonZeroTrCue;
        
        %% added by YIngxue on 2/14/2021
        noStimGoodTrSim = triu(spikeTrainSimVPRun{n}(trialNoNonStimGood,trialNoNonStimGood),1);
        noStimGoodTrSim = noStimGoodTrSim(:);
        nNoStimGood = length(trialNoNonStimGood);
        nElemNoStimGood = (nNoStimGood*nNoStimGood-nNoStimGood)/2;
        meanSimVPRun.meanNoStimGood(n) = sum(noStimGoodTrSim(isnan(noStimGoodTrSim) == 0))/nElemNoStimGood;
        indNoStimGoodNonZeroTrRun = intersect(trialNoNonStimGood, find(nonZeroTrRun{n} == 1));
        nNoStimGoodNonZeroTrRun = length(indNoStimGoodNonZeroTrRun);
        nElemNoStimGoodNonZeroTrRun = (nNoStimGoodNonZeroTrRun*nNoStimGoodNonZeroTrRun-nNoStimGoodNonZeroTrRun)/2;
        meanSimVPRun.meanNoStimGoodNZ(n) = sum(noStimGoodTrSim(isnan(noStimGoodTrSim) == 0))/nElemNoStimGoodNonZeroTrRun;
        meanSimVPRun.nNoStimGoodNonZeroTr(n) = nNoStimGoodNonZeroTrRun;
        
        noStimBadTrSim = triu(spikeTrainSimVPRun{n}(trialNoNonStimBad,trialNoNonStimBad),1);
        noStimBadTrSim = noStimBadTrSim(:);
        nNoStimBad = length(trialNoNonStimBad);
        nElemNoStimBad = (nNoStimBad*nNoStimBad-nNoStimBad)/2;
        meanSimVPRun.meanNoStimBad(n) = sum(noStimBadTrSim(isnan(noStimBadTrSim) == 0))/nElemNoStimBad;
        indNoStimBadNonZeroTrRun = intersect(trialNoNonStimBad, find(nonZeroTrRun{n} == 1));
        nNoStimBadNonZeroTrRun = length(indNoStimBadNonZeroTrRun);
        nElemNoStimBadNonZeroTrRun = (nNoStimBadNonZeroTrRun*nNoStimBadNonZeroTrRun-nNoStimBadNonZeroTrRun)/2;
        meanSimVPRun.meanNoStimBadNZ(n) = sum(noStimBadTrSim(isnan(noStimBadTrSim) == 0))/nElemNoStimBadNonZeroTrRun;
        meanSimVPRun.nNoStimBadNonZeroTr(n) = nNoStimBadNonZeroTrRun;
        
        meanSimVPRun.meanStim = [];
        meanSimVPRun.meanStimNZ = [];
        meanSimVPRun.nStimNonZeroTr = [];
        meanSimVPRun.meanStimCtrl = [];
        meanSimVPRun.meanStimCtrlNZ = [];
        meanSimVPRun.nStimCtrlNonZeroTr = [];
        for i = 1:length(pulseMeth) 
            stimTrSim = triu(spikeTrainSimVPRun{n}(trialNoStim{i},trialNoStim{i}),1);
            stimTrSim = stimTrSim(:);
            nStim = length(trialNoStim{i});
            nElemStim = (nStim*nStim-nStim)/2;
            meanSimVPRun.meanStim{i}(n) = sum(stimTrSim(isnan(stimTrSim) == 0))/nElemStim;
            indStimNonZeroTrRun = intersect(trialNoStim{i}, find(nonZeroTrRun{n} == 1));
            nStimNonZeroTrRun = length(indStimNonZeroTrRun);
            nElemStimNonZeroTrRun = (nStimNonZeroTrRun*nStimNonZeroTrRun-nStimNonZeroTrRun)/2;
            meanSimVPRun.meanStimNZ{i}(n) = sum(stimTrSim(isnan(stimTrSim) == 0))/nElemStimNonZeroTrRun;
            meanSimVPRun.nStimNonZeroTr{i}(n) = nStimNonZeroTrRun;

            stimCtrlTrSim = triu(spikeTrainSimVPRun{n}(trialNoStimCtrl{i},trialNoStimCtrl{i}),1);
            stimCtrlTrSim = stimCtrlTrSim(:);
            nStimCtrl = length(trialNoStimCtrl{i});
            nElemStimCtrl = (nStimCtrl*nStimCtrl-nStimCtrl)/2;
            meanSimVPRun.meanStimCtrl{i}(n) = sum(stimCtrlTrSim(isnan(stimCtrlTrSim) == 0))/nElemStimCtrl;
            indStimCtrlNonZeroTrRun = intersect(trialNoStimCtrl{i}, find(nonZeroTrRun{n} == 1));
            nStimCtrlNonZeroTrRun = length(indStimCtrlNonZeroTrRun);
            nElemStimCtrlNonZeroTrRun = (nStimCtrlNonZeroTrRun*nStimCtrlNonZeroTrRun-nStimCtrlNonZeroTrRun)/2;
            meanSimVPRun.meanStimCtrlNZ{i}(n) = sum(stimCtrlTrSim(isnan(stimCtrlTrSim) == 0))/nElemStimCtrlNonZeroTrRun;
            meanSimVPRun.nStimCtrlNonZeroTr{i}(n) = nStimCtrlNonZeroTrRun;
        end
        %%
    end
    
    % compare single neuron correlation between good and bad trials --
    % aligned to run
    plotCompCorr(meanSimVPRun.meanGood,meanSimVPRun.meanBad,...
        'Neu Tr SimVP Run - good trials','Neu Tr SimVP Run - bad trials');
    
%     % compare single neuron correlation between good and bad trials --
%     % aligned to reward
%     plotCompCorr(meanSimVPRew.meanGood,meanSimVPRew.meanBad,...
%         'Neu Tr SimVP Rew - good trials','Neu Tr SimVP Rew - bad trials');
%         
%     % compare single neuron correlation between good and bad trials --
%     % aligned to cue
%     plotCompCorr(meanSimVPCue.meanGood,meanSimVPCue.meanBad,...
%         'Neu Tr SimVP Cue - good trials','Neu Tr SimVP Cue - bad trials');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward
%     plotCompCorr(meanSimVPRun.mean,meanSimVPRew.mean,...
%         'Neu Tr SimVP Run','Neu Tr SimVP Rew');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue
%     plotCompCorr(meanSimVPRun.mean,meanSimVPCue.mean,...
%         'Neu Tr SimVP Run','Neu Tr SimVP Cue');
    
    % compare single neuron correlation between aligned to run and aligned
    % to reward for good trials
    plotCompCorr(meanSimVPRun.meanGood,meanSimVPRew.meanGood,...
        'Neu Tr SimVP Run - good trials','Neu Tr SimVP Rew - good trials');
   
    % compare single neuron correlation between aligned to run and aligned
    % to cue for good trials
    plotCompCorr(meanSimVPRun.meanGood,meanSimVPCue.meanGood,...
        'Neu Tr SimVP Run - good trials','Neu Tr SimVP Cue - good trials');
    
%     %% only consider nonzero trials
%     % compare single neuron correlation between good and bad trials --
%     % aligned to run
%     indNeu = meanSimVPRun.nGoodNonZeroTr > 10 & meanSimVPRun.nBadNonZeroTr > 10;
%     plotCompCorr(meanSimVPRun.meanGoodNZ(indNeu),meanSimVPRun.meanBadNZ(indNeu),...
%         'Neu None-zero Tr SimVP Run - good trials','Neu None-zero Tr SimVP Run - bad trials');
%     
%     % compare single neuron correlation between good and bad trials --
%     % aligned to reward
%     indNeu = meanSimVPRew.nGoodNonZeroTr > 10 & meanSimVPRew.nBadNonZeroTr > 10;
%     plotCompCorr(meanSimVPRew.meanGoodNZ(indNeu),meanSimVPRew.meanBadNZ(indNeu),...
%         'Neu None-zero Tr SimVP Rew - good trials','Neu None-zero Tr SimVP Rew - bad trials');
%         
%     % compare single neuron correlation between good and bad trials --
%     % aligned to cue
%     indNeu = meanSimVPCue.nGoodNonZeroTr > 10 & meanSimVPCue.nBadNonZeroTr > 10;
%     plotCompCorr(meanSimVPCue.meanGoodNZ(indNeu),meanSimVPCue.meanBadNZ(indNeu),...
%         'Neu None-zero Tr SimVP Cue - good trials','Neu None-zero Tr SimVP Cue - bad trials');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward
%     indNeu = meanSimVPRun.nNonZeroTr > 10 & meanSimVPRew.nNonZeroTr > 10;
%     plotCompCorr(meanSimVPRun.meanNZ(indNeu),meanSimVPRew.meanNZ(indNeu),...
%         'Neu None-zero Tr SimVP Run ','Neu None-zero Tr SimVP Rew');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue
%     indNeu = meanSimVPRun.nNonZeroTr > 10 & meanSimVPCue.nNonZeroTr > 10;
%     plotCompCorr(meanSimVPRun.meanNZ(indNeu),meanSimVPCue.meanNZ(indNeu),...
%         'Neu None-zero Tr SimVP Run','Neu None-zero Tr SimVP Cue');
%     
    % compare single neuron correlation between aligned to run and aligned
    % to reward for good trials
    indNeu = meanSimVPRun.nGoodNonZeroTr > 5 & meanSimVPRew.nGoodNonZeroTr > 5;
    plotCompCorr(meanSimVPRun.meanGoodNZ(indNeu),meanSimVPRew.meanGoodNZ(indNeu),...
        'Neu None-zero Tr SimVP Run - good trials','Neu None-zero Tr SimVP Rew - good trials');
   
    % compare single neuron correlation between aligned to run and aligned
    % to cue for good trials
    indNeu = meanSimVPRun.nGoodNonZeroTr > 5 & meanSimVPCue.nGoodNonZeroTr > 5;
    plotCompCorr(meanSimVPRun.meanGoodNZ(indNeu),meanSimVPCue.meanGoodNZ(indNeu),...
        'Neu None-zero Tr SimVP Run - good trials','Neu None-zero Tr SimVP Cue - good trials');
    
    save([path fileNameCorr],'meanSimVPRun','meanSimVPRew','meanSimVPCue');
end

function plotCompCorr(x,y,xlab,ylab)
    figure
    plot(x,y,'ro');
    hold on;
    corrAll = [x y];
    if(~isempty(corrAll))
        maxCorr = max(corrAll);
        minCorr = min(corrAll);
        if(minCorr == maxCorr)
            minCorr = minCorr-0.1;
        end
    else
        maxCorr = 1;
        minCorr = 0;
    end
    plot([minCorr maxCorr],[minCorr maxCorr],'k:');
    xlim(gca,[minCorr maxCorr]);
    ylim(gca,[minCorr maxCorr]);
    xlabel(xlab);
    ylabel(ylab);
end
