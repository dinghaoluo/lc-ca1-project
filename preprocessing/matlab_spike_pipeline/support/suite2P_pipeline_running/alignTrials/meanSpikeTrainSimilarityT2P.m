function meanSpikeTrainSimilarityT2P(path,fileName,onlyRun,mazeSess,intervalT)
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
    
    fullPath = [path fileName '_PeakFRAlignedRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];    
    if(exist(fullPath) == 0)
        disp('The _PeakFRAligned file does not exist');
        return;
    end
    load(fullPath,'trialNoNonStimGood','trialNoNonStimBad','trialNoStim','trialNoStimCtrl','pulseMeth');
    
    fileNameCorr = [fileName '_meanSpikeTrainSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    
    GlobalConst2P;
        
    neuronNo = length(spikeTrainSimTRun);
    trialNo = size(spikeTrainSimTRun{1},1);
    indTr = 1:trialNo;
    indSelTr = indTr > startTrNo;
    trialNo = sum(indSelTr);
    nElem = (trialNo*trialNo-trialNo)/2;
    
    indBadBeh = behPar.indTrBadBeh;
    indGoodTr = indBadBeh == 0 & indSelTr;
    indBadTr = indBadBeh == 1 & indSelTr;
    nGoodTr = sum(indGoodTr);
    nBadTr = sum(indBadTr);
    nElemGood = (nGoodTr*nGoodTr-nGoodTr)/2;
    nElemBad = (nBadTr*nBadTr-nBadTr)/2;
    
    for n = 1:neuronNo
        simArr = triu(spikeTrainSimTRun{n}(indSelTr,indSelTr),1);
        simArr = simArr(:);
        meanSimTRun.mean(n) = sum(simArr(isnan(simArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrRun{n} == 1 & indSelTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimTRun.meanNZ(n) = sum(simArr(isnan(simArr) == 0))/nElemNonZeroRun;
        meanSimTRun.nNonZeroTr(n) = nNonZeroTr;
                
        simArr = triu(spikeTrainSimTRew{n}(indSelTr,indSelTr),1);
        simArr = simArr(:);
        meanSimTRew.mean(n) = sum(simArr(isnan(simArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrRew{n} == 1 & indSelTr);
        nElemNonZeroRew = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimTRew.meanNZ(n) = sum(simArr(isnan(simArr) == 0))/nElemNonZeroRew;
        meanSimTRew.nNonZeroTr(n) = nNonZeroTr;
        
        simArr = triu(spikeTrainSimTCue{n}(indSelTr,indSelTr),1);
        simArr = simArr(:);
        meanSimTCue.mean(n) = sum(simArr(isnan(simArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrCue{n} == 1 & indSelTr);
        nElemNonZeroCue = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimTCue.meanNZ(n) = sum(simArr(isnan(simArr) == 0))/nElemNonZeroCue;
        meanSimTCue.nNonZeroTr(n) = nNonZeroTr;
        
        goodTrSim = spikeTrainSimTRun{n}(indGoodTr,indGoodTr);
        goodTrSim = triu(goodTrSim,1);
        goodTrSim = goodTrSim(:);
        meanSimTRun.meanGood(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGood;
        indGoodNonZeroTrRun = indGoodTr == 1 & nonZeroTrRun{n} == 1;
        nGoodNonZeroTrRun = sum(indGoodNonZeroTrRun);
        nElemGoodNonZeroTrRun = (nGoodNonZeroTrRun*nGoodNonZeroTrRun-nGoodNonZeroTrRun)/2;
        meanSimTRun.meanGoodNZ(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGoodNonZeroTrRun;
        meanSimTRun.nGoodNonZeroTr(n) = nGoodNonZeroTrRun;
        
        goodTrSim = spikeTrainSimTRew{n}(indGoodTr,indGoodTr);
        goodTrSim = triu(goodTrSim,1);
        goodTrSim = goodTrSim(:);
        meanSimTRew.meanGood(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGood;
        indGoodNonZeroTrRew = indGoodTr == 1 & nonZeroTrRew{n} == 1;
        nGoodNonZeroTrRew = sum(indGoodNonZeroTrRew);
        nElemGoodNonZeroTrRew = (nGoodNonZeroTrRew*nGoodNonZeroTrRew-nGoodNonZeroTrRew)/2;
        meanSimTRew.meanGoodNZ(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGoodNonZeroTrRew;
        meanSimTRew.nGoodNonZeroTr(n) = nGoodNonZeroTrRew;
        
        goodTrSim = spikeTrainSimTCue{n}(indGoodTr,indGoodTr);
        goodTrSim = triu(goodTrSim,1);
        goodTrSim = goodTrSim(:);
        meanSimTCue.meanGood(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGood;
        indGoodNonZeroTrCue = indGoodTr == 1 & nonZeroTrCue{n} == 1;
        nGoodNonZeroTrCue = sum(indGoodNonZeroTrCue);
        nElemGoodNonZeroTrCue = (nGoodNonZeroTrCue*nGoodNonZeroTrCue-nGoodNonZeroTrCue)/2;
        meanSimTCue.meanGoodNZ(n) = sum(goodTrSim(isnan(goodTrSim) == 0))/nElemGoodNonZeroTrCue;
        meanSimTCue.nGoodNonZeroTr(n) = nGoodNonZeroTrCue;
        
        badTrSim = spikeTrainSimTRun{n}(indBadTr,indBadTr);
        badTrSim = triu(badTrSim,1);
        badTrSim = badTrSim(:);
        meanSimTRun.meanBad(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBad;
        indBadNonZeroTrRun = indBadTr == 1 & nonZeroTrRun{n} == 1;
        nBadNonZeroTrRun = sum(indBadNonZeroTrRun);
        nElemBadNonZeroTrRun = (nBadNonZeroTrRun*nBadNonZeroTrRun-nBadNonZeroTrRun)/2;
        meanSimTRun.meanBadNZ(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBadNonZeroTrRun;
        meanSimTRun.nBadNonZeroTr(n) = nBadNonZeroTrRun;
        
        badTrSim = spikeTrainSimTRew{n}(indBadTr,indBadTr);
        badTrSim = triu(badTrSim,1);
        badTrSim = badTrSim(:);
        meanSimTRew.meanBad(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBad;
        indBadNonZeroTrRew = indBadTr == 1 & nonZeroTrRew{n} == 1;
        nBadNonZeroTrRew = sum(indBadNonZeroTrRew);
        nElemBadNonZeroTrRew = (nBadNonZeroTrRew*nBadNonZeroTrRew-nBadNonZeroTrRew)/2;
        meanSimTRew.meanBadNZ(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBadNonZeroTrRew;
        meanSimTRew.nBadNonZeroTr(n) = nBadNonZeroTrRew;
        
        badTrSim = spikeTrainSimTCue{n}(indBadTr,indBadTr);
        badTrSim = triu(badTrSim,1);
        badTrSim = badTrSim(:);
        meanSimTCue.meanBad(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBad;
        indBadNonZeroTrCue = indBadTr == 1 & nonZeroTrCue{n} == 1;
        nBadNonZeroTrCue = sum(indBadNonZeroTrCue);
        nElemBadNonZeroTrCue = (nBadNonZeroTrCue*nBadNonZeroTrCue-nBadNonZeroTrCue)/2;
        meanSimTCue.meanBadNZ(n) = sum(badTrSim(isnan(badTrSim) == 0))/nElemBadNonZeroTrCue;
        meanSimTCue.nBadNonZeroTr(n) = nBadNonZeroTrCue;
        
        %% added by YIngxue on 2/14/2021
        trialNoNonStimGood = setdiff(trialNoNonStimGood,1:startTrNo);
        noStimGoodTrSim = triu(spikeTrainSimTRun{n}(trialNoNonStimGood,trialNoNonStimGood),1);
        noStimGoodTrSim = noStimGoodTrSim(:);
        nNoStimGood = length(trialNoNonStimGood);
        nElemNoStimGood = (nNoStimGood*nNoStimGood-nNoStimGood)/2;
        meanSimTRun.meanNoStimGood(n) = sum(noStimGoodTrSim(isnan(noStimGoodTrSim) == 0))/nElemNoStimGood;
        indNoStimGoodNonZeroTrRun = intersect(trialNoNonStimGood, find(nonZeroTrRun{n} == 1));
        nNoStimGoodNonZeroTrRun = length(indNoStimGoodNonZeroTrRun);
        nElemNoStimGoodNonZeroTrRun = (nNoStimGoodNonZeroTrRun*nNoStimGoodNonZeroTrRun-nNoStimGoodNonZeroTrRun)/2;
        meanSimTRun.meanNoStimGoodNZ(n) = sum(noStimGoodTrSim(isnan(noStimGoodTrSim) == 0))/nElemNoStimGoodNonZeroTrRun;
        meanSimTRun.nNoStimGoodNonZeroTr(n) = nNoStimGoodNonZeroTrRun;
        
        trialNoNonStimBad = setdiff(trialNoNonStimBad,1:startTrNo);
        noStimBadTrSim = triu(spikeTrainSimTRun{n}(trialNoNonStimBad,trialNoNonStimBad),1);
        noStimBadTrSim = noStimBadTrSim(:);
        nNoStimBad = length(trialNoNonStimBad);
        nElemNoStimBad = (nNoStimBad*nNoStimBad-nNoStimBad)/2;
        meanSimTRun.meanNoStimBad(n) = sum(noStimBadTrSim(isnan(noStimBadTrSim) == 0))/nElemNoStimBad;
        indNoStimBadNonZeroTrRun = intersect(trialNoNonStimBad, find(nonZeroTrRun{n} == 1));
        nNoStimBadNonZeroTrRun = length(indNoStimBadNonZeroTrRun);
        nElemNoStimBadNonZeroTrRun = (nNoStimBadNonZeroTrRun*nNoStimBadNonZeroTrRun-nNoStimBadNonZeroTrRun)/2;
        meanSimTRun.meanNoStimBadNZ(n) = sum(noStimBadTrSim(isnan(noStimBadTrSim) == 0))/nElemNoStimBadNonZeroTrRun;
        meanSimTRun.nNoStimBadNonZeroTr(n) = nNoStimBadNonZeroTrRun;
        
        meanSimTRun.meanStim = [];
        meanSimTRun.meanStimNZ = [];
        meanSimTRun.nStimNonZeroTr = [];
        meanSimTRun.meanStimCtrl = [];
        meanSimTRun.meanStimCtrlNZ = [];
        meanSimTRun.nStimCtrlNonZeroTr = [];
        for i = 1:length(pulseMeth) 
            trialNoStim{i} = setdiff(trialNoStim{i},1:startTrNo);
            stimTrSim = triu(spikeTrainSimTRun{n}(trialNoStim{i},trialNoStim{i}),1);
            stimTrSim = stimTrSim(:);
            nStim = length(trialNoStim{i});
            nElemStim = (nStim*nStim-nStim)/2;
            meanSimTRun.meanStim{i}(n) = sum(stimTrSim(isnan(stimTrSim) == 0))/nElemStim;
            indStimNonZeroTrRun = intersect(trialNoStim{i}, find(nonZeroTrRun{n} == 1));
            nStimNonZeroTrRun = length(indStimNonZeroTrRun);
            nElemStimNonZeroTrRun = (nStimNonZeroTrRun*nStimNonZeroTrRun-nStimNonZeroTrRun)/2;
            meanSimTRun.meanStimNZ{i}(n) = sum(stimTrSim(isnan(stimTrSim) == 0))/nElemStimNonZeroTrRun;
            meanSimTRun.nStimNonZeroTr{i}(n) = nStimNonZeroTrRun;

            trialNoStimCtrl{i} = setdiff(trialNoStimCtrl{i},1:startTrNo);
            stimCtrlTrSim = triu(spikeTrainSimTRun{n}(trialNoStimCtrl{i},trialNoStimCtrl{i}),1);
            stimCtrlTrSim = stimCtrlTrSim(:);
            nStimCtrl = length(trialNoStimCtrl{i});
            nElemStimCtrl = (nStimCtrl*nStimCtrl-nStimCtrl)/2;
            meanSimTRun.meanStimCtrl{i}(n) = sum(stimCtrlTrSim(isnan(stimCtrlTrSim) == 0))/nElemStimCtrl;
            indStimCtrlNonZeroTrRun = intersect(trialNoStimCtrl{i}, find(nonZeroTrRun{n} == 1));
            nStimCtrlNonZeroTrRun = length(indStimCtrlNonZeroTrRun);
            nElemStimCtrlNonZeroTrRun = (nStimCtrlNonZeroTrRun*nStimCtrlNonZeroTrRun-nStimCtrlNonZeroTrRun)/2;
            meanSimTRun.meanStimCtrlNZ{i}(n) = sum(stimCtrlTrSim(isnan(stimCtrlTrSim) == 0))/nElemStimCtrlNonZeroTrRun;
            meanSimTRun.nStimCtrlNonZeroTr{i}(n) = nStimCtrlNonZeroTrRun;
        end
        %%
    end
    
%     % compare single neuron correlation between good and bad trials --
%     % aligned to run
    plotCompCorr(meanSimTRun.meanGood,meanSimTRun.meanBad,...
        'Neu Tr SimT Run - good trials','Neu Tr SimT Run - bad trials');
%     
%     % compare single neuron correlation between good and bad trials --
%     % aligned to reward
%     plotCompCorr(meanSimTRew.meanGood,meanSimTRew.meanBad,...
%         'Neu Tr SimT Rew - good trials','Neu Tr SimT Rew - bad trials');
%         
%     % compare single neuron correlation between good and bad trials --
%     % aligned to cue
%     plotCompCorr(meanSimTCue.meanGood,meanSimTCue.meanBad,...
%         'Neu Tr SimT Cue - good trials','Neu Tr SimT Cue - bad trials');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward
%     plotCompCorr(meanSimTRun.mean,meanSimTRew.mean,...
%         'Neu Tr SimT Run','Neu Tr SimT Rew');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue
%     plotCompCorr(meanSimTRun.mean,meanSimTCue.mean,...
%         'Neu Tr SimT Run','Neu Tr SimT Cue');
%     
    % compare single neuron correlation between aligned to run and aligned
    % to reward for good trials
    plotCompCorr(meanSimTRun.meanGood,meanSimTRew.meanGood,...
        'Neu Tr SimT Run - good trials','Neu Tr SimT Rew - good trials');
   
    % compare single neuron correlation between aligned to run and aligned
    % to cue for good trials
    plotCompCorr(meanSimTRun.meanGood,meanSimTCue.meanGood,...
        'Neu Tr SimT Run - good trials','Neu Tr SimT Cue - good trials');
    
    %% only consider nonzero trials
    % compare single neuron correlation between good and bad trials --
    % aligned to run
%     indNeu = meanSimTRun.nGoodNonZeroTr > 10 & meanSimTRun.nBadNonZeroTr > 5;
%     plotCompCorr(meanSimTRun.meanGoodNZ(indNeu),meanSimTRun.meanBadNZ(indNeu),...
%         'Neu None-zero Tr SimT Run - good trials','Neu None-zero Tr SimT Run - bad trials');
    
%     % compare single neuron correlation between good and bad trials --
%     % aligned to reward
%     indNeu = meanSimTRew.nGoodNonZeroTr > 10 & meanSimTRew.nBadNonZeroTr > 10;
%     plotCompCorr(meanSimTRew.meanGoodNZ(indNeu),meanSimTRew.meanBadNZ(indNeu),...
%         'Neu None-zero Tr SimT Rew - good trials','Neu None-zero Tr SimT Rew - bad trials');
%         
%     % compare single neuron correlation between good and bad trials --
%     % aligned to cue
%     indNeu = meanSimTCue.nGoodNonZeroTr > 10 & meanSimTCue.nBadNonZeroTr > 10;
%     plotCompCorr(meanSimTCue.meanGoodNZ(indNeu),meanSimTCue.meanBadNZ(indNeu),...
%         'Neu None-zero Tr SimT Cue - good trials','Neu None-zero Tr SimT Cue - bad trials');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward
%     indNeu = meanSimTRun.nNonZeroTr > 10 & meanSimTRew.nNonZeroTr > 10;
%     plotCompCorr(meanSimTRun.meanNZ(indNeu),meanSimTRew.meanNZ(indNeu),...
%         'Neu None-zero Tr SimT Run ','Neu None-zero Tr SimT Rew');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue
%     indNeu = meanSimTRun.nNonZeroTr > 10 & meanSimTCue.nNonZeroTr > 10;
%     plotCompCorr(meanSimTRun.meanNZ(indNeu),meanSimTCue.meanNZ(indNeu),...
%         'Neu None-zero Tr SimT Run','Neu None-zero Tr SimT Cue');
    
    % compare single neuron correlation between aligned to run and aligned
    % to reward for good trials
    indNeu = meanSimTRun.nGoodNonZeroTr > 5 & meanSimTRew.nGoodNonZeroTr > 5;
    plotCompCorr(meanSimTRun.meanGoodNZ(indNeu),meanSimTRew.meanGoodNZ(indNeu),...
        'Neu None-zero Tr SimT Run - good trials','Neu None-zero Tr SimT Rew - good trials');
   
    % compare single neuron correlation between aligned to run and aligned
    % to cue for good trials
    indNeu = meanSimTRun.nGoodNonZeroTr > 5 & meanSimTCue.nGoodNonZeroTr > 5;
    plotCompCorr(meanSimTRun.meanGoodNZ(indNeu),meanSimTCue.meanGoodNZ(indNeu),...
        'Neu None-zero Tr SimT Run - good trials','Neu None-zero Tr SimT Cue - good trials');
    
    save([path fileNameCorr],'meanSimTRun','meanSimTRew','meanSimTCue');
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
