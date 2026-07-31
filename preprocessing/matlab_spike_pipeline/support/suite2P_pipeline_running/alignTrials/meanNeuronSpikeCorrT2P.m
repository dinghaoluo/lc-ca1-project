function meanNeuronSpikeCorrT2P(path,fileName,onlyRun,mazeSess,intervalT)
% single neuron level mean spike correlation across trials
% meanNeuronSpikeCorrT('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,10)

    if(nargin == 4)
        intervalT = 0;
    end
    fullPath = [path fileName '_spikesCorrTAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    fileNameCorr = [fileName '_meanSpikesCorrTAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
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
    
    fullPath = [path fileName '_PeakFRAlignedRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];    
    if(exist(fullPath) == 0)
        disp('The _PeakFRAligned file does not exist');
        return;
    end
    load(fullPath,'trialNoNonStimGood','trialNoNonStimBad','trialNoStim','trialNoStimCtrl','pulseMeth');
    
    GlobalConst2P;
    
    neuronNo = length(spikeCorrTRun);
    trialNo = size(spikeCorrTRun{1},1);
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
        corrArr = triu(spikeCorrTRun{n}(indSelTr,indSelTr),1);
        corrArr = corrArr(:);        
        meanCorrTRun.mean(n) = sum(corrArr(isnan(corrArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrRun{n} == 1 & indSelTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanCorrTRun.meanNZ(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanCorrTRun.nNonZeroTr(n) = nNonZeroTr;
        
        corrArr = triu(spikeCorrTRew{n}(indSelTr,indSelTr),1);
        corrArr = corrArr(:);  
        meanCorrTRew.mean(n) = sum(corrArr(isnan(corrArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrRew{n} == 1 & indSelTr);
        nElemNonZeroRew = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanCorrTRew.meanNZ(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRew;
        meanCorrTRew.nNonZeroTr(n) = nNonZeroTr;
        
        corrArr = triu(spikeCorrTCue{n}(indSelTr,indSelTr),1);
        corrArr = corrArr(:);
        meanCorrTCue.mean(n) = sum(corrArr(isnan(corrArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTrCue{n} == 1 & indSelTr);
        nElemNonZeroCue = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanCorrTCue.meanNZ(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroCue;
        meanCorrTCue.nNonZeroTr(n) = nNonZeroTr;
        
        goodTrCorr = spikeCorrTRun{n}(indGoodTr,indGoodTr);
        goodTrCorr = triu(goodTrCorr,1);
        goodTrCorr = goodTrCorr(:);
        meanCorrTRun.meanGood(n) = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGood;
        indGoodNonZeroTrRun = indGoodTr == 1 & nonZeroTrRun{n} == 1;
        nGoodNonZeroTrRun = sum(indGoodNonZeroTrRun);
        nElemGoodNonZeroTrRun = (nGoodNonZeroTrRun*nGoodNonZeroTrRun-nGoodNonZeroTrRun)/2;
        meanCorrTRun.meanGoodNZ(n) = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGoodNonZeroTrRun;
        meanCorrTRun.nGoodNonZeroTr(n) = nGoodNonZeroTrRun;
        
        goodTrCorr = triu(spikeCorrTRew{n}(indGoodTr,indGoodTr),1);
        goodTrCorr = goodTrCorr(:);
        meanCorrTRew.meanGood(n) = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGood;
        indGoodNonZeroTrRew = indGoodTr == 1 & nonZeroTrRew{n} == 1;
        nGoodNonZeroTrRew = sum(indGoodNonZeroTrRew);
        nElemGoodNonZeroTrRew = (nGoodNonZeroTrRew*nGoodNonZeroTrRew-nGoodNonZeroTrRew)/2;
        meanCorrTRew.meanGoodNZ(n) = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGoodNonZeroTrRew;
        meanCorrTRew.nGoodNonZeroTr(n) = nGoodNonZeroTrRew;
        
        goodTrCorr = triu(spikeCorrTCue{n}(indGoodTr,indGoodTr),1);
        goodTrCorr = goodTrCorr(:);
        meanCorrTCue.meanGood(n) = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGood;
        indGoodNonZeroTrCue = indGoodTr == 1 & nonZeroTrCue{n} == 1;
        nGoodNonZeroTrCue = sum(indGoodNonZeroTrCue);
        nElemGoodNonZeroTrCue = (nGoodNonZeroTrCue*nGoodNonZeroTrCue-nGoodNonZeroTrCue)/2;
        meanCorrTCue.meanGoodNZ(n) = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGoodNonZeroTrCue;
        meanCorrTCue.nGoodNonZeroTr(n) = nGoodNonZeroTrCue;
        
        badTrCorr = triu(spikeCorrTRun{n}(indBadTr,indBadTr),1);
        badTrCorr = badTrCorr(:);
        meanCorrTRun.meanBad(n) = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBad;
        indBadNonZeroTrRun = indBadTr == 1 & nonZeroTrRun{n} == 1;
        nBadNonZeroTrRun = sum(indBadNonZeroTrRun);
        nElemBadNonZeroTrRun = (nBadNonZeroTrRun*nBadNonZeroTrRun-nBadNonZeroTrRun)/2;
        meanCorrTRun.meanBadNZ(n) = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBadNonZeroTrRun;
        meanCorrTRun.nBadNonZeroTr(n) = nBadNonZeroTrRun;
        
        badTrCorr = triu(spikeCorrTRew{n}(indBadTr,indBadTr),1);
        badTrCorr = badTrCorr(:);
        meanCorrTRew.meanBad(n) = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBad;
        indBadNonZeroTrRew = indBadTr == 1 & nonZeroTrRew{n} == 1;
        nBadNonZeroTrRew = sum(indBadNonZeroTrRew);
        nElemBadNonZeroTrRew = (nBadNonZeroTrRew*nBadNonZeroTrRew-nBadNonZeroTrRew)/2;
        meanCorrTRew.meanBadNZ(n) = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBadNonZeroTrRew;
        meanCorrTRew.nBadNonZeroTr(n) = nBadNonZeroTrRew;
        
        badTrCorr = triu(spikeCorrTCue{n}(indBadTr,indBadTr),1);
        badTrCorr = badTrCorr(:);
        meanCorrTCue.meanBad(n) = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBad;
        indBadNonZeroTrCue = indBadTr == 1 & nonZeroTrCue{n} == 1;
        nBadNonZeroTrCue = sum(indBadNonZeroTrCue);
        nElemBadNonZeroTrCue = (nBadNonZeroTrCue*nBadNonZeroTrCue-nBadNonZeroTrCue)/2;
        meanCorrTCue.meanBadNZ(n) = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBadNonZeroTrCue;
        meanCorrTCue.nBadNonZeroTr(n) = nBadNonZeroTrCue;
        
        %% added by YIngxue on 2/14/2021
        trialNoNonStimGood = setdiff(trialNoNonStimGood,1:startTrNo);
        noStimGoodTrCorr = triu(spikeCorrTRun{n}(trialNoNonStimGood,trialNoNonStimGood),1);
        noStimGoodTrCorr = noStimGoodTrCorr(:);
        nNoStimGood = length(trialNoNonStimGood);
        nElemNoStimGood = (nNoStimGood*nNoStimGood-nNoStimGood)/2;
        meanCorrTRun.meanNoStimGood(n) = sum(noStimGoodTrCorr(isnan(noStimGoodTrCorr) == 0))/nElemNoStimGood;
        indNoStimGoodNonZeroTrRun = intersect(trialNoNonStimGood, find(nonZeroTrRun{n} == 1));
        nNoStimGoodNonZeroTrRun = length(indNoStimGoodNonZeroTrRun);
        nElemNoStimGoodNonZeroTrRun = (nNoStimGoodNonZeroTrRun*nNoStimGoodNonZeroTrRun-nNoStimGoodNonZeroTrRun)/2;
        meanCorrTRun.meanNoStimGoodNZ(n) = sum(noStimGoodTrCorr(isnan(noStimGoodTrCorr) == 0))/nElemNoStimGoodNonZeroTrRun;
        meanCorrTRun.nNoStimGoodNonZeroTr(n) = nNoStimGoodNonZeroTrRun;
        
        trialNoNonStimBad = setdiff(trialNoNonStimBad,1:startTrNo);
        noStimBadTrCorr = triu(spikeCorrTRun{n}(trialNoNonStimBad,trialNoNonStimBad),1);
        noStimBadTrCorr = noStimBadTrCorr(:);
        nNoStimBad = length(trialNoNonStimBad);
        nElemNoStimBad = (nNoStimBad*nNoStimBad-nNoStimBad)/2;
        meanCorrTRun.meanNoStimBad(n) = sum(noStimBadTrCorr(isnan(noStimBadTrCorr) == 0))/nElemNoStimBad;
        indNoStimBadNonZeroTrRun = intersect(trialNoNonStimBad, find(nonZeroTrRun{n} == 1));
        nNoStimBadNonZeroTrRun = length(indNoStimBadNonZeroTrRun);
        nElemNoStimBadNonZeroTrRun = (nNoStimBadNonZeroTrRun*nNoStimBadNonZeroTrRun-nNoStimBadNonZeroTrRun)/2;
        meanCorrTRun.meanNoStimBadNZ(n) = sum(noStimBadTrCorr(isnan(noStimBadTrCorr) == 0))/nElemNoStimBadNonZeroTrRun;
        meanCorrTRun.nNoStimBadNonZeroTr(n) = nNoStimBadNonZeroTrRun;
        
        meanCorrTRun.meanStim = [];
        meanCorrTRun.meanStimNZ = [];
        meanCorrTRun.nStimNonZeroTr = [];
        meanCorrTRun.meanStimCtrl = [];
        meanCorrTRun.meanStimCtrlNZ = [];
        meanCorrTRun.nStimCtrlNonZeroTr = [];
        for i = 1:length(pulseMeth) 
            trialNoStim{i} = setdiff(trialNoStim{i},1:startTrNo);
            stimTrCorr = triu(spikeCorrTRun{n}(trialNoStim{i},trialNoStim{i}),1);
            stimTrCorr = stimTrCorr(:);
            nStim = length(trialNoStim{i});
            nElemStim = (nStim*nStim-nStim)/2;
            meanCorrTRun.meanStim{i}(n) = sum(stimTrCorr(isnan(stimTrCorr) == 0))/nElemStim;
            indStimNonZeroTrRun = intersect(trialNoStim{i}, find(nonZeroTrRun{n} == 1));
            nStimNonZeroTrRun = length(indStimNonZeroTrRun);
            nElemStimNonZeroTrRun = (nStimNonZeroTrRun*nStimNonZeroTrRun-nStimNonZeroTrRun)/2;
            meanCorrTRun.meanStimNZ{i}(n) = sum(stimTrCorr(isnan(stimTrCorr) == 0))/nElemStimNonZeroTrRun;
            meanCorrTRun.nStimNonZeroTr{i}(n) = nStimNonZeroTrRun;

            trialNoStimCtrl{i} = setdiff(trialNoStimCtrl{i},1:startTrNo);
            stimCtrlTrCorr = triu(spikeCorrTRun{n}(trialNoStimCtrl{i},trialNoStimCtrl{i}),1);
            stimCtrlTrCorr = stimCtrlTrCorr(:);
            nStimCtrl = length(trialNoStimCtrl{i});
            nElemStimCtrl = (nStimCtrl*nStimCtrl-nStimCtrl)/2;
            meanCorrTRun.meanStimCtrl{i}(n) = sum(stimCtrlTrCorr(isnan(stimCtrlTrCorr) == 0))/nElemStimCtrl;
            indStimCtrlNonZeroTrRun = intersect(trialNoStimCtrl{i}, find(nonZeroTrRun{n} == 1));
            nStimCtrlNonZeroTrRun = length(indStimCtrlNonZeroTrRun);
            nElemStimCtrlNonZeroTrRun = (nStimCtrlNonZeroTrRun*nStimCtrlNonZeroTrRun-nStimCtrlNonZeroTrRun)/2;
            meanCorrTRun.meanStimCtrlNZ{i}(n) = sum(stimCtrlTrCorr(isnan(stimCtrlTrCorr) == 0))/nElemStimCtrlNonZeroTrRun;
            meanCorrTRun.nStimCtrlNonZeroTr{i}(n) = nStimCtrlNonZeroTrRun;
        end
        %%
        
        %% added by YIngxue on 5/29/2021
        noStimGoodTrCorr = triu(spikeCorrTCue{n}(trialNoNonStimGood,trialNoNonStimGood),1);
        noStimGoodTrCorr = noStimGoodTrCorr(:);
        nNoStimGood = length(trialNoNonStimGood);
        nElemNoStimGood = (nNoStimGood*nNoStimGood-nNoStimGood)/2;
        meanCorrTCue.meanNoStimGood(n) = sum(noStimGoodTrCorr(isnan(noStimGoodTrCorr) == 0))/nElemNoStimGood;
        indNoStimGoodNonZeroTrCue = intersect(trialNoNonStimGood, find(nonZeroTrCue{n} == 1));
        nNoStimGoodNonZeroTrCue = length(indNoStimGoodNonZeroTrCue);
        nElemNoStimGoodNonZeroTrCue = (nNoStimGoodNonZeroTrCue*nNoStimGoodNonZeroTrCue-nNoStimGoodNonZeroTrCue)/2;
        meanCorrTCue.meanNoStimGoodNZ(n) = sum(noStimGoodTrCorr(isnan(noStimGoodTrCorr) == 0))/nElemNoStimGoodNonZeroTrCue;
        meanCorrTCue.nNoStimGoodNonZeroTr(n) = nNoStimGoodNonZeroTrCue;
        
        noStimBadTrCorr = triu(spikeCorrTCue{n}(trialNoNonStimBad,trialNoNonStimBad),1);
        noStimBadTrCorr = noStimBadTrCorr(:);
        nNoStimBad = length(trialNoNonStimBad);
        nElemNoStimBad = (nNoStimBad*nNoStimBad-nNoStimBad)/2;
        meanCorrTCue.meanNoStimBad(n) = sum(noStimBadTrCorr(isnan(noStimBadTrCorr) == 0))/nElemNoStimBad;
        indNoStimBadNonZeroTrCue = intersect(trialNoNonStimBad, find(nonZeroTrCue{n} == 1));
        nNoStimBadNonZeroTrCue = length(indNoStimBadNonZeroTrCue);
        nElemNoStimBadNonZeroTrCue = (nNoStimBadNonZeroTrCue*nNoStimBadNonZeroTrCue-nNoStimBadNonZeroTrCue)/2;
        meanCorrTCue.meanNoStimBadNZ(n) = sum(noStimBadTrCorr(isnan(noStimBadTrCorr) == 0))/nElemNoStimBadNonZeroTrCue;
        meanCorrTCue.nNoStimBadNonZeroTr(n) = nNoStimBadNonZeroTrCue;
        
        meanCorrTCue.meanStim = [];
        meanCorrTCue.meanStimNZ = [];
        meanCorrTCue.nStimNonZeroTr = [];
        meanCorrTCue.meanStimCtrl = [];
        meanCorrTCue.meanStimCtrlNZ = [];
        meanCorrTCue.nStimCtrlNonZeroTr = [];
        for i = 1:length(pulseMeth) 
            stimTrCorr = triu(spikeCorrTCue{n}(trialNoStim{i},trialNoStim{i}),1);
            stimTrCorr = stimTrCorr(:);
            nStim = length(trialNoStim{i});
            nElemStim = (nStim*nStim-nStim)/2;
            meanCorrTCue.meanStim{i}(n) = sum(stimTrCorr(isnan(stimTrCorr) == 0))/nElemStim;
            indStimNonZeroTrCue = intersect(trialNoStim{i}, find(nonZeroTrCue{n} == 1));
            nStimNonZeroTrCue = length(indStimNonZeroTrCue);
            nElemStimNonZeroTrCue = (nStimNonZeroTrCue*nStimNonZeroTrCue-nStimNonZeroTrCue)/2;
            meanCorrTCue.meanStimNZ{i}(n) = sum(stimTrCorr(isnan(stimTrCorr) == 0))/nElemStimNonZeroTrCue;
            meanCorrTCue.nStimNonZeroTr{i}(n) = nStimNonZeroTrCue;

            stimCtrlTrCorr = triu(spikeCorrTCue{n}(trialNoStimCtrl{i},trialNoStimCtrl{i}),1);
            stimCtrlTrCorr = stimCtrlTrCorr(:);
            nStimCtrl = length(trialNoStimCtrl{i});
            nElemStimCtrl = (nStimCtrl*nStimCtrl-nStimCtrl)/2;
            meanCorrTCue.meanStimCtrl{i}(n) = sum(stimCtrlTrCorr(isnan(stimCtrlTrCorr) == 0))/nElemStimCtrl;
            indStimCtrlNonZeroTrCue = intersect(trialNoStimCtrl{i}, find(nonZeroTrCue{n} == 1));
            nStimCtrlNonZeroTrCue = length(indStimCtrlNonZeroTrCue);
            nElemStimCtrlNonZeroTrCue = (nStimCtrlNonZeroTrCue*nStimCtrlNonZeroTrCue-nStimCtrlNonZeroTrCue)/2;
            meanCorrTCue.meanStimCtrlNZ{i}(n) = sum(stimCtrlTrCorr(isnan(stimCtrlTrCorr) == 0))/nElemStimCtrlNonZeroTrCue;
            meanCorrTCue.nStimCtrlNonZeroTr{i}(n) = nStimCtrlNonZeroTrCue;
        end
        %%
    end
    
    % compare single neuron correlation between good and bad trials --
    % aligned to run
    plotCompCorr(meanCorrTRun.meanGood,meanCorrTRun.meanBad,...
        'Neu Tr CorrT Run - good trials','Neu Tr CorrT Run - bad trials');
    
%     % compare single neuron correlation between good and bad trials --
%     % aligned to reward
%     plotCompCorr(meanCorrTRew.meanGood,meanCorrTRew.meanBad,...
%         'Neu Tr CorrT Rew - good trials','Neu Tr CorrT Rew - bad trials');
%         
%     % compare single neuron correlation between good and bad trials --
%     % aligned to cue
%     plotCompCorr(meanCorrTCue.meanGood,meanCorrTCue.meanBad,...
%         'Neu Tr CorrT Cue - good trials','Neu Tr CorrT Cue - bad trials');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward
%     plotCompCorr(meanCorrTRun.mean,meanCorrTRew.mean,...
%         'Neu Tr CorrT Run','Neu Tr CorrT Rew');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue
%     plotCompCorr(meanCorrTRun.mean,meanCorrTCue.mean,...
%         'Neu Tr CorrT Run','Neu Tr CorrT Cue');
%     
    % compare single neuron correlation between aligned to run and aligned
    % to reward for good trials
    plotCompCorr(meanCorrTRun.meanGood,meanCorrTRew.meanGood,...
        'Neu Tr CorrT Run - good trials','Neu Tr CorrT Rew - good trials');
   
    % compare single neuron correlation between aligned to run and aligned
    % to cue for good trials
    plotCompCorr(meanCorrTRun.meanGood,meanCorrTCue.meanGood,...
        'Neu Tr CorrT Run - good trials','Neu Tr CorrT Cue - good trials');
%     
%     %% only consider nonzero trials
%     % compare single neuron correlation between good and bad trials --
%     % aligned to run
%     indNeu = meanCorrTRun.nGoodNonZeroTr > 10 & meanCorrTRun.nBadNonZeroTr > 5;
%     plotCompCorr(meanCorrTRun.meanGoodNZ(indNeu),meanCorrTRun.meanBadNZ(indNeu),...
%         'Neu None-zero Tr CorrT Run - good trials','Neu None-zero Tr CorrT Run - bad trials');
%     
%     % compare single neuron correlation between good and bad trials --
%     % aligned to reward
%     indNeu = meanCorrTRew.nGoodNonZeroTr > 10 & meanCorrTRew.nBadNonZeroTr > 5;
%     plotCompCorr(meanCorrTRew.meanGoodNZ(indNeu),meanCorrTRew.meanBadNZ(indNeu),...
%         'Neu None-zero Tr CorrT Rew - good trials','Neu None-zero Tr CorrT Rew - bad trials');
%         
%     % compare single neuron correlation between good and bad trials --
%     % aligned to cue
%     indNeu = meanCorrTCue.nGoodNonZeroTr > 10 & meanCorrTCue.nBadNonZeroTr > 5;
%     plotCompCorr(meanCorrTCue.meanGoodNZ(indNeu),meanCorrTCue.meanBadNZ(indNeu),...
%         'Neu None-zero Tr CorrT Cue - good trials','Neu None-zero Tr CorrT Cue - bad trials');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to reward
%     indNeu = meanCorrTRun.nNonZeroTr > 10 & meanCorrTRew.nNonZeroTr > 10;
%     plotCompCorr(meanCorrTRun.meanNZ(indNeu),meanCorrTRew.meanNZ(indNeu),...
%         'Neu None-zero Tr CorrT Run ','Neu None-zero Tr CorrT Rew');
%     
%     % compare single neuron correlation between aligned to run and aligned
%     % to cue
%     indNeu = meanCorrTRun.nNonZeroTr > 10 & meanCorrTCue.nNonZeroTr > 10;
%     plotCompCorr(meanCorrTRun.meanNZ(indNeu),meanCorrTCue.meanNZ(indNeu),...
%         'Neu None-zero Tr CorrT Run','Neu None-zero Tr CorrT Cue');
%     
    % compare single neuron correlation between aligned to run and aligned
    % to reward for good trials
    indNeu = meanCorrTRun.nGoodNonZeroTr > 5 & meanCorrTRew.nGoodNonZeroTr > 5;
    if(sum(indNeu) > 0)
        plotCompCorr(meanCorrTRun.meanGoodNZ(indNeu),meanCorrTRew.meanGoodNZ(indNeu),...
            'Neu None-zero Tr CorrT Run - good trials','Neu None-zero Tr CorrT Rew - good trials');
    end
   
    % compare single neuron correlation between aligned to run and aligned
    % to cue for good trials
    indNeu = meanCorrTRun.nGoodNonZeroTr > 5 & meanCorrTCue.nGoodNonZeroTr > 5;
    if(sum(indNeu) > 0)
        plotCompCorr(meanCorrTRun.meanGoodNZ(indNeu),meanCorrTCue.meanGoodNZ(indNeu),...
            'Neu None-zero Tr CorrT Run - good trials','Neu None-zero Tr CorrT Cue - good trials');
    end
    
    save([path fileNameCorr],'meanCorrTRun','meanCorrTRew','meanCorrTCue');
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
