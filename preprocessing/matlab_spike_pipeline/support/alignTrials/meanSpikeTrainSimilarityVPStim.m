function meanSpikeTrainSimilarityVPStim(path,fileName,onlyRun,mazeSess,cost,intervalT)
% single neuron level mean spike correlation across trials
% meanSpikeTrainSimilarityVPStim('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,0.016,8)

    s = num2str(cost);
    ind = strfind(s,'.');
    s(ind) = 'p';
    fullPath = [path fileName '_spikeTrainSimVP_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_q' s '_intT' ...
            num2str(intervalT) '.mat'];
    fileNameCorr = [fileName '_meanSpikeTrainSimVPStim_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_q' s '_intT' ...
            num2str(intervalT) '.mat'];
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
        
    neuronNo = length(spikeTrainSimVPRun);
    nPulseMeth = length(pulseMeth);
    trialNoNonStim = [trialNoNonStimGood; trialNoNonStimBad];
    nTrialNonStim = length(trialNoNonStim);
    nElemNonStim = (nTrialNonStim*nTrialNonStim-nTrialNonStim)/2;
    
    nTrialNonStimGood = length(trialNoNonStimGood);
    nElemNonStimGood = (nTrialNonStimGood*nTrialNonStimGood-nTrialNonStimGood)/2;
    
    nTrialNonStimBad = length(trialNoNonStimBad);
    nElemNonStimBad = (nTrialNonStimBad*nTrialNonStimBad-nTrialNonStimBad)/2;
    
    nTrialStim = zeros(1,length(pulseMeth));
    nTrialStimCtrl = zeros(1,length(pulseMeth));
    nElemStim = zeros(1,length(pulseMeth));
    nElemStimCtrl = zeros(1,length(pulseMeth));
    for i = 1:length(pulseMeth)
        nTrialStim(i) = length(trialNoStim{i});
        nElemStim(i) = (nTrialStim(i)*nTrialStim(i)-nTrialStim(i))/2;
        
        nTrialStimCtrl(i) = length(trialNoStimCtrl{i});
        nElemStimCtrl(i) = (nTrialStimCtrl(i)*nTrialStimCtrl(i)-nTrialStimCtrl(i))/2;
    end
    
    meanSimVP = struct('meanNonStim',zeros(1,neuronNo),...
                       'meanNZNonStim',zeros(1,neuronNo),...
                       'nNonZeroTrNonStim',zeros(1,neuronNo),...
                       ...
                       'meanNonStimGood',zeros(1,neuronNo),...
                       'meanNZNonStimGood',zeros(1,neuronNo),...
                       'nNonZeroTrNonStimGood',zeros(1,neuronNo),...
                       ...
                       'meanNonStimBad',zeros(1,neuronNo),...
                       'meanNZNonStimBad',zeros(1,neuronNo),...
                       'nNonZeroTrNonStimBad',zeros(1,neuronNo),...
                       ...
                       'meanStim',zeros(nPulseMeth,neuronNo),...
                       'meanNZStim',zeros(nPulseMeth,neuronNo),...
                       'nNonZeroTrStim',zeros(nPulseMeth,neuronNo),...
                       ...
                       'meanStimCtrl',zeros(nPulseMeth,neuronNo),...
                       'meanNZStimCtrl',zeros(nPulseMeth,neuronNo),...
                       'nNonZeroTrStimCtrl',zeros(nPulseMeth,neuronNo));
    
    meanSimVPRun = meanSimVP;
    meanSimVPRew = meanSimVP;
    meanSimVPCue = meanSimVP;
    for n = 1:neuronNo
        %% non-stim trials
        corrArr = triu(spikeTrainSimVPRun{n}(trialNoNonStim,trialNoNonStim),1);
        corrArr = corrArr(:);        
        meanSimVPRun.meanNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStim;
        nonZeroTr = intersect(find(nonZeroTrRun{n} == 1),trialNoNonStim);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPRun.meanNZNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimVPRun.nNonZeroTrNonStim(n) = nNonZeroTr;
        
        corrArr = triu(spikeTrainSimVPRew{n}(trialNoNonStim,trialNoNonStim),1);
        corrArr = corrArr(:);        
        meanSimVPRew.meanNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStim;
        nonZeroTr = intersect(find(nonZeroTrRew{n} == 1),trialNoNonStim);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPRew.meanNZNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimVPRew.nNonZeroTrNonStim(n) = nNonZeroTr;
        
        corrArr = triu(spikeTrainSimVPCue{n}(trialNoNonStim,trialNoNonStim),1);
        corrArr = corrArr(:);        
        meanSimVPCue.meanNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStim;
        nonZeroTr = intersect(find(nonZeroTrCue{n} == 1),trialNoNonStim);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPCue.meanNZNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimVPCue.nNonZeroTrNonStim(n) = nNonZeroTr;
        
        %% non-stim good trials
        corrArr = triu(spikeTrainSimVPRun{n}(trialNoNonStimGood,trialNoNonStimGood),1);
        corrArr = corrArr(:);        
        meanSimVPRun.meanNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimGood;
        nonZeroTr = intersect(find(nonZeroTrRun{n} == 1),trialNoNonStimGood);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPRun.meanNZNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimVPRun.nNonZeroTrNonStimGood(n) = nNonZeroTr;
        
        corrArr = triu(spikeTrainSimVPRew{n}(trialNoNonStimGood,trialNoNonStimGood),1);
        corrArr = corrArr(:);        
        meanSimVPRew.meanNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimGood;
        nonZeroTr = intersect(find(nonZeroTrRew{n} == 1),trialNoNonStimGood);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPRew.meanNZNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimVPRew.nNonZeroTrNonStimGood(n) = nNonZeroTr;
        
        corrArr = triu(spikeTrainSimVPCue{n}(trialNoNonStimGood,trialNoNonStimGood),1);
        corrArr = corrArr(:);        
        meanSimVPCue.meanNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimGood;
        nonZeroTr = intersect(find(nonZeroTrCue{n} == 1),trialNoNonStimGood);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPCue.meanNZNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimVPCue.nNonZeroTrNonStimGood(n) = nNonZeroTr;
        
        %% non-stim bad trials
        corrArr = triu(spikeTrainSimVPRun{n}(trialNoNonStimBad,trialNoNonStimBad),1);
        corrArr = corrArr(:);        
        meanSimVPRun.meanNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimBad;
        nonZeroTr = intersect(find(nonZeroTrRun{n} == 1),trialNoNonStimBad);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPRun.meanNZNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimVPRun.nNonZeroTrNonStimBad(n) = nNonZeroTr;
        
        corrArr = triu(spikeTrainSimVPRew{n}(trialNoNonStimBad,trialNoNonStimBad),1);
        corrArr = corrArr(:);        
        meanSimVPRew.meanNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimBad;
        nonZeroTr = intersect(find(nonZeroTrRew{n} == 1),trialNoNonStimBad);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPRew.meanNZNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimVPRew.nNonZeroTrNonStimBad(n) = nNonZeroTr;
        
        corrArr = triu(spikeTrainSimVPCue{n}(trialNoNonStimBad,trialNoNonStimBad),1);
        corrArr = corrArr(:);        
        meanSimVPCue.meanNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimBad;
        nonZeroTr = intersect(find(nonZeroTrCue{n} == 1),trialNoNonStimBad);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPCue.meanNZNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimVPCue.nNonZeroTrNonStimBad(n) = nNonZeroTr;
        
        for i = 1:length(pulseMeth)
            %% stim trials
            corrArr = triu(spikeTrainSimVPRun{n}(trialNoStim{i},trialNoStim{i}),1);
            corrArr = corrArr(:);        
            meanSimVPRun.meanStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStim(i);
            nonZeroTr = intersect(find(nonZeroTrRun{n} == 1),trialNoStim{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanSimVPRun.meanNZStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanSimVPRun.nNonZeroTrStim(i,n) = nNonZeroTr;

            corrArr = triu(spikeTrainSimVPRew{n}(trialNoStim{i},trialNoStim{i}),1);
            corrArr = corrArr(:);        
            meanSimVPRew.meanStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStim(i);
            nonZeroTr = intersect(find(nonZeroTrRew{n} == 1),trialNoStim{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanSimVPRew.meanNZStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanSimVPRew.nNonZeroTrStim(i,n) = nNonZeroTr;

            corrArr = triu(spikeTrainSimVPCue{n}(trialNoStim{i},trialNoStim{i}),1);
            corrArr = corrArr(:);        
            meanSimVPCue.meanStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStim(i);
            nonZeroTr = intersect(find(nonZeroTrCue{n} == 1),trialNoStim{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanSimVPCue.meanNZStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanSimVPCue.nNonZeroTrStim(i,n) = nNonZeroTr;
            
            %% stim ctrl trials
            corrArr = triu(spikeTrainSimVPRun{n}(trialNoStimCtrl{i},trialNoStimCtrl{i}),1);
            corrArr = corrArr(:);        
            meanSimVPRun.meanStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStimCtrl(i);
            nonZeroTr = intersect(find(nonZeroTrRun{n} == 1),trialNoStimCtrl{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanSimVPRun.meanNZStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanSimVPRun.nNonZeroTrStimCtrl(i,n) = nNonZeroTr;

            corrArr = triu(spikeTrainSimVPRew{n}(trialNoStimCtrl{i},trialNoStimCtrl{i}),1);
            corrArr = corrArr(:);        
            meanSimVPRew.meanStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStimCtrl(i);
            nonZeroTr = intersect(find(nonZeroTrRew{n} == 1),trialNoStimCtrl{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanSimVPRew.meanNZStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanSimVPRew.nNonZeroTrStimCtrl(i,n) = nNonZeroTr;

            corrArr = triu(spikeTrainSimVPCue{n}(trialNoStimCtrl{i},trialNoStimCtrl{i}),1);
            corrArr = corrArr(:);        
            meanSimVPCue.meanStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStimCtrl(i);
            nonZeroTr = intersect(find(nonZeroTrCue{n} == 1),trialNoStimCtrl{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanSimVPCue.meanNZStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanSimVPCue.nNonZeroTrStimCtrl(i,n) = nNonZeroTr;
        end
        
    end
    
    save([path fileNameCorr],'meanSimVPRun','meanSimVPRew','meanSimVPCue');
end
