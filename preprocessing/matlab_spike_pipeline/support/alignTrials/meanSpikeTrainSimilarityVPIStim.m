function meanSpikeTrainSimilarityVPIStim(path,fileName,onlyRun,mazeSess,cost,intervalT)
% single neuron level mean spike correlation across trials
% meanSpikeTrainSimilarityVPIStim('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,0.016,8)

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
    load(fullPath,'spikeTrainSimVPIRun','spikeTrainSimVPIRew','spikeTrainSimVPICue',...
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
        
    neuronNo = length(spikeTrainSimVPIRun);
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
    
    meanSimVPI = struct('meanNonStim',zeros(1,neuronNo),...
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
    
    meanSimVPIRun = meanSimVPI;
    meanSimVPIRew = meanSimVPI;
    meanSimVPICue = meanSimVPI;
    for n = 1:neuronNo
        %% non-stim trials
        corrArr = triu(spikeTrainSimVPIRun{n}(trialNoNonStim,trialNoNonStim),1);
        corrArr = corrArr(:);        
        meanSimVPIRun.meanNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStim;
        nonZeroTr = intersect(find(nonZeroTrRun{n} == 1),trialNoNonStim);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPIRun.meanNZNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimVPIRun.nNonZeroTrNonStim(n) = nNonZeroTr;
        
        corrArr = triu(spikeTrainSimVPIRew{n}(trialNoNonStim,trialNoNonStim),1);
        corrArr = corrArr(:);        
        meanSimVPIRew.meanNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStim;
        nonZeroTr = intersect(find(nonZeroTrRew{n} == 1),trialNoNonStim);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPIRew.meanNZNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimVPIRew.nNonZeroTrNonStim(n) = nNonZeroTr;
        
        corrArr = triu(spikeTrainSimVPICue{n}(trialNoNonStim,trialNoNonStim),1);
        corrArr = corrArr(:);        
        meanSimVPICue.meanNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStim;
        nonZeroTr = intersect(find(nonZeroTrCue{n} == 1),trialNoNonStim);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPICue.meanNZNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimVPICue.nNonZeroTrNonStim(n) = nNonZeroTr;
        
        %% non-stim good trials
        corrArr = triu(spikeTrainSimVPIRun{n}(trialNoNonStimGood,trialNoNonStimGood),1);
        corrArr = corrArr(:);        
        meanSimVPIRun.meanNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimGood;
        nonZeroTr = intersect(find(nonZeroTrRun{n} == 1),trialNoNonStimGood);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPIRun.meanNZNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimVPIRun.nNonZeroTrNonStimGood(n) = nNonZeroTr;
        
        corrArr = triu(spikeTrainSimVPIRew{n}(trialNoNonStimGood,trialNoNonStimGood),1);
        corrArr = corrArr(:);        
        meanSimVPIRew.meanNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimGood;
        nonZeroTr = intersect(find(nonZeroTrRew{n} == 1),trialNoNonStimGood);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPIRew.meanNZNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimVPIRew.nNonZeroTrNonStimGood(n) = nNonZeroTr;
        
        corrArr = triu(spikeTrainSimVPICue{n}(trialNoNonStimGood,trialNoNonStimGood),1);
        corrArr = corrArr(:);        
        meanSimVPICue.meanNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimGood;
        nonZeroTr = intersect(find(nonZeroTrCue{n} == 1),trialNoNonStimGood);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPICue.meanNZNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimVPICue.nNonZeroTrNonStimGood(n) = nNonZeroTr;
        
        %% non-stim bad trials
        corrArr = triu(spikeTrainSimVPIRun{n}(trialNoNonStimBad,trialNoNonStimBad),1);
        corrArr = corrArr(:);        
        meanSimVPIRun.meanNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimBad;
        nonZeroTr = intersect(find(nonZeroTrRun{n} == 1),trialNoNonStimBad);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPIRun.meanNZNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimVPIRun.nNonZeroTrNonStimBad(n) = nNonZeroTr;
        
        corrArr = triu(spikeTrainSimVPIRew{n}(trialNoNonStimBad,trialNoNonStimBad),1);
        corrArr = corrArr(:);        
        meanSimVPIRew.meanNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimBad;
        nonZeroTr = intersect(find(nonZeroTrRew{n} == 1),trialNoNonStimBad);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPIRew.meanNZNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimVPIRew.nNonZeroTrNonStimBad(n) = nNonZeroTr;
        
        corrArr = triu(spikeTrainSimVPICue{n}(trialNoNonStimBad,trialNoNonStimBad),1);
        corrArr = corrArr(:);        
        meanSimVPICue.meanNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimBad;
        nonZeroTr = intersect(find(nonZeroTrCue{n} == 1),trialNoNonStimBad);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVPICue.meanNZNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimVPICue.nNonZeroTrNonStimBad(n) = nNonZeroTr;
        
        for i = 1:length(pulseMeth)
            %% stim trials
            corrArr = triu(spikeTrainSimVPIRun{n}(trialNoStim{i},trialNoStim{i}),1);
            corrArr = corrArr(:);        
            meanSimVPIRun.meanStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStim(i);
            nonZeroTr = intersect(find(nonZeroTrRun{n} == 1),trialNoStim{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanSimVPIRun.meanNZStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanSimVPIRun.nNonZeroTrStim(i,n) = nNonZeroTr;

            corrArr = triu(spikeTrainSimVPIRew{n}(trialNoStim{i},trialNoStim{i}),1);
            corrArr = corrArr(:);        
            meanSimVPIRew.meanStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStim(i);
            nonZeroTr = intersect(find(nonZeroTrRew{n} == 1),trialNoStim{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanSimVPIRew.meanNZStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanSimVPIRew.nNonZeroTrStim(i,n) = nNonZeroTr;

            corrArr = triu(spikeTrainSimVPICue{n}(trialNoStim{i},trialNoStim{i}),1);
            corrArr = corrArr(:);        
            meanSimVPICue.meanStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStim(i);
            nonZeroTr = intersect(find(nonZeroTrCue{n} == 1),trialNoStim{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanSimVPICue.meanNZStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanSimVPICue.nNonZeroTrStim(i,n) = nNonZeroTr;
            
            %% stim ctrl trials
            corrArr = triu(spikeTrainSimVPIRun{n}(trialNoStimCtrl{i},trialNoStimCtrl{i}),1);
            corrArr = corrArr(:);        
            meanSimVPIRun.meanStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStimCtrl(i);
            nonZeroTr = intersect(find(nonZeroTrRun{n} == 1),trialNoStimCtrl{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanSimVPIRun.meanNZStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanSimVPIRun.nNonZeroTrStimCtrl(i,n) = nNonZeroTr;

            corrArr = triu(spikeTrainSimVPIRew{n}(trialNoStimCtrl{i},trialNoStimCtrl{i}),1);
            corrArr = corrArr(:);        
            meanSimVPIRew.meanStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStimCtrl(i);
            nonZeroTr = intersect(find(nonZeroTrRew{n} == 1),trialNoStimCtrl{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanSimVPIRew.meanNZStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanSimVPIRew.nNonZeroTrStimCtrl(i,n) = nNonZeroTr;

            corrArr = triu(spikeTrainSimVPICue{n}(trialNoStimCtrl{i},trialNoStimCtrl{i}),1);
            corrArr = corrArr(:);        
            meanSimVPICue.meanStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStimCtrl(i);
            nonZeroTr = intersect(find(nonZeroTrCue{n} == 1),trialNoStimCtrl{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanSimVPICue.meanNZStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanSimVPICue.nNonZeroTrStimCtrl(i,n) = nNonZeroTr;
        end
        
    end
    
    save([path fileNameCorr],'meanSimVPIRun','meanSimVPIRew','meanSimVPICue');
end
