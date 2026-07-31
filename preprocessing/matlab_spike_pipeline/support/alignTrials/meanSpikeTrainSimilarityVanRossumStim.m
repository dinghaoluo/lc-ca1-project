function meanSpikeTrainSimilarityVanRossumStim(path,fileName,onlyRun,mazeSess,tc,intervalT)
% single neuron level mean spike correlation across trials
% meanSpikeTrainSimilarityVPIStim('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,100,8)

    fullPath = [path fileName '_spikeTrainSimiVanR_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_tc' num2str(tc) '_intT' ...
            num2str(intervalT) '.mat'];
    fileNameCorr = [fileName '_meanSpikeTrainSimVanRStim_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_tc' num2str(tc) '_intT' ...
            num2str(intervalT) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _spikeTrainSimiVR_Run file does not exist');
        return;
    end
    load(fullPath);
    
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
        
    neuronNo = length(spikeTrainSimVRRun);
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
    
    meanSimVR = struct('meanNonStim',zeros(1,neuronNo),...
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
    
    meanSimVRRun = meanSimVR;
    meanSimVRRew = meanSimVR;
    meanSimVRCue = meanSimVR;
    for n = 1:neuronNo
        %% non-stim trials
        corrArr = triu(spikeTrainSimVRRun{n}(trialNoNonStim,trialNoNonStim),1);
        corrArr = corrArr(:);        
        meanSimVRRun.meanNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStim;
        nonZeroTr = intersect(find(nonZeroTrRun{n} == 1),trialNoNonStim);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVRRun.meanNZNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimVRRun.nNonZeroTrNonStim(n) = nNonZeroTr;
        
        corrArr = triu(spikeTrainSimVRRew{n}(trialNoNonStim,trialNoNonStim),1);
        corrArr = corrArr(:);        
        meanSimVRRew.meanNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStim;
        nonZeroTr = intersect(find(nonZeroTrRew{n} == 1),trialNoNonStim);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVRRew.meanNZNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimVRRew.nNonZeroTrNonStim(n) = nNonZeroTr;
        
        corrArr = triu(spikeTrainSimVRCue{n}(trialNoNonStim,trialNoNonStim),1);
        corrArr = corrArr(:);        
        meanSimVRCue.meanNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStim;
        nonZeroTr = intersect(find(nonZeroTrCue{n} == 1),trialNoNonStim);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVRCue.meanNZNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimVRCue.nNonZeroTrNonStim(n) = nNonZeroTr;
        
        %% non-stim good trials
        corrArr = triu(spikeTrainSimVRRun{n}(trialNoNonStimGood,trialNoNonStimGood),1);
        corrArr = corrArr(:);        
        meanSimVRRun.meanNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimGood;
        nonZeroTr = intersect(find(nonZeroTrRun{n} == 1),trialNoNonStimGood);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVRRun.meanNZNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimVRRun.nNonZeroTrNonStimGood(n) = nNonZeroTr;
        
        corrArr = triu(spikeTrainSimVRRew{n}(trialNoNonStimGood,trialNoNonStimGood),1);
        corrArr = corrArr(:);        
        meanSimVRRew.meanNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimGood;
        nonZeroTr = intersect(find(nonZeroTrRew{n} == 1),trialNoNonStimGood);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVRRew.meanNZNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimVRRew.nNonZeroTrNonStimGood(n) = nNonZeroTr;
        
        corrArr = triu(spikeTrainSimVRCue{n}(trialNoNonStimGood,trialNoNonStimGood),1);
        corrArr = corrArr(:);        
        meanSimVRCue.meanNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimGood;
        nonZeroTr = intersect(find(nonZeroTrCue{n} == 1),trialNoNonStimGood);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVRCue.meanNZNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimVRCue.nNonZeroTrNonStimGood(n) = nNonZeroTr;
        
        %% non-stim bad trials
        corrArr = triu(spikeTrainSimVRRun{n}(trialNoNonStimBad,trialNoNonStimBad),1);
        corrArr = corrArr(:);        
        meanSimVRRun.meanNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimBad;
        nonZeroTr = intersect(find(nonZeroTrRun{n} == 1),trialNoNonStimBad);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVRRun.meanNZNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimVRRun.nNonZeroTrNonStimBad(n) = nNonZeroTr;
        
        corrArr = triu(spikeTrainSimVRRew{n}(trialNoNonStimBad,trialNoNonStimBad),1);
        corrArr = corrArr(:);        
        meanSimVRRew.meanNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimBad;
        nonZeroTr = intersect(find(nonZeroTrRew{n} == 1),trialNoNonStimBad);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVRRew.meanNZNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimVRRew.nNonZeroTrNonStimBad(n) = nNonZeroTr;
        
        corrArr = triu(spikeTrainSimVRCue{n}(trialNoNonStimBad,trialNoNonStimBad),1);
        corrArr = corrArr(:);        
        meanSimVRCue.meanNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimBad;
        nonZeroTr = intersect(find(nonZeroTrCue{n} == 1),trialNoNonStimBad);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimVRCue.meanNZNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimVRCue.nNonZeroTrNonStimBad(n) = nNonZeroTr;
        
        for i = 1:length(pulseMeth)
            %% stim trials
            corrArr = triu(spikeTrainSimVRRun{n}(trialNoStim{i},trialNoStim{i}),1);
            corrArr = corrArr(:);        
            meanSimVRRun.meanStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStim(i);
            nonZeroTr = intersect(find(nonZeroTrRun{n} == 1),trialNoStim{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanSimVRRun.meanNZStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanSimVRRun.nNonZeroTrStim(i,n) = nNonZeroTr;

            corrArr = triu(spikeTrainSimVRRew{n}(trialNoStim{i},trialNoStim{i}),1);
            corrArr = corrArr(:);        
            meanSimVRRew.meanStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStim(i);
            nonZeroTr = intersect(find(nonZeroTrRew{n} == 1),trialNoStim{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanSimVRRew.meanNZStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanSimVRRew.nNonZeroTrStim(i,n) = nNonZeroTr;

            corrArr = triu(spikeTrainSimVRCue{n}(trialNoStim{i},trialNoStim{i}),1);
            corrArr = corrArr(:);        
            meanSimVRCue.meanStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStim(i);
            nonZeroTr = intersect(find(nonZeroTrCue{n} == 1),trialNoStim{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanSimVRCue.meanNZStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanSimVRCue.nNonZeroTrStim(i,n) = nNonZeroTr;
            
            %% stim ctrl trials
            corrArr = triu(spikeTrainSimVRRun{n}(trialNoStimCtrl{i},trialNoStimCtrl{i}),1);
            corrArr = corrArr(:);        
            meanSimVRRun.meanStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStimCtrl(i);
            nonZeroTr = intersect(find(nonZeroTrRun{n} == 1),trialNoStimCtrl{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanSimVRRun.meanNZStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanSimVRRun.nNonZeroTrStimCtrl(i,n) = nNonZeroTr;

            corrArr = triu(spikeTrainSimVRRew{n}(trialNoStimCtrl{i},trialNoStimCtrl{i}),1);
            corrArr = corrArr(:);        
            meanSimVRRew.meanStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStimCtrl(i);
            nonZeroTr = intersect(find(nonZeroTrRew{n} == 1),trialNoStimCtrl{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanSimVRRew.meanNZStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanSimVRRew.nNonZeroTrStimCtrl(i,n) = nNonZeroTr;

            corrArr = triu(spikeTrainSimVRCue{n}(trialNoStimCtrl{i},trialNoStimCtrl{i}),1);
            corrArr = corrArr(:);        
            meanSimVRCue.meanStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStimCtrl(i);
            nonZeroTr = intersect(find(nonZeroTrCue{n} == 1),trialNoStimCtrl{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanSimVRCue.meanNZStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanSimVRCue.nNonZeroTrStimCtrl(i,n) = nNonZeroTr;
        end
        
    end
    
    save([path fileNameCorr],'meanSimVRRun','meanSimVRRew','meanSimVRCue');
end
