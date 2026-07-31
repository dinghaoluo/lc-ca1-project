function meanSpikeTrainSimilarityTStim(path,fileName,onlyRun,mazeSess,intervalT)
% single neuron level mean spike correlation across trials
% meanNeuronSpikeCorrT('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,10)

    if(nargin == 4)
        intervalT = 0;
    end
    fullPath = [path fileName '_spikeTrainSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    fileNameCorr = [fileName '_meanSpikeTrainSimTStim_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
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
    
    fullPath = [path fileName '_PeakFRAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];    
    if(exist(fullPath) == 0)
        disp('The _PeakFRAligned file does not exist');
        return;
    end
    load(fullPath,'trialNoNonStimGood','trialNoNonStimBad','trialNoStim','trialNoStimCtrl','pulseMeth');
        
    neuronNo = length(spikeTrainSimTRun);
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
    
    meanSimT = struct('meanNonStim',zeros(1,neuronNo),...
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
    
    meanSimTRun = meanSimT;
    meanSimTRew = meanSimT;
    meanSimTCue = meanSimT;
    for n = 1:neuronNo
        %% non-stim trials
        corrArr = triu(spikeTrainSimTRun{n}(trialNoNonStim,trialNoNonStim),1);
        corrArr = corrArr(:);        
        meanSimTRun.meanNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStim;
        nonZeroTr = intersect(find(nonZeroTrRun{n} == 1),trialNoNonStim);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimTRun.meanNZNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimTRun.nNonZeroTrNonStim(n) = nNonZeroTr;
        
        corrArr = triu(spikeTrainSimTRew{n}(trialNoNonStim,trialNoNonStim),1);
        corrArr = corrArr(:);        
        meanSimTRew.meanNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStim;
        nonZeroTr = intersect(find(nonZeroTrRew{n} == 1),trialNoNonStim);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimTRew.meanNZNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimTRew.nNonZeroTrNonStim(n) = nNonZeroTr;
        
        corrArr = triu(spikeTrainSimTCue{n}(trialNoNonStim,trialNoNonStim),1);
        corrArr = corrArr(:);        
        meanSimTCue.meanNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStim;
        nonZeroTr = intersect(find(nonZeroTrCue{n} == 1),trialNoNonStim);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimTCue.meanNZNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimTCue.nNonZeroTrNonStim(n) = nNonZeroTr;
        
        %% non-stim good trials
        corrArr = triu(spikeTrainSimTRun{n}(trialNoNonStimGood,trialNoNonStimGood),1);
        corrArr = corrArr(:);        
        meanSimTRun.meanNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimGood;
        nonZeroTr = intersect(find(nonZeroTrRun{n} == 1),trialNoNonStimGood);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimTRun.meanNZNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimTRun.nNonZeroTrNonStimGood(n) = nNonZeroTr;
        
        corrArr = triu(spikeTrainSimTRew{n}(trialNoNonStimGood,trialNoNonStimGood),1);
        corrArr = corrArr(:);        
        meanSimTRew.meanNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimGood;
        nonZeroTr = intersect(find(nonZeroTrRew{n} == 1),trialNoNonStimGood);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimTRew.meanNZNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimTRew.nNonZeroTrNonStimGood(n) = nNonZeroTr;
        
        corrArr = triu(spikeTrainSimTCue{n}(trialNoNonStimGood,trialNoNonStimGood),1);
        corrArr = corrArr(:);        
        meanSimTCue.meanNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimGood;
        nonZeroTr = intersect(find(nonZeroTrCue{n} == 1),trialNoNonStimGood);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimTCue.meanNZNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimTCue.nNonZeroTrNonStimGood(n) = nNonZeroTr;
        
        %% non-stim bad trials
        corrArr = triu(spikeTrainSimTRun{n}(trialNoNonStimBad,trialNoNonStimBad),1);
        corrArr = corrArr(:);        
        meanSimTRun.meanNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimBad;
        nonZeroTr = intersect(find(nonZeroTrRun{n} == 1),trialNoNonStimBad);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimTRun.meanNZNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimTRun.nNonZeroTrNonStimBad(n) = nNonZeroTr;
        
        corrArr = triu(spikeTrainSimTRew{n}(trialNoNonStimBad,trialNoNonStimBad),1);
        corrArr = corrArr(:);        
        meanSimTRew.meanNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimBad;
        nonZeroTr = intersect(find(nonZeroTrRew{n} == 1),trialNoNonStimBad);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimTRew.meanNZNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimTRew.nNonZeroTrNonStimBad(n) = nNonZeroTr;
        
        corrArr = triu(spikeTrainSimTCue{n}(trialNoNonStimBad,trialNoNonStimBad),1);
        corrArr = corrArr(:);        
        meanSimTCue.meanNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimBad;
        nonZeroTr = intersect(find(nonZeroTrCue{n} == 1),trialNoNonStimBad);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanSimTCue.meanNZNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanSimTCue.nNonZeroTrNonStimBad(n) = nNonZeroTr;
        
        for i = 1:length(pulseMeth)
            %% stim trials
            corrArr = triu(spikeTrainSimTRun{n}(trialNoStim{i},trialNoStim{i}),1);
            corrArr = corrArr(:);        
            meanSimTRun.meanStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStim(i);
            nonZeroTr = intersect(find(nonZeroTrRun{n} == 1),trialNoStim{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanSimTRun.meanNZStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanSimTRun.nNonZeroTrStim(i,n) = nNonZeroTr;

            corrArr = triu(spikeTrainSimTRew{n}(trialNoStim{i},trialNoStim{i}),1);
            corrArr = corrArr(:);        
            meanSimTRew.meanStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStim(i);
            nonZeroTr = intersect(find(nonZeroTrRew{n} == 1),trialNoStim{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanSimTRew.meanNZStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanSimTRew.nNonZeroTrStim(i,n) = nNonZeroTr;

            corrArr = triu(spikeTrainSimTCue{n}(trialNoStim{i},trialNoStim{i}),1);
            corrArr = corrArr(:);        
            meanSimTCue.meanStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStim(i);
            nonZeroTr = intersect(find(nonZeroTrCue{n} == 1),trialNoStim{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanSimTCue.meanNZStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanSimTCue.nNonZeroTrStim(i,n) = nNonZeroTr;
            
            %% stim ctrl trials
            corrArr = triu(spikeTrainSimTRun{n}(trialNoStimCtrl{i},trialNoStimCtrl{i}),1);
            corrArr = corrArr(:);        
            meanSimTRun.meanStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStimCtrl(i);
            nonZeroTr = intersect(find(nonZeroTrRun{n} == 1),trialNoStimCtrl{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanSimTRun.meanNZStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanSimTRun.nNonZeroTrStimCtrl(i,n) = nNonZeroTr;

            corrArr = triu(spikeTrainSimTRew{n}(trialNoStimCtrl{i},trialNoStimCtrl{i}),1);
            corrArr = corrArr(:);        
            meanSimTRew.meanStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStimCtrl(i);
            nonZeroTr = intersect(find(nonZeroTrRew{n} == 1),trialNoStimCtrl{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanSimTRew.meanNZStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanSimTRew.nNonZeroTrStimCtrl(i,n) = nNonZeroTr;

            corrArr = triu(spikeTrainSimTCue{n}(trialNoStimCtrl{i},trialNoStimCtrl{i}),1);
            corrArr = corrArr(:);        
            meanSimTCue.meanStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStimCtrl(i);
            nonZeroTr = intersect(find(nonZeroTrCue{n} == 1),trialNoStimCtrl{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanSimTCue.meanNZStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanSimTCue.nNonZeroTrStimCtrl(i,n) = nNonZeroTr;
        end
        
    end
    
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
