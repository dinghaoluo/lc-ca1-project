function meanNeuronSpikeCorrTStim(path,fileName,onlyRun,mazeSess,intervalT)
% single neuron level mean spike correlation across trials
% meanNeuronSpikeCorrT('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,10)

    if(nargin == 4)
        intervalT = 0;
    end
    fullPath = [path fileName '_spikesCorrTAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    fileNameCorr = [fileName '_meanSpikesCorrTAlignedStim_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
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
        
    neuronNo = length(spikeCorrTRun);
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
    
    meanCorrT = struct('meanNonStim',zeros(1,neuronNo),...
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
    
    meanCorrTRun = meanCorrT;
    meanCorrTRew = meanCorrT;
    meanCorrTCue = meanCorrT;
    for n = 1:neuronNo
        %% non-stim trials
        corrArr = triu(spikeCorrTRun{n}(trialNoNonStim,trialNoNonStim),1);
        corrArr = corrArr(:);        
        meanCorrTRun.meanNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStim;
        nonZeroTr = intersect(find(nonZeroTrRun{n} == 1),trialNoNonStim);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanCorrTRun.meanNZNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanCorrTRun.nNonZeroTrNonStim(n) = nNonZeroTr;
        
        corrArr = triu(spikeCorrTRew{n}(trialNoNonStim,trialNoNonStim),1);
        corrArr = corrArr(:);        
        meanCorrTRew.meanNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStim;
        nonZeroTr = intersect(find(nonZeroTrRew{n} == 1),trialNoNonStim);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanCorrTRew.meanNZNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanCorrTRew.nNonZeroTrNonStim(n) = nNonZeroTr;
        
        corrArr = triu(spikeCorrTCue{n}(trialNoNonStim,trialNoNonStim),1);
        corrArr = corrArr(:);        
        meanCorrTCue.meanNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStim;
        nonZeroTr = intersect(find(nonZeroTrCue{n} == 1),trialNoNonStim);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanCorrTCue.meanNZNonStim(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanCorrTCue.nNonZeroTrNonStim(n) = nNonZeroTr;
        
        %% non-stim good trials
        corrArr = triu(spikeCorrTRun{n}(trialNoNonStimGood,trialNoNonStimGood),1);
        corrArr = corrArr(:);        
        meanCorrTRun.meanNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimGood;
        nonZeroTr = intersect(find(nonZeroTrRun{n} == 1),trialNoNonStimGood);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanCorrTRun.meanNZNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanCorrTRun.nNonZeroTrNonStimGood(n) = nNonZeroTr;
        
        corrArr = triu(spikeCorrTRew{n}(trialNoNonStimGood,trialNoNonStimGood),1);
        corrArr = corrArr(:);        
        meanCorrTRew.meanNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimGood;
        nonZeroTr = intersect(find(nonZeroTrRew{n} == 1),trialNoNonStimGood);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanCorrTRew.meanNZNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanCorrTRew.nNonZeroTrNonStimGood(n) = nNonZeroTr;
        
        corrArr = triu(spikeCorrTCue{n}(trialNoNonStimGood,trialNoNonStimGood),1);
        corrArr = corrArr(:);        
        meanCorrTCue.meanNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimGood;
        nonZeroTr = intersect(find(nonZeroTrCue{n} == 1),trialNoNonStimGood);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanCorrTCue.meanNZNonStimGood(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanCorrTCue.nNonZeroTrNonStimGood(n) = nNonZeroTr;
        
        %% non-stim bad trials
        corrArr = triu(spikeCorrTRun{n}(trialNoNonStimBad,trialNoNonStimBad),1);
        corrArr = corrArr(:);        
        meanCorrTRun.meanNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimBad;
        nonZeroTr = intersect(find(nonZeroTrRun{n} == 1),trialNoNonStimBad);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanCorrTRun.meanNZNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanCorrTRun.nNonZeroTrNonStimBad(n) = nNonZeroTr;
        
        corrArr = triu(spikeCorrTRew{n}(trialNoNonStimBad,trialNoNonStimBad),1);
        corrArr = corrArr(:);        
        meanCorrTRew.meanNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimBad;
        nonZeroTr = intersect(find(nonZeroTrRew{n} == 1),trialNoNonStimBad);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanCorrTRew.meanNZNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanCorrTRew.nNonZeroTrNonStimBad(n) = nNonZeroTr;
        
        corrArr = triu(spikeCorrTCue{n}(trialNoNonStimBad,trialNoNonStimBad),1);
        corrArr = corrArr(:);        
        meanCorrTCue.meanNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimBad;
        nonZeroTr = intersect(find(nonZeroTrCue{n} == 1),trialNoNonStimBad);
        nNonZeroTr = length(nonZeroTr);
        nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanCorrTCue.meanNZNonStimBad(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
        meanCorrTCue.nNonZeroTrNonStimBad(n) = nNonZeroTr;
        
        for i = 1:length(pulseMeth)
            %% stim trials
            corrArr = triu(spikeCorrTRun{n}(trialNoStim{i},trialNoStim{i}),1);
            corrArr = corrArr(:);        
            meanCorrTRun.meanStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStim(i);
            nonZeroTr = intersect(find(nonZeroTrRun{n} == 1),trialNoStim{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanCorrTRun.meanNZStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanCorrTRun.nNonZeroTrStim(i,n) = nNonZeroTr;

            corrArr = triu(spikeCorrTRew{n}(trialNoStim{i},trialNoStim{i}),1);
            corrArr = corrArr(:);        
            meanCorrTRew.meanStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStim(i);
            nonZeroTr = intersect(find(nonZeroTrRew{n} == 1),trialNoStim{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanCorrTRew.meanNZStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanCorrTRew.nNonZeroTrStim(i,n) = nNonZeroTr;

            corrArr = triu(spikeCorrTCue{n}(trialNoStim{i},trialNoStim{i}),1);
            corrArr = corrArr(:);        
            meanCorrTCue.meanStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStim(i);
            nonZeroTr = intersect(find(nonZeroTrCue{n} == 1),trialNoStim{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanCorrTCue.meanNZStim(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanCorrTCue.nNonZeroTrStim(i,n) = nNonZeroTr;
            
            %% stim ctrl trials
            corrArr = triu(spikeCorrTRun{n}(trialNoStimCtrl{i},trialNoStimCtrl{i}),1);
            corrArr = corrArr(:);        
            meanCorrTRun.meanStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStimCtrl(i);
            nonZeroTr = intersect(find(nonZeroTrRun{n} == 1),trialNoStimCtrl{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanCorrTRun.meanNZStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanCorrTRun.nNonZeroTrStimCtrl(i,n) = nNonZeroTr;

            corrArr = triu(spikeCorrTRew{n}(trialNoStimCtrl{i},trialNoStimCtrl{i}),1);
            corrArr = corrArr(:);        
            meanCorrTRew.meanStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStimCtrl(i);
            nonZeroTr = intersect(find(nonZeroTrRew{n} == 1),trialNoStimCtrl{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanCorrTRew.meanNZStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanCorrTRew.nNonZeroTrStimCtrl(i,n) = nNonZeroTr;

            corrArr = triu(spikeCorrTCue{n}(trialNoStimCtrl{i},trialNoStimCtrl{i}),1);
            corrArr = corrArr(:);        
            meanCorrTCue.meanStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemStimCtrl(i);
            nonZeroTr = intersect(find(nonZeroTrCue{n} == 1),trialNoStimCtrl{i});
            nNonZeroTr = length(nonZeroTr);
            nElemNonZeroRun = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanCorrTCue.meanNZStimCtrl(i,n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZeroRun;
            meanCorrTCue.nNonZeroTrStimCtrl(i,n) = nNonZeroTr;
        end
        
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
