function populationSimilarityVP(path,fileName,onlyRun,costq,costk,intervalT)
% single neuron level spike correlation across trials
% Victor & Purpura spike time distance and Victor & Purpura spike time
% interval
% e.g.: spikeTrainSimilarityVP('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,0.016,10,8)
% cost: cost per unit time to move a spike
% intervalT: spike train max time in sec
    
    fullPath = [path fileName '_alignedSpikesPerNPerT_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _alignedSpikesPerNPerT_Run file does not exist');
        return;
    end
    load(fullPath,'trialsRunSpikes','trialsRewSpikes','trialsCueSpikes');
    
    s = num2str(costq);
    ind = findstr(s,'.');
    s(ind) = 'p';
    fileNameCorr = [fileName '_popSimiVP_Run' num2str(onlyRun) '_q' ...
            s '_intT' num2str(intervalT) '.mat'];
    
    fullPath = [path fileName '_behPar.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    GlobalConst;
    intervalT1 = intervalT*sampleFq;
    
    neuronNo = size(trialsRunSpikes.Time,1);
    nTrials = size(trialsRunSpikes.Time,2);
    maxSize = 100000000;
    minFR = 0.2;
    
    thrSimVR = 5;
    [indSelSimVP,meanSimVP] = selNeuronsSimVP(path,fileName,onlyRun,intervalT,costq,thrSimVR,minFR);
    
%     fullPath = [path fileName '_convSpikesAligned_Run' num2str(onlyRun) '.mat'];
%     if(exist(fullPath) == 0)
%         disp('The _convSpikesAligned_Run file does not exist');
%         return;
%     end
%     load(fullPath,'filteredSpikeArrayRun');
%     ind = find(indSelSimVP);
%     for i = 1:sum(indSelSimVP)
%         imagesc(filteredSpikeArrayRun{ind(i)});
%         title(['meanSimVP = ' num2str(meanSimVP(ind(i))) ' neuronNo = ' num2str(ind(i))]);
%         pause;
%     end
    
    disp('Convert population spike trains into a spike time vector for run onset aligned trials')
    [spikeTimeArrPerTr,spikeLabelArrPerTr] = convSpikeTrainToArr(indSelSimVP,nTrials,trialsRunSpikes.Time,intervalT1);
    disp('Calculate multi-neuron Victor & Purpura spike time distance for trials aligned to run onset')
    popSimVPRun = calSpikeTrainSimVP(nTrials,spikeTimeArrPerTr,spikeLabelArrPerTr,costq,costk,maxSize); 
    
    disp('Convert population spike trains into a spike time vector for reward onset aligned trials')
    [spikeTimeArrPerTr,spikeLabelArrPerTr] = convSpikeTrainToArr(indSelSimVP,nTrials,trialsRewSpikes.Time,intervalT1);
    disp('Calculate multi-neuron Victor & Purpura spike time distance for trials aligned to reward onset')
    popSimVPRew = calSpikeTrainSimVP(nTrials,spikeTimeArrPerTr,spikeLabelArrPerTr,costq,costk,maxSize);
    
    disp('Convert population spike trains into a spike time vector for cue onset aligned trials')
    [spikeTimeArrPerTr,spikeLabelArrPerTr] = convSpikeTrainToArr(indSelSimVP,nTrials,trialsCueSpikes.Time,intervalT1);
    disp('Calculate multi-neuron Victor & Purpura spike time distance for trials aligned to cue onset')
    popSimVPCue = calSpikeTrainSimVP(nTrials,spikeTimeArrPerTr,spikeLabelArrPerTr,costq,costk,maxSize);
    
    save([path fileNameCorr], ...
        'popSimVPRun','popSimVPRew','popSimVPCue',...
        'costq','costk','intervalT');
    
end

function [spikeTimeArrPerTr,spikeLabelArrPerTr] = convSpikeTrainToArr(indSelSimVP,nTrials,spikeTrain,intervalT)
    neuronNo = sum(indSelSimVP);
    indNeurons = find(indSelSimVP == 1);
    spikeTimeArrPerTr = cell(1,nTrials);
    spikeLabelArrPerTr = cell(1,nTrials);
    for j = 1:nTrials
        spTime = [];
        spLabel = [];
        for i = 1:neuronNo   
            neu = indNeurons(i);
            ind = spikeTrain{neu,j} < intervalT;
            nSp = sum(ind);
            if(nSp > 0)
                spTime = [spTime spikeTrain{neu,j}(ind)'];
                spLabel = [spLabel i*ones(1,nSp)];
            end
        end
        spikeTimeArrPerTr{j} = spTime;
        spikeLabelArrPerTr{j} = spLabel;
    end
end

function spikeTrainSimVP = calSpikeTrainSimVP(nTrials,spikeTimeArrPerTr,spikeLabelArrPerTr,costq,costk,maxSize)
    spikeTrainSimVP = zeros(nTrials,nTrials);
    for n = 1:nTrials
        if(isempty(spikeTimeArrPerTr{n}))
            continue;
        end
        spikeTrainX = spikeTimeArrPerTr{n};
        spikeTrainLX = spikeLabelArrPerTr{n};
        for m = n+1:nTrials
            if(isempty(spikeTimeArrPerTr{m}))
                continue;
            end
            spikeTrainY = spikeTimeArrPerTr{m};
            spikeTrainLY = spikeLabelArrPerTr{m};
            spikeTrainSimVP(n,m) = ...
                labdist_faster_qkpara_opt(spikeTrainX,spikeTrainLX,...
                    spikeTrainY,spikeTrainLY,costq,costk,maxSize);
        end
    end
    disp(['Max similarity = ' num2str(max(abs(spikeTrainSimVP(:))))]);
end

function [indSelCorrT,meanSimVP] = selNeuronsSimVP(path,fileName,onlyRun,intervalT,cost,thr,minFR)
% select neurons with high similarity based on single neuron correlation

    indSelCorrT = [];
    s = num2str(cost);
    ind = findstr(s,'.');
    s(ind) = 'p';
    fileNameCorr = [fileName '_meanSpikeTrainSimVP_Run' num2str(onlyRun) '_q' ...
            s '_intT' num2str(intervalT) '.mat'];
    fullPath = [path fileNameCorr];
    if(exist(fullPath) == 0)
        disp('The _meanSpikeTrainSimVR_Run file does not exist');
        return;
    end    
    load(fullPath,'meanSimVPRun');
    
    fileNameInfo = [fileName '_Info.mat'];
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        disp('_Info File does not exist.');
        return;
    end
    load(fullPath,'autoCorr');
    
    fileNameFR = [fileName '_FR_Run' num2str(onlyRun) '.mat'];
    fullPath = [path fileNameFR];
    if(exist(fullPath) == 0)
        disp('_FR_Run File does not exist.');
        return;
    end
    load(fullPath,'mFRStruct');
    
    indPyr = autoCorr.isPyrneuron;
    indCorrT = meanSimVPRun.meanGood < thr;
    indMFR = mFRStruct.mFR > minFR;
    indSelCorrT = indPyr & indCorrT & indMFR;
    meanSimVP = meanSimVPRun.meanGood;
    
end