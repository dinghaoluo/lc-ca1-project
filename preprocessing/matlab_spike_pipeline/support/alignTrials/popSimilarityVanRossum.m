function popSimilarityVanRossum(path,fileName,onlyRun,tc,intervalT)
% population level spike train similarity across trials, calculated using van
% Rossum distance
% e.g.: popSimilarityVanRossum('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,100)
    
    fullPath = [path fileName '_alignedSpikesPerNPerT_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _alignedSpikesPerNPerT_Run file does not exist');
        return;
    end
    load(fullPath,'spikeTrainRun','spikeTrainRew','spikeTrainCue');
    
    fileNameCorr = [fileName '_spikeTrainSimiVanR_Run' num2str(onlyRun) '_tc' ...
            num2str(tc) '_intT' num2str(intervalT) '.mat'];
    
    fullPath = [path fileName '_behPar.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    GlobalConst;
    intervalT1 = intervalT*sampleFq;
    neuronNo = size(spikeTrainRun{2},1);
    nTrials = length(spikeTrainRun);
    tcSample = tc/1000*sampleFq;
    
    disp('Convert spiek trains into spike arrays for run onset aligned trials')
    spikeArrPerNeuron = convSpikeTrainToArr(neuronNo,nTrials,spikeTrainRun,intervalT1);
    disp('Calculate population van Rossum distance for trials aligned to run onset')
    popSimVRRun = calPopSimVR(nTrials,spikeArrPerNeuron,tcSample); 
    
    disp('Convert spiek trains into spike arrays for reward onset aligned trials')
    spikeArrPerNeuron = convSpikeTrainToArr(neuronNo,nTrials,spikeTrainRew,intervalT1);
    disp('Calculate population van Rossum distance for trials aligned to reward onset')
    popSimVRRew = calPopSimVR(nTrials,spikeArrPerNeuron,tcSample);
    
    disp('Convert spiek trains into spike arrays for cue onset aligned trials')
    spikeArrPerNeuron = convSpikeTrainToArr(neuronNo,nTrials,spikeTrainCue,intervalT1);
    disp('Calculate population van Rossum distance for trials aligned to cue onset')
    popSimVRCue = calPopSimVR(nTrials,spikeArrPerNeuron,tcSample);
    
    save([path fileNameCorr], ...
        'popSimVRRun','popSimVRRew','spikeTrainSimVRCue',...
        'nonZeroTrRun','nonZeroTrRew','nonZeroTrCue','tc');
    
end

function spikeArrPerNeuron = convSpikeTrainToArr(neuronNo,nTrials,spikeTrain,intervalT)
    nonZeroTr = cell(1,neuronNo);
    nSpPerNeu = zeros(neuronNo,nTrials);
    maxNSpikes = zeros(1,neuronNo);
    for i = 1:neuronNo
        nonZeroTrTmp = zeros(1,nTrials);
        for j = 1:nTrials
            ind = spikeTrain{i,j} < intervalT;
            nSp = sum(ind);
            if(nSp > 0)
                nonZeroTrTmp(j) = 1;
                nSpPerNeu(i,j) = nSp;
            end
        end
        maxNSpikes(i) = max(nSpPerNeu(i,:));
        nonZeroTr{i} = nonZeroTrTmp;
    end
    
    spikeArrPerNeuron = cell(1,neuronNo);
    for i = 1:neuronNo
        spikeArrTmp = zeros(nTrials,maxNSpikes(i));
        for j = 1:nTrials
            if(nSpPerNeu(i,j) == 0)
                continue;
            end
            spikeArrTmp(j,1:nSpPerNeu(i,j)) = spikeTrain{i,j}(1:nSpPerNeu(i,j));
        end
        spikeArrPerNeuron{i} = spikeArrTmp;
    end
end

% function spikeArrPerNeuron = convSpikeTrainToArr(neuronNo,nTrials,spikeTrain,intervalT)
%     spikeArrPerNeuron = zeros(nTrials,neuronNo,intervalT);
%     for i = 1:nTrials     
%         if(isempty(spikeTrain{i}))
%             continue;
%         end
%         for j = 1:neuronNo
%             spikeTmp =  full(spikeTrain{i}(j,:));
%             if(sum(spikeTmp) == 0)
%                 continue;
%             end
%             if(length(spikeTmp) > intervalT)
%                 spikeArrPerNeuron(i,j,:) = spikeTmp(1:intervalT);
%             else
%                 spikeArrPerNeuron(i,j,1:length(spikeTmp)) = spikeTmp;
%             end
%         end
%     end
% end

function [popSimVR,nonZeroTr] = calPopSimVR(nTrials,spikes,tcSample)
    popSimVR = zeros(nTrials,nTrials);
    nonZeroTr = zeros(1,nTrials);
    for n = 1:nTrials
        popSpikeArrTmp = squeeze(spikes(n,:,:));
        nonZeroTrTmp = sum(sum(popSpikeArrTmp));
        nonZeroTr(n) = (nonZeroTrTmp > 0);
    end
    
    spikeTrainSimVRTmp = vanRossumMNPW(spikes(nonZeroTr==1,:,:),tcSample);
    popSimVR(nonZeroTr,nonZeroTr) = spikeTrainSimVRTmp;
end
