function spikeTrainSimilarityVanRossum(path,fileName,onlyRun,mazeSess,tc,intervalT,intervalTMin)
% single neuron level spike correlation across trials
% van Rossum distance
% e.g.: spikeTrainSimilarityVanRossum('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,100,10)
% tc: time constant in ms
% intervalT: spike train max time in sec
    
    fullPath = [path fileName '_alignedSpikesPerNPerT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _alignedSpikesPerNPerT_Run file does not exist');
        return;
    end
    load(fullPath,'trialsRunSpikes','trialsRewSpikes','trialsCueSpikes');
    
    fileNameCorr = [fileName '_spikeTrainSimiVanR_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_tc' ...
            num2str(tc) '_intT' num2str(intervalT) '.mat'];
    
    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    GlobalConst;
    intervalT1 = intervalT*sampleFq;
    
    neuronNo = size(trialsRunSpikes.Time,1);
    nTrials = size(trialsRunSpikes.Time,2);
    tcSample = tc/1000*sampleFq;
        
    indMin = intervalTMin*sampleFq;
    indMax = intervalT*sampleFq;
    indInt = [indMin indMax];
    
    disp('Convert spike trains into spike arrays for run onset aligned trials')
    [spikeArrPerNeuron,nonZeroTrRun] = convSpikeTrainToArr(neuronNo,nTrials,trialsRunSpikes.Time,intervalT1);
    disp('Calculate van Rossum distance for trials aligned to run onset')
    spikeTrainSimVRRun = calSpikeTrainSimVR(neuronNo,nTrials,nonZeroTrRun,spikeArrPerNeuron,tcSample); 
    
    disp('Convert spike trains into spike arrays for reward onset aligned trials')
    [spikeArrPerNeuron,nonZeroTrRew] = convSpikeTrainToArr(neuronNo,nTrials,trialsRewSpikes.Time,intervalT1);
    disp('Calculate van Rossum distance for trials aligned to reward onset')
    spikeTrainSimVRRew = calSpikeTrainSimVR(neuronNo,nTrials,nonZeroTrRew,spikeArrPerNeuron,tcSample);
    
    disp('Convert spike trains into spike arrays for cue onset aligned trials')
    [spikeArrPerNeuron,nonZeroTrCue] = convSpikeTrainToArr(neuronNo,nTrials,trialsCueSpikes.Time,intervalT1);
    disp('Calculate van Rossum distance for trials aligned to cue onset')
    spikeTrainSimVRCue = calSpikeTrainSimVR(neuronNo,nTrials,nonZeroTrCue,spikeArrPerNeuron,tcSample);
    
    disp('Convert spike trains into spike arrays for run onset aligned trials (all spikes from last reward)')
    [spikeArrPerNeuron,nonZeroTrRun_LasttoCurTr] = convSpikeTrainToArr(neuronNo,nTrials,trialsRunSpikes.Time_LasttoCurTr,indInt);
    disp('Calculate van Rossum distance for trials aligned to run onset (all spikes from last reward)')
    spikeTrainSimVRRun_LasttoCurTr = calSpikeTrainSimVR(neuronNo,nTrials,nonZeroTrRun_LasttoCurTr,spikeArrPerNeuron,tcSample); 
    
    disp('Convert spike trains into spike arrays for reward onset aligned trials (all spikes from last reward)')
    [spikeArrPerNeuron,nonZeroTrRew_LasttoCurTr] = convSpikeTrainToArr(neuronNo,nTrials,trialsRewSpikes.Time_LasttoCurTr,indInt);
    disp('Calculate van Rossum distance for trials aligned to reward onset (all spikes from last reward)')
    spikeTrainSimVRRew_LasttoCurTr = calSpikeTrainSimVR(neuronNo,nTrials,nonZeroTrRew_LasttoCurTr,spikeArrPerNeuron,tcSample);
    
    disp('Convert spike trains into spike arrays for cue onset aligned trials (all spikes from last reward)')
    [spikeArrPerNeuron,nonZeroTrCue_LasttoCurTr] = convSpikeTrainToArr(neuronNo,nTrials,trialsCueSpikes.Time_LasttoCurTr,indInt);
    disp('Calculate van Rossum distance for trials aligned to cue onset (all spikes from last reward)')
    spikeTrainSimVRCue_LasttoCurTr = calSpikeTrainSimVR(neuronNo,nTrials,nonZeroTrCue_LasttoCurTr,spikeArrPerNeuron,tcSample);
    
    save([path fileNameCorr], ...
        'spikeTrainSimVRRun','spikeTrainSimVRRew','spikeTrainSimVRCue',...
        'nonZeroTrRun','nonZeroTrRew','nonZeroTrCue',...
        'spikeTrainSimVRRun_LasttoCurTr','spikeTrainSimVRRew_LasttoCurTr','spikeTrainSimVRCue_LasttoCurTr',...
        'nonZeroTrRun_LasttoCurTr','nonZeroTrRew_LasttoCurTr','nonZeroTrCue_LasttoCurTr','tc','intervalT');
    
end

function [spikeArrPerNeuron,nonZeroTr] = convSpikeTrainToArr(neuronNo,nTrials,spikeTrain,intervalT)
    nonZeroTr = cell(1,neuronNo);
    nSpPerNeu = zeros(neuronNo,nTrials);
    maxNSpikes = zeros(1,neuronNo);
    for i = 1:neuronNo
        nonZeroTrTmp = zeros(1,nTrials);
        for j = 1:nTrials
            if(i == 6 & j == 159)
                a = 1;
            end
            if(length(intervalT) == 1)
                ind = spikeTrain{i,j} < intervalT;
            else
                ind = spikeTrain{i,j} >= intervalT(1) ...
                    & spikeTrain{i,j} <= intervalT(2);
            end
            nSp = sum(ind);
            if(nSp == 1 & spikeTrain{i,j}(1) == 0)
                nSp = 0;
            end
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
            if(i == 6 & j == 159)
                a = 1;
            end
            if(nSpPerNeu(i,j) == 0)
                continue;
            end
            spikeArrTmp(j,1:nSpPerNeu(i,j)) = spikeTrain{i,j}(1:nSpPerNeu(i,j));
        end
        spikeArrPerNeuron{i} = spikeArrTmp;
    end
end

% function spikeArrPerNeuron = convSpikeTrainToArr(neuronNo,nTrials,spikeTrain,intervalT)
%     spikeArrPerNeuron = cell(1,neuronNo);
%     for i = 1:neuronNo
%         spikeArrPerNeuron{i} = zeros(nTrials,intervalT);
%     end
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
%                 spikeArrPerNeuron{j}(i,:) = spikeTmp(1:intervalT);
%             else
%                 spikeArrPerNeuron{j}(i,1:length(spikeTmp)) = spikeTmp;
%             end
%         end
%     end
% end

function spikeTrainSimVR = calSpikeTrainSimVR(neuronNo,nTrials,nonZeroTr,spikes,tcSample)
    spikeTrainSimVR = cell(1,neuronNo);
    for n = 1:neuronNo
        if(n == 6)
            a = 1;
        end
        disp(['Neuron no. ' num2str(n)]);
        spikeTrainSimVR{n} = zeros(nTrials,nTrials);

        if(sum(nonZeroTr{n}) > 1)
            ind = nonZeroTr{n}==1;
            spikeTrainSimVRTmp = vanRossumPW(spikes{n}(ind,:),tcSample);
            spikeTrainSimVR{n}(ind,ind) = spikeTrainSimVRTmp;
            disp(['Max similarity = ' num2str(max(abs(spikeTrainSimVR{n}(:))))]);
        end
    end
end

%         spikeTrainSimTmp = zeros(nTrials,nTrials);
%         nonZeroTrTmp = zeros(1,nTrials);
%         for i = 1:nTrials 
%             if(isempty(spikes{i}))
%                 continue;
%             end
%             spikesX = full(spikes{i}(n,:));
%             if(sum(spikesX) == 0)
%                 continue;
%             end
%             nonZeroTrTmp(i) = 1;
%             for j = i+1:nTrials
%                 spikesY = full(spikes{j}(n,:));
%                 if(sum(spikesY) == 0)
%                     continue;
%                 end
%                 spikeTrainSimTmp(i,j) = vanRossum(spikesX,spikesY,tcSample);
%             end
%         end