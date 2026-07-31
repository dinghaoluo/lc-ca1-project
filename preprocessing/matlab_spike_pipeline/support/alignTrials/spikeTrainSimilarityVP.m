function spikeTrainSimilarityVP(path,fileName,onlyRun,mazeSess,cost,intervalT,intervalTMin)
% single neuron level spike correlation across trials
% Victor & Purpura spike time distance and Victor & Purpura spike time
% interval
% e.g.: spikeTrainSimilarityVP('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,0.016,8)
% cost: cost per unit time to move a spike
% intervalT: spike train max time in sec
    
    fullPath = [path fileName '_alignedSpikesPerNPerT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _alignedSpikesPerNPerT_Run file does not exist');
        return;
    end
    load(fullPath,'trialsRunSpikes','trialsRewSpikes','trialsCueSpikes');
       
    s = num2str(cost);
    ind = findstr(s,'.');
    s(ind) = 'p';
    fileNameCorr = [fileName '_spikeTrainSimVP_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_q' ...
            s '_intT' num2str(intervalT) '.mat'];
    
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
    
    indMin = intervalTMin*sampleFq;
    indMax = intervalT*sampleFq;
    indInt = [indMin indMax];
    
    disp('Convert spike trains into spike arrays for run onset aligned trials')
    [spikeArrPerNeuron,nonZeroTrRun,nSpPerNeuRun] = convSpikeTrainToArr(neuronNo,nTrials,trialsRunSpikes.Time,intervalT1);
    disp('Calculate Victor & Purpura spike time distance for trials aligned to run onset')
    spikeTrainSimVPRun = calSpikeTrainSimVP(neuronNo,nTrials,nonZeroTrRun,spikeArrPerNeuron,nSpPerNeuRun,cost); 
    disp('Calculate Victor & Purpura spike time interval for trials aligned to run onset')
    spikeTrainSimVPIRun = calSpikeTrainSimVPI(neuronNo,nTrials,nonZeroTrRun,spikeArrPerNeuron,nSpPerNeuRun,cost,intervalT1); 
    
    disp('Convert spike trains into spike arrays for reward onset aligned trials')
    [spikeArrPerNeuron,nonZeroTrRew,nSpPerNeuRew] = convSpikeTrainToArr(neuronNo,nTrials,trialsRewSpikes.Time,intervalT1);
    disp('Calculate Victor & Purpura spike time distance for trials aligned to reward onset')
    spikeTrainSimVPRew = calSpikeTrainSimVP(neuronNo,nTrials,nonZeroTrRew,spikeArrPerNeuron,nSpPerNeuRew,cost);
    disp('Calculate Victor & Purpura spike time interval for trials aligned to reward onset')
    spikeTrainSimVPIRew = calSpikeTrainSimVPI(neuronNo,nTrials,nonZeroTrRew,spikeArrPerNeuron,nSpPerNeuRew,cost,intervalT1);
    
    disp('Convert spike trains into spike arrays for cue onset aligned trials')
    [spikeArrPerNeuron,nonZeroTrCue,nSpPerNeuCue] = convSpikeTrainToArr(neuronNo,nTrials,trialsCueSpikes.Time,intervalT1);
    disp('Calculate Victor & Purpura spike time distance for trials aligned to cue onset')
    spikeTrainSimVPCue = calSpikeTrainSimVP(neuronNo,nTrials,nonZeroTrCue,spikeArrPerNeuron,nSpPerNeuCue,cost);
    disp('Calculate Victor & Purpura spike time interval for trials aligned to cue onset')
    spikeTrainSimVPICue = calSpikeTrainSimVPI(neuronNo,nTrials,nonZeroTrCue,spikeArrPerNeuron,nSpPerNeuCue,cost,intervalT1);
    
    disp('Convert spike trains into spike arrays for run onset aligned trials (all spikes from last reward)')
    [spikeArrPerNeuron,nonZeroTrRun_LasttoCurTr,nSpPerNeuRun_LasttoCurTr] = convSpikeTrainToArr(neuronNo,nTrials,trialsRunSpikes.Time_LasttoCurTr,indInt);
    disp('Calculate Victor & Purpura spike time distance for trials aligned to run onset (all spikes from last reward)')
    spikeTrainSimVPRun_LasttoCurTr = calSpikeTrainSimVP(neuronNo,nTrials,nonZeroTrRun_LasttoCurTr,spikeArrPerNeuron,nSpPerNeuRun_LasttoCurTr,cost); 
    disp('Calculate Victor & Purpura spike time interval for trials aligned to run onset (all spikes from last reward)')
    spikeTrainSimVPIRun_LasttoCurTr = calSpikeTrainSimVPI(neuronNo,nTrials,nonZeroTrRun_LasttoCurTr,spikeArrPerNeuron,nSpPerNeuRun_LasttoCurTr,cost,indInt); 
    
    disp('Convert spike trains into spike arrays for reward onset aligned trials (all spikes from last reward)')
    [spikeArrPerNeuron,nonZeroTrRew_LasttoCurTr,nSpPerNeuRew_LasttoCurTr] = convSpikeTrainToArr(neuronNo,nTrials,trialsRewSpikes.Time_LasttoCurTr,indInt);
    disp('Calculate Victor & Purpura spike time distance for trials aligned to reward onset (all spikes from last reward)')
    spikeTrainSimVPRew_LasttoCurTr = calSpikeTrainSimVP(neuronNo,nTrials,nonZeroTrRew_LasttoCurTr,spikeArrPerNeuron,nSpPerNeuRew_LasttoCurTr,cost);
    disp('Calculate Victor & Purpura spike time interval for trials aligned to reward onset (all spikes from last reward)')
    spikeTrainSimVPIRew_LasttoCurTr = calSpikeTrainSimVPI(neuronNo,nTrials,nonZeroTrRew_LasttoCurTr,spikeArrPerNeuron,nSpPerNeuRew_LasttoCurTr,cost,indInt);
    
    disp('Convert spike trains into spike arrays for cue onset aligned trials (all spikes from last reward)')
    [spikeArrPerNeuron,nonZeroTrCue_LasttoCurTr,nSpPerNeuCue_LasttoCurTr] = convSpikeTrainToArr(neuronNo,nTrials,trialsCueSpikes.Time_LasttoCurTr,indInt);
    disp('Calculate Victor & Purpura spike time distance for trials aligned to cue onset (all spikes from last reward)')
    spikeTrainSimVPCue_LasttoCurTr = calSpikeTrainSimVP(neuronNo,nTrials,nonZeroTrCue_LasttoCurTr,spikeArrPerNeuron,nSpPerNeuCue_LasttoCurTr,cost);
    disp('Calculate Victor & Purpura spike time interval for trials aligned to cue onset (all spikes from last reward)')
    spikeTrainSimVPICue_LasttoCurTr = calSpikeTrainSimVPI(neuronNo,nTrials,nonZeroTrCue_LasttoCurTr,spikeArrPerNeuron,nSpPerNeuCue_LasttoCurTr,cost,indInt);
    
    save([path fileNameCorr], ...
        'spikeTrainSimVPRun','spikeTrainSimVPRew','spikeTrainSimVPCue',...
        'spikeTrainSimVPIRun','spikeTrainSimVPIRew','spikeTrainSimVPICue',...
        'nonZeroTrRun','nonZeroTrRew','nonZeroTrCue',...
        'nSpPerNeuRun','nSpPerNeuRew','nSpPerNeuCue',...
        'spikeTrainSimVPRun_LasttoCurTr','spikeTrainSimVPRew_LasttoCurTr','spikeTrainSimVPCue_LasttoCurTr',...
        'spikeTrainSimVPIRun_LasttoCurTr','spikeTrainSimVPIRew_LasttoCurTr','spikeTrainSimVPICue_LasttoCurTr',...
        'nonZeroTrRun_LasttoCurTr','nonZeroTrRew_LasttoCurTr','nonZeroTrCue_LasttoCurTr',...
        'nSpPerNeuRun_LasttoCurTr','nSpPerNeuRew_LasttoCurTr','nSpPerNeuCue_LasttoCurTr',...
        'cost','intervalT','indInt');
    
end

function [spikeArrPerNeuron,nonZeroTr,nSpPerNeu] = convSpikeTrainToArr(neuronNo,nTrials,spikeTrain,intervalT)
    nonZeroTr = cell(1,neuronNo);
    nSpPerNeu = zeros(neuronNo,nTrials);
    maxNSpikes = zeros(1,neuronNo);
    for i = 1:neuronNo
        nonZeroTrTmp = zeros(1,nTrials);
        for j = 1:nTrials
            if(length(intervalT) == 1)
                ind = spikeTrain{i,j} < intervalT;
            else
                ind = spikeTrain{i,j} >= intervalT(1) ...
                    & spikeTrain{i,j} <= intervalT(2);
            end
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

function spikeTrainSimVP = calSpikeTrainSimVP(neuronNo,nTrials,nonZeroTr,spikes,nSpPerNeu,cost)
    spikeTrainSimVP = cell(1,neuronNo);
    parfor n = 1:neuronNo
        disp(['Neuron no. ' num2str(n)]);
        spikeTrainSimVP{n} = zeros(nTrials,nTrials);

        if(sum(nonZeroTr{n}) > 2)
            ind = find(nonZeroTr{n}==1);
            for tr = 1:length(ind)
                spikeTrainX = spikes{n}(ind(tr),1:nSpPerNeu(n,ind(tr)));
                for tr1 = tr+1:length(ind)
                    spikeTrainY = spikes{n}(ind(tr1),1:nSpPerNeu(n,ind(tr1)));
                    spikeTrainSimVP{n}(ind(tr),ind(tr1)) = ...
                        spkd(spikeTrainX,spikeTrainY,cost);
                end
            end
            disp(['Max similarity = ' num2str(max(abs(spikeTrainSimVP{n}(:))))]);
        end
    end
end

function spikeTrainSimVPI = calSpikeTrainSimVPI(neuronNo,nTrials,nonZeroTr,spikes,nSpPerNeu,cost,intervalT)
    if(length(intervalT)>1)
        intervalT = intervalT(end) - intervalT(1);
    end
    spikeTrainSimVPI = cell(1,neuronNo);
    for n = 1:neuronNo
        disp(['Neuron no. ' num2str(n)]);
        spikeTrainSimVPI{n} = zeros(nTrials,nTrials);

        if(sum(nonZeroTr{n}) > 2)
            ind = find(nonZeroTr{n}==1);
            for tr = 1:length(ind)
                spikeTrainX = spikes{n}(ind(tr),1:nSpPerNeu(n,ind(tr)));
                for tr1 = tr+1:length(ind)
                    spikeTrainY = spikes{n}(ind(tr1),1:nSpPerNeu(n,ind(tr1)));
                    spikeTrainSimVPI{n}(ind(tr),ind(tr1)) = ...
                        spkd_int_FAST_post(spikeTrainX,spikeTrainY,cost,intervalT);
                end
            end
            disp(['Max similarity = ' num2str(max(abs(spikeTrainSimVPI{n}(:))))]);
        end
    end
end
