function popSimilarityT2P(path,fileName,onlyRun,mazeSess,intervalT)
% single neuron level spike correlation across trials
% e.g.: spikeTrainSimilarityT('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,10)
    
    fileNameCorr = [fileName '_meanSpikeTrainSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    fullPath = [path fileNameCorr];
    if(exist(fullPath) == 0)
        disp('The _meanSpikesCorrTAligned_Run file does not exist');
        return;
    end    
    load(fullPath,'meanSimTRun');
    
    fileNameFR = [fileName '_FR_Run' num2str(onlyRun) '.mat'];
    fullPath = [path fileNameFR];
    if(exist(fullPath) == 0)
        disp('_FR_Run File does not exist.');
        return;
    end
    load(fullPath,'mFRStruct');
       
    fileNameCorr = [fileName '_popSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    
    GlobalConst2P;
    
    thrSimT = 0.3;
    
    disp('Getting the population spike array per trial aligned to run onset');
    fullPath = [path fileName '_convSpikesAlignedRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _convSpikesAligned_Run file does not exist');
        return;
    end
    load(fullPath,'filteredSpikeArrayRun','paramC');
    
    if(intervalT > paramC.trialLenT | intervalT == 0)
        intervalT = paramC.trialLenT;
    end
    intervalT1 = floor(intervalT/paramC.timeBin);
    indSelCorrT = meanSimTRun.meanGood < thrSimT & mFRStruct.mFR > minFR;

    nTrials = size(filteredSpikeArrayRun{1},1);
        
    popSpikePerTr = convSpikeTrain(indSelCorrT,nTrials,filteredSpikeArrayRun,intervalT1);
    disp('calculating population cosine similarity aligned to run onset');
    popSimTRun = squareform(pdist(popSpikePerTr,'cosine'));
    clear filteredSpikeArrayRun
    
    disp('Getting the population spike array per trial aligned to reward onset');
    fullPath = [path fileName '_convSpikesAlignedRew_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _convSpikesAlignedRun_Run file does not exist');
        return;
    end
    load(fullPath,'filteredSpikeArrayRew');
    popSpikePerTr = convSpikeTrain(indSelCorrT,nTrials,filteredSpikeArrayRew,intervalT1);
    disp('calculating population cosine similarity aligned to run onset');
    popSimTRew = squareform(pdist(popSpikePerTr,'cosine'));
    clear filteredSpikeArrayRew
    
    disp('Getting the population spike array per trial aligned to cue onset');
    fullPath = [path fileName '_convSpikesAlignedCue_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _convSpikesAlignedRun_Run file does not exist');
        return;
    end
    load(fullPath,'filteredSpikeArrayCue');
    popSpikePerTr = convSpikeTrain(indSelCorrT,nTrials,filteredSpikeArrayCue,intervalT1);
    disp('calculating population cosine similarity aligned to cue onset');
    popSimTCue = squareform(pdist(popSpikePerTr,'cosine'));
    clear filteredSpikeArrayCue
    
    save([path fileNameCorr], 'indSelCorrT','popSimTRun','popSimTRew','popSimTCue','intervalT');
    
end

function popSpikePerTr = convSpikeTrain(indSelCorrT,nTrials,filteredSpikeArrayRun,intervalT)
    neuronNo = sum(indSelCorrT);
    indNeurons = find(indSelCorrT == 1);
    popSpikePerTr = zeros(nTrials,neuronNo*intervalT);
    for i = 1:nTrials
        popSpikePerTRTmp = zeros(neuronNo,intervalT);
        for j = 1:neuronNo
            neu = indNeurons(j);
            totSample = size(filteredSpikeArrayRun{neu},2);
            if(intervalT <= totSample)
                maxSp = max(filteredSpikeArrayRun{neu}(i,1:intervalT));
                if(maxSp > 0)
                    popSpikePerTRTmp(j,:) = filteredSpikeArrayRun{neu}(i,1:intervalT)/maxSp;
                else
                    popSpikePerTRTmp(j,:) = filteredSpikeArrayRun{neu}(i,1:intervalT);
                end
            else
                maxSp = max(filteredSpikeArrayRun{neu}(i,:));
                if(maxSp > 0)
                    popSpikePerTRTmp(j,1:totSample) = filteredSpikeArrayRun{neu}(i,:)/maxSp;
                else
                    popSpikePerTRTmp(j,1:totSample) = filteredSpikeArrayRun{neu}(i,:);
                end
            end
        end
        popSpikePerTRTmp = popSpikePerTRTmp';
        popSpikePerTr(i,:) = popSpikePerTRTmp(:);
    end
end
