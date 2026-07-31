function popSimilarityT(path,fileName,onlyRun,mazeSess,intervalT,intervalTMin)
% single neuron level spike correlation across trials
% e.g.: spikeTrainSimilarityT('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,10)
    
    fullPath = [path fileName '_convSpikesAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _convSpikesAligned_Run file does not exist');
        return;
    end
    load(fullPath,'filteredSpikeArrayRun','filteredSpikeArrayRew','filteredSpikeArrayCue',...
        'filteredSpikeArrayRun_LasttoCurTr','filteredSpikeArrayRew_LasttoCurTr','filteredSpikeArrayCue_LasttoCurTr',...
            'paramC','timeStep1');
    
    fileNameCorr = [fileName '_popSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    
    GlobalConst;
    
    thrSimT = 0.3;
    
    if(intervalT > paramC.trialLenT | intervalT == 0)
        intervalT = paramC.trialLenT;
    end
    intervalT1 = floor(intervalT/paramC.timeBin);
    indSelCorrT = selNeuronsSimT(path,fileName,onlyRun,intervalT,thrSimT,minFR);
%     ind = find(indSelCorrT);
%     for i = 1:sum(indSelCorrT)
%         imagesc(filteredSpikeArrayRun{ind(i)});
%         pause;
%     end
    
    nTrials = size(filteredSpikeArrayRun{1},1);
    
    indMin = find(timeStep1 >= intervalTMin*sampleFq,1);
    indMax = find(timeStep1 <= intervalT*sampleFq,1,'last');
    if(isempty(indMin))
        indMin = 1;
    end
    indInt = [indMin indMax];
    
    disp('Getting the population spike array per trial aligned to run onset');
    popSpikePerTr = convSpikeTrain(indSelCorrT,nTrials,filteredSpikeArrayRun,intervalT1);
    disp('calculating population cosine similarity aligned to run onset');
    popSimTRun = squareform(pdist(popSpikePerTr,'cosine'));
    
    disp('Getting the population spike array per trial aligned to reward onset');
    popSpikePerTr = convSpikeTrain(indSelCorrT,nTrials,filteredSpikeArrayRew,intervalT1);
    disp('calculating population cosine similarity aligned to run onset');
    popSimTRew = squareform(pdist(popSpikePerTr,'cosine'));
    
    disp('Getting the population spike array per trial aligned to cue onset');
    popSpikePerTr = convSpikeTrain(indSelCorrT,nTrials,filteredSpikeArrayCue,intervalT1);
    disp('calculating population cosine similarity aligned to cue onset');
    popSimTCue = squareform(pdist(popSpikePerTr,'cosine'));
    
    disp('Getting the population spike array per trial aligned to run onset (all spikes from last reward)');
    popSpikePerTr = convSpikeTrain1(indSelCorrT,nTrials,filteredSpikeArrayRun_LasttoCurTr,indInt);
    disp('calculating population cosine similarity aligned to run onset (all spikes from last reward)');
    popSimTRun_LasttoCurTr = squareform(pdist(popSpikePerTr,'cosine'));
    
    disp('Getting the population spike array per trial aligned to reward onset (all spikes from last reward)');
    popSpikePerTr = convSpikeTrain1(indSelCorrT,nTrials,filteredSpikeArrayRew_LasttoCurTr,indInt);
    disp('calculating population cosine similarity aligned to run onset (all spikes from last reward)');
    popSimTRew_LasttoCurTr = squareform(pdist(popSpikePerTr,'cosine'));
    
    disp('Getting the population spike array per trial aligned to cue onset (all spikes from last reward)');
    popSpikePerTr = convSpikeTrain1(indSelCorrT,nTrials,filteredSpikeArrayCue_LasttoCurTr,indInt);
    disp('calculating population cosine similarity aligned to cue onset (all spikes from last reward)');
    popSimTCue_LasttoCurTr = squareform(pdist(popSpikePerTr,'cosine'));
        
    save([path fileNameCorr], 'popSimTRun','popSimTRew','popSimTCue',...
        'popSimTRun_LasttoCurTr','popSimTRew_LasttoCurTr','popSimTCue_LasttoCurTr','intervalT','indInt');
    
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
%                 maxSp = max(filteredSpikeArrayRun{neu}(i,1:intervalT));
%                 if(maxSp > 0)
%                     popSpikePerTRTmp(j,:) = filteredSpikeArrayRun{neu}(i,1:intervalT)/maxSp;
%                 else
                    popSpikePerTRTmp(j,:) = filteredSpikeArrayRun{neu}(i,1:intervalT);
%                 end
            else
%                 maxSp = max(filteredSpikeArrayRun{neu}(i,:));
%                 if(maxSp > 0)
%                     popSpikePerTRTmp(j,1:totSample) = filteredSpikeArrayRun{neu}(i,:)/maxSp;
%                 else
                    popSpikePerTRTmp(j,1:totSample) = filteredSpikeArrayRun{neu}(i,:);
%                 end
            end
        end
        popSpikePerTr(i,:) = popSpikePerTRTmp(:)';
    end
end

function popSpikePerTr = convSpikeTrain1(indSelCorrT,nTrials,filteredSpikeArrayRun,intervalT)
    interval = (intervalT(2)-intervalT(1)+1);
    neuronNo = sum(indSelCorrT);
    indNeurons = find(indSelCorrT == 1);
    popSpikePerTr = zeros(nTrials,neuronNo*interval);
    for i = 1:nTrials
        popSpikePerTRTmp = zeros(neuronNo,interval);
        for j = 1:neuronNo
            neu = indNeurons(j);
            totSample = size(filteredSpikeArrayRun{neu},2);
            if(intervalT(2) <= totSample)
                maxSp = max(filteredSpikeArrayRun{neu}(i,intervalT(1):intervalT(2)));
                if(maxSp > 0)
                    popSpikePerTRTmp(j,:) = filteredSpikeArrayRun{neu}(i,intervalT(1):intervalT(2))/maxSp;
                else
                    popSpikePerTRTmp(j,:) = filteredSpikeArrayRun{neu}(i,intervalT(1):intervalT(2));
                end
            else
                maxSp = max(filteredSpikeArrayRun{neu}(i,:));
                if(maxSp > 0)
                    popSpikePerTRTmp(j,intervalT(1):totSample) = filteredSpikeArrayRun{neu}(i,:)/maxSp;
                else
                    popSpikePerTRTmp(j,intervalT(1):totSample) = filteredSpikeArrayRun{neu}(i,:);
                end
            end
        end
        popSpikePerTr(i,:) = popSpikePerTRTmp(:)';
    end
end
