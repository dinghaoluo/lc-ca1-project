function spikeTrainSimilarityT(path,fileName,onlyRun,mazeSess,intervalT,intervalTMin)
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
    
    fileNameCorr = [fileName '_spikeTrainSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    
    GlobalConst;
    
    if(intervalT > paramC.trialLenT | intervalT == 0)
        intervalT = paramC.trialLenT;
    end
    intervalT1 = floor(intervalT/paramC.timeBin);
    neuronNo = length(filteredSpikeArrayRun);
    
    indMin = find(timeStep1 >= intervalTMin*sampleFq,1);
    indMax = find(timeStep1 <= intervalT*sampleFq,1,'last');
    indInt = [indMin indMax];
    
    disp('calculating single neuron cosine similarity aligned to run onset');
    [spikeTrainSimTRun,nonZeroTrRun] = calSpikeTrainSimT(neuronNo,filteredSpikeArrayRun,intervalT1); 
    disp('calculating single neuron cosine similarity aligned to reward onset');
    [spikeTrainSimTRew,nonZeroTrRew] = calSpikeTrainSimT(neuronNo,filteredSpikeArrayRew,intervalT1);
    disp('calculating single neuron cosine similarity aligned to cue onset');
    [spikeTrainSimTCue,nonZeroTrCue] = calSpikeTrainSimT(neuronNo,filteredSpikeArrayCue,intervalT1);
    
    disp('calculating single neuron cosine similarity aligned to run onset (all spikes from last reward)');
    [spikeTrainSimTRun_LasttoCurTr,nonZeroTrRun_LasttoCurTr] = calSpikeTrainSimT(neuronNo,filteredSpikeArrayRun_LasttoCurTr,indInt); 
    disp('calculating single neuron cosine similarity aligned to reward onset (all spikes from last reward)');
    [spikeTrainSimTRew_LasttoCurTr,nonZeroTrRew_LasttoCurTr] = calSpikeTrainSimT(neuronNo,filteredSpikeArrayRew_LasttoCurTr,indInt);
    disp('calculating single neuron cosine similarity aligned to cue onset (all spikes from last reward)');
    [spikeTrainSimTCue_LasttoCurTr,nonZeroTrCue_LasttoCurTr] = calSpikeTrainSimT(neuronNo,filteredSpikeArrayCue_LasttoCurTr,indInt);
    
    save([path fileNameCorr], 'spikeTrainSimTRun','spikeTrainSimTRew','spikeTrainSimTCue',...
        'nonZeroTrRun','nonZeroTrRew','nonZeroTrCue','intervalT',...
        'spikeTrainSimTRun_LasttoCurTr','spikeTrainSimTRew_LasttoCurTr','spikeTrainSimTCue_LasttoCurTr',...
        'nonZeroTrRun_LasttoCurTr','nonZeroTrRew_LasttoCurTr','nonZeroTrCue_LasttoCurTr','intervalT','indInt');
    
end

function [spikeTrainSimT,nonZeroTr] = calSpikeTrainSimT(neuronNo,filteredSpikeArray,intervalT)
    parfor n = 1:neuronNo
        spikeArr = filteredSpikeArray{n};
        if(length(intervalT) == 1)
            spikeArr = spikeArr(:,1:intervalT);
        else
            spikeArr = spikeArr(:,intervalT(1):intervalT(2));
        end
        nonZeroTr{n} = sum(spikeArr') > 0;
        spikeTrainSimTTmp = pdist(spikeArr,'cosine');
        spikeTrainSimTTmp = squareform(spikeTrainSimTTmp);
%         for i = 1:trialNo
%             if(sum(spikeArr(i,:)) ~= 0)
%                 for j = i+1:trialNo
%                     if(sum(spikeArr(j,:)) ~= 0)
%                         spikesCorrTmp(i,j) = ...
%                             corr(spikeArr(i,:)',spikeArr(j,:)','Type','Spearman');
%                     end
%                 end
%             end            
%         end
        spikeTrainSimT{n} = spikeTrainSimTTmp;
    end
end
