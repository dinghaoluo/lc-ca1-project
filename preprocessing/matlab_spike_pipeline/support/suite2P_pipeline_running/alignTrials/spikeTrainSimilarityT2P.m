function spikeTrainSimilarityT2P(path,fileName,onlyRun,mazeSess,intervalT)
% single neuron level spike correlation across trials
% e.g.: spikeTrainSimilarityT('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,10)
   
    fileNameCorr = [fileName '_spikeTrainSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    
    GlobalConst2P;
    
    disp('calculating single neuron cosine similarity aligned to run onset');
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
    neuronNo = length(filteredSpikeArrayRun); 
    
    [spikeTrainSimTRun,nonZeroTrRun] = calSpikeTrainSimT(neuronNo,filteredSpikeArrayRun,intervalT1); 
    clear filteredSpikeArrayRun
    
    disp('calculating single neuron cosine similarity aligned to reward onset');
    fullPath = [path fileName '_convSpikesAlignedRew_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _convSpikesAlignedRun_Run file does not exist');
        return;
    end
    load(fullPath,'filteredSpikeArrayRew');
    [spikeTrainSimTRew,nonZeroTrRew] = calSpikeTrainSimT(neuronNo,filteredSpikeArrayRew,intervalT1);
    clear filteredSpikeArrayRew
    
    disp('calculating single neuron cosine similarity aligned to cue onset');
    fullPath = [path fileName '_convSpikesAlignedCue_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _convSpikesAlignedRun_Run file does not exist');
        return;
    end
    load(fullPath,'filteredSpikeArrayCue');
    [spikeTrainSimTCue,nonZeroTrCue] = calSpikeTrainSimT(neuronNo,filteredSpikeArrayCue,intervalT1);
    clear filteredSpikeArrayCue
    
    save([path fileNameCorr], 'spikeTrainSimTRun','spikeTrainSimTRew','spikeTrainSimTCue',...
        'nonZeroTrRun','nonZeroTrRew','nonZeroTrCue','intervalT');
    
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
