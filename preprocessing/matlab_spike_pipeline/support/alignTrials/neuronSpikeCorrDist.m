function neuronSpikeCorrDist(path,fileName,onlyRun,mazeSess,intervalD)
% single neuron level spike correlation over distance across trials
% e.g.: neuronSpikeCorrDist('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,500)

    if(nargin == 4)
        intervalD = 0;
    end
    
    fullPath = [path fileName '_convSpikesDistAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _convSpikesDistAligned_Run file does not exist');
        return;
    end
    load(fullPath,'filteredSpikeDistNormTArrayRun','filteredSpikeDistNormTArrayRew',...
        'filteredSpikeDistNormTArrayCue','paramC');
    
    fileNameCorr = [fileName '_spikesCorrDistAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intD' ...
            num2str(intervalD) '.mat'];
    
    GlobalConst;
    
    if(intervalD > paramC.trackLen | intervalD == 0)
        intervalD = paramC.trackLen;
    end
    intervalD1 = floor(intervalD/spaceMergeBin);
    neuronNo = length(filteredSpikeDistNormTArrayRun);
    trialNo = size(filteredSpikeDistNormTArrayRun{1},1);
    
    [spikeCorrDistRun,nonZeroTrRun] = calSpikeCorr(neuronNo,filteredSpikeDistNormTArrayRun,intervalD1); 
    [spikeCorrDistRew,nonZeroTrRew] = calSpikeCorr(neuronNo,filteredSpikeDistNormTArrayRew,intervalD1);
    [spikeCorrDistCue,nonZeroTrCue] = calSpikeCorr(neuronNo,filteredSpikeDistNormTArrayCue,intervalD1);
    
    save([path fileNameCorr], 'spikeCorrDistRun','spikeCorrDistRew','spikeCorrDistCue',...
        'nonZeroTrRun','nonZeroTrRew','nonZeroTrCue','intervalD');
    
end

function [spikeCorr,nonZeroTr] = calSpikeCorr(neuronNo,filteredSpikeArray,intervalD)
    parfor n = 1:neuronNo
        spikeArr = filteredSpikeArray{n}';
        spikeArr = spikeArr(1:intervalD,:);
        nonZeroTr{n} = sum(spikeArr) > 0;
        spikesCorrTmp = corr(spikeArr,'Type','Spearman');
        
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
        spikeCorr{n} = spikesCorrTmp;
    end
end