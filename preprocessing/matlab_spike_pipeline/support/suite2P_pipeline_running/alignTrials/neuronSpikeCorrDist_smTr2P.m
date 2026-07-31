function neuronSpikeCorrDist_smTr2P(path,fileName,onlyRun,mazeSess,spaceBin,intervalD)
% single neuron level spike correlation over distance across trials
% e.g.: neuronSpikeCorrDist('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,500)

    if(nargin == 4)
        intervalD = 0;
    end
    
    fileNameConv = [fileName '_convSpikesDist' num2str(spaceBin) ...
                    'mm_Run' num2str(onlyRun) '.mat'];
    fullPath = [path fileNameConv];
    if(exist(fullPath) == 0)
        disp('The ConvSpikeTrainDistPar_smTr2P file does not exist');
        return;
    end
    load(fullPath,'filteredSpikeArrayNormT','paramC');
    
    fileNameCorr = [fileName '_spikesCorrDist_Run' num2str(onlyRun) '_intD' ...
            num2str(intervalD) '.mat'];
    
    GlobalConst2P;
    
    if(intervalD > max(paramC.spaceSteps{mazeSess}) | intervalD == 0)
        intervalD = max(paramC.spaceSteps{mazeSess});
    end
    intervalD1 = floor(intervalD/spaceMergeBin);
    neuronNo = size(filteredSpikeArrayNormT{end},1);
    trialNo = length(filteredSpikeArrayNormT);
    
    [spikeCorrDist,nonZeroTr] = calSpikeCorr(neuronNo,trialNo,filteredSpikeArrayNormT,intervalD1); 
    
    save([path fileNameCorr], 'spikeCorrDist','nonZeroTr','intervalD');
    
end

function [spikeCorr,nonZeroTr] = calSpikeCorr(neuronNo,trialNo,filteredSpikeArray,intervalD)
    lenTr = size(filteredSpikeArray{end},2);
    parfor n = 1:neuronNo
        spikeArr = zeros(trialNo,lenTr);
        for i = 1:trialNo
            if(~isempty(filteredSpikeArray{i})) % added by Yingxue on 02/05/2022, avoid trials longer than 100 s
                spikeArr(i,:) = filteredSpikeArray{i}(n,:);
            end
        end
    
        spikeArr = spikeArr';
        spikeArr = spikeArr(1:intervalD,:);
        nonZeroTr{n} = sum(spikeArr) > 0;
        spikeCorr{n} = corr(spikeArr,'Type','Spearman');
    end
end