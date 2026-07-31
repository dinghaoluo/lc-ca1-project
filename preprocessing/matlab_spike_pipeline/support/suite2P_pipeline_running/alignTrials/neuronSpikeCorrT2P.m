function neuronSpikeCorrT2P(path,fileName,onlyRun,mazeSess,intervalT)
% single neuron level spike correlation across trials
% e.g.: neuronSpikeCorrT('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,20)
    
    if(nargin == 4)
        intervalT = 0;
    end
    
    fileNameCorr = [fileName '_spikesCorrTAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    
    GlobalConst2P;
    
    disp('calculating single neuron correlation aligned to run onset');
    fullPath = [path fileName '_convSpikesAlignedRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _convSpikesAlignedRun_Run file does not exist');
        return;
    end
    load(fullPath,'filteredSpikeArrayRun','paramC');
    if(intervalT > paramC.trialLenT | intervalT == 0)
        intervalT = paramC.trialLenT;
    end
    intervalT1 = floor(intervalT/paramC.timeBin);
    neuronNo = length(filteredSpikeArrayRun);  
    
    [spikeCorrTRun,nonZeroTrRun] = calSpikeCorrT(neuronNo,filteredSpikeArrayRun,intervalT1); 
    clear filteredSpikeArrayRun
    
    disp('calculating single neuron correlation aligned to reward onset');
    fullPath = [path fileName '_convSpikesAlignedRew_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _convSpikesAlignedRun_Run file does not exist');
        return;
    end
    load(fullPath,'filteredSpikeArrayRew');
    [spikeCorrTRew,nonZeroTrRew] = calSpikeCorrT(neuronNo,filteredSpikeArrayRew,intervalT1);
    clear filteredSpikeArrayRew
    
    disp('calculating single neuron correlation aligned to cue onset');
    fullPath = [path fileName '_convSpikesAlignedCue_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _convSpikesAlignedRun_Run file does not exist');
        return;
    end
    load(fullPath,'filteredSpikeArrayCue');
    [spikeCorrTCue,nonZeroTrCue] = calSpikeCorrT(neuronNo,filteredSpikeArrayCue,intervalT1);
    clear filteredSpikeArrayCue
    
    save([path fileNameCorr], 'spikeCorrTRun','spikeCorrTRew','spikeCorrTCue',...
        'nonZeroTrRun','nonZeroTrRew','nonZeroTrCue',...
        'intervalT');    

end

function [spikeCorrT,nonZeroTr] = calSpikeCorrT(neuronNo,filteredSpikeArray,intervalT)
    parfor n = 1:neuronNo
        spikeArr = filteredSpikeArray{n}';
        if(length(intervalT) == 1)
            spikeArr = spikeArr(1:intervalT,:);
        else
            spikeArr = spikeArr(intervalT(1):intervalT(2),:);
        end
        nonZeroTr{n} = sum(spikeArr) > 0;
        spikeCorrT{n} = corr(spikeArr,'Type','Spearman');
    end
end