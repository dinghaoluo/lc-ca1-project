function neuronSpikeCorrDist_CtrlOnly(path,fileName,onlyRun,spaceBin,mazeSess)
% single neuron level spike correlation over distance across trials
% e.g.: neuronSpikeCorrDist_CtrlOnly('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,20)

    if(nargin == 4)
        spaceBin = 20;
    end
    
    fullPath = [path fileName '_convSpikesDist' num2str(spaceBin) ...
                    'mm_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _convSpikesDist file does not exist');
        return;
    end
    load(fullPath,'filteredSpikeArrayNormT','paramC');
    
    fileNameCorr = [fileName '_spikesCorrDist_Ctrl_Run' num2str(onlyRun) ...
        '_mazeSess' num2str(mazeSess) '.mat'];
    
    fileNameFR = [fileName '_FR_Ctrl_Run' num2str(onlyRun) '_mazeSess' ...
        num2str(mazeSess) '.mat'];
    fullPath = [path fileNameFR];
    if(exist(fullPath,'file') == 0)
        disp('_FR_Ctrl file does not exist. Please run "MeanFiringRateCtrlOnly" function first');
    end
    load(fullPath,'mFRStructSessCtrl');
    
    GlobalConst;
    fileNameInfo = [fileName '_Info.mat'];
        
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo_smTr(path,fileName);
    end
    load(fullPath);
    
    if(~isempty(filteredSpikeArrayNormT{2}))
        neuronNo = size(filteredSpikeArrayNormT{2},1);
        interval = size(filteredSpikeArrayNormT{2},2);
    else
        disp('filteredSpikeArrayNormT{2} is empty. Please check.');
        return;
    end
    
    disp('Calculate spike correlation over distance for ctrl trials')
    %%% calculate peak firing rate for good, bad and ok trials
    %%% within each session
    %%% added by Yingxue on 12/16/2020
    indLaps = mFRStructSessCtrl.indLapList; 
    if(~isempty(indLaps))
        [spikeCorrDistCtrl,nonZeroCtrl] = calSpikeCorr(neuronNo,filteredSpikeArrayNormT,indLaps,interval);
    end
        
    save([path fileNameCorr], 'spikeCorrDistCtrl','nonZeroCtrl','indLaps');
    
end

function [spikeCorr,nonZeroTr] = calSpikeCorr(neuronNo,filteredSpikeArray,indLaps,interval)
    nonZeroTr = cell(length(neuronNo),1);
    spikeCorr = cell(length(neuronNo),1);
    for n = 1:neuronNo
        spikeArr = zeros(length(indLaps),interval);
        m = 1;
        for i = 1:length(indLaps)        
            spikeArr(m,:) =  filteredSpikeArray{indLaps(i)}(n,1:interval);
            m = m+1;
        end
        nonZeroTr{n} = sum(spikeArr') > 0;
        spikesCorrTmp = corr(spikeArr','Type','Spearman');
        spikeCorr{n} = spikesCorrTmp;
    end
end