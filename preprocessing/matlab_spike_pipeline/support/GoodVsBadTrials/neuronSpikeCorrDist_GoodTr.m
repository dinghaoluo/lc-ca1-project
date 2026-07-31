function neuronSpikeCorrDist_GoodTr(path,fileName,onlyRun,spaceBin)
% single neuron level spike correlation over distance across trials
% e.g.: neuronSpikeCorrDist_GoodTr('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,20)

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
    
    fileNameCorr = [fileName '_spikesCorrDist_GoodTr_Run' num2str(onlyRun) '.mat'];
    
    GlobalConst;
    fileNameInfo = [fileName '_Info.mat'];
        
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo_smTr(path,fileName);
    end
    load(fullPath);
    mazeSess = beh.mazeSessAll;
    
    spikeCorrDistGoodTr = cell(length(mazeSess),1);
    spikeCorrDistBadTr = cell(length(mazeSess),1);
    spikeCorrDistOKTr = cell(length(mazeSess),1);
    
    nonZeroGoodTr = cell(length(mazeSess),1);
    nonZeroBadTr = cell(length(mazeSess),1);
    nonZeroOKTr = cell(length(mazeSess),1);
    
    if(~isempty(filteredSpikeArrayNormT{2}))
        neuronNo = size(filteredSpikeArrayNormT{2},1);
        interval = size(filteredSpikeArrayNormT{2},2);
    else
        disp('filteredSpikeArrayNormT{2} is empty. Please check.');
        return;
    end
    if(length(mazeSess) > 1)
        disp('Calculate spike correlation over distance for each session')
        for i = 1:length(mazeSess) 
            fprintf('\nSession %d\n',i);
            %%% calculate peak firing rate for good, bad and ok trials
            %%% within each session
            %%% added by Yingxue on 12/16/2020
            indLaps = find(beh.mazeSess == mazeSess(i));
            indLaps = intersect(indLaps,beh.indTrCtrl); 
            indLapsGoodTr{i} = intersect(indLaps,beh.indGoodTrCtrl); 
            if(~isempty(indLapsGoodTr{i}))
                [spikeCorrDistGoodTr{i},nonZeroGoodTr{i}] = calSpikeCorr(neuronNo,filteredSpikeArrayNormT,indLapsGoodTr{i},interval);
            end
            
            indLapsBadTr{i} = intersect(indLaps,beh.indBadTrCtrl); 
            if(~isempty(indLapsBadTr{i}))
                [spikeCorrDistBadTr{i},nonZeroBadTr{i}] = calSpikeCorr(neuronNo, filteredSpikeArrayNormT,indLapsBadTr{i},interval);
            end
            
            indLapsOKTr{i} = setdiff(indLaps,[beh.indGoodTrCtrl,beh.indBadTrCtrl]); 
            if(~isempty(indLapsOKTr{i}))
                [spikeCorrDistOKTr{i},nonZeroOKTr{i}] = calSpikeCorr(neuronNo, filteredSpikeArrayNormT,indLapsOKTr{i},interval);
            end
        end
    else
        disp('Calculate spike correlation over distance')
        indLapsGoodTr{1} = beh.indGoodTrCtrl; 
        if(~isempty(indLapsGoodTr{1}))
            [spikeCorrDistGoodTr{1},nonZeroGoodTr{1}] = calSpikeCorr(neuronNo,filteredSpikeArrayNormT,indLapsGoodTr{1},interval);
        end
        
        indLapsBadTr{1} = beh.indBadTrCtrl; 
        if(~isempty(indLapsBadTr{1}))
            [spikeCorrDistBadTr{1},nonZeroBadTr{1}] = calSpikeCorr(neuronNo,filteredSpikeArrayNormT,indLapsBadTr{1},interval);
        end
        
        indLapsOKTr{1} = setdiff(beh.indTrCtrl,[beh.indGoodTrCtrl,beh.indBadTrCtrl]);
        if(~isempty(indLapsOKTr{1}))
            [spikeCorrDistOKTr{1},nonZeroOKTr{1}] = calSpikeCorr(neuronNo,filteredSpikeArrayNormT,indLapsOKTr{1},interval);
        end
    end
        
    save([path fileNameCorr], 'spikeCorrDistGoodTr','spikeCorrDistBadTr','spikeCorrDistOKTr',...
        'nonZeroGoodTr','nonZeroBadTr','nonZeroOKTr','indLapsGoodTr','indLapsBadTr','indLapsOKTr');
    
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