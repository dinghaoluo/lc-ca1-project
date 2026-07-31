function meanNeuronSpikeCorrDist_smTr2P(path,fileName,onlyRun,mazeSess,intervalD)
% single neuron level mean spike correlation across trials
% meanNeuronSpikeCorrDist('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,500)

    if(nargin == 4)
        intervalD = 0;
    end
    fullPath = [path fileName '_spikesCorrDist_Run' num2str(onlyRun) '_intD' ...
            num2str(intervalD) '.mat'];
    if(exist(fullPath) == 0)
        disp('The neuronSpikeCorrDist_smTr2P file does not exist');
        return;
    end
    load(fullPath);
        
    fileNameInfo = [fileName '_Info.mat'];
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        disp('The BasicInfo_smTr2P file does not exist');
        return;
    end
    load(fullPath);
    
    fileNameCorr = [fileName '_meanSpikesCorrDist_Run' num2str(onlyRun) '_intD' ...
            num2str(intervalD) '.mat'];
        
    GlobalConst2P;
    
    neuronNo = length(spikeCorrDist);
    trialNo = size(spikeCorrDist{1},1);
    indTr = 1:trialNo;
    indSelTr = indTr > startTrNo & beh.mazeSess' == mazeSess;
    trialNo = sum(indSelTr);
    nElem = (trialNo*trialNo-trialNo)/2;
    
    for n = 1:neuronNo
        corrArr = triu(spikeCorrDist{n}(indSelTr,indSelTr),1);
        corrArr = corrArr(:); 
        meanCorrDist.mean(n) = sum(corrArr(isnan(corrArr) == 0))/nElem;
        nNonZeroTr = sum(nonZeroTr{n} == 1 & indSelTr);
        nElemNonZero = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
        meanCorrDist.meanNZ(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZero;
        meanCorrDist.nNonZeroTr(n) = nNonZeroTr;
        
    end
    
    save([path fileNameCorr],'meanCorrDist');
end
