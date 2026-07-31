function meanNeuronSpikeCorrDist_CtrlOnly(path,fileName,onlyRun,mazeSess)
% single neuron level mean spike correlation across trials
% meanNeuronSpikeCorrDist_CtrlOnly('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1)

    fullPath = [path fileName '_spikesCorrDist_Ctrl_Run' num2str(onlyRun) ...
        '_mazeSess' num2str(mazeSess) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _spikesCorrDist_Ctrl_Run file does not exist');
        return;
    end
    load(fullPath);
    
    fileNameFR = [fileName '_FR_Ctrl_Run' num2str(onlyRun) '_mazeSess' ...
        num2str(mazeSess) '.mat'];
    fullPath = [path fileNameFR];
    if(exist(fullPath,'file') == 0)
        disp('_FR_Ctrl file does not exist. Please run "MeanFiringRateCtrlOnly" function first');
    end
    load(fullPath,'mFRStructSessCtrl');
    indLaps = mFRStructSessCtrl.indLapList; 
                
    fileNameCorr = [fileName '_meanSpikesCorrDist_Ctrl_Run' num2str(onlyRun) ...
        '_mazeSess' num2str(mazeSess) '.mat'];
    
    fileNameInfo = [fileName '_Info.mat'];     
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo_smTr(path,fileName);
    end
    load(fullPath);
    
    disp('Calculate mean spike correlation over distance for ctrl trials')
    trialNoCtrl = length(indLaps);
    nElemCtrl = (trialNoCtrl*trialNoCtrl-trialNoCtrl)/2;
        
    if(trialNoCtrl > 0)
        neuronNo = size(spikeCorrDistCtrl,2);
        for n = 1:neuronNo            
            corrArr = triu(spikeCorrDistCtrl{n},1);
            corrArr = corrArr(:); 
            meanCorrDistCtrl.mean(n) = sum(corrArr(isnan(corrArr) == 0))/nElemCtrl;
            nNonZeroTr = sum(nonZeroCtrl{n} == 1);
            nElemNonZero = (nNonZeroTr*nNonZeroTr-nNonZeroTr)/2;
            meanCorrDistCtrl.meanNZ(n) = sum(corrArr(isnan(corrArr) == 0))/nElemNonZero;
            meanCorrDistCtrl.nNonZeroTr(n) = nNonZeroTr;
        end
    end
    
    save([path fileNameCorr],'meanCorrDistCtrl');
    
end
