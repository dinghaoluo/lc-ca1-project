function indSelCorrT = selNeuronsSimT(path,fileName,onlyRun,intervalT,thr,minFR)
% select neurons with high similarity based on single neuron correlation

    indSelCorrT = [];
    fileNameCorr = [fileName '_meanSpikeTrainSimT_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    fullPath = [path fileNameCorr];
    if(exist(fullPath) == 0)
        disp('The _meanSpikeTrainSimT_Run file does not exist');
        return;
    end    
    load(fullPath,'meanSimTRun');
    
    fileNameInfo = [fileName '_Info.mat'];
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        disp('_Info File does not exist.');
        return;
    end
    load(fullPath,'autoCorr');
    
    fileNameFR = [fileName '_FR_Run' num2str(onlyRun) '.mat'];
    fullPath = [path fileNameFR];
    if(exist(fullPath) == 0)
        disp('_FR_Run File does not exist.');
        return;
    end
    load(fullPath,'mFRStruct');
    
    indPyr = autoCorr.isPyrneuron;
    indCorrT = meanSimTRun.meanGood < thr;
    indMFR = mFRStruct.mFR > minFR;
    indSelCorrT = indPyr & indCorrT & indMFR;
    
end
