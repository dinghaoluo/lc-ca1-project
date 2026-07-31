function indNeuSel = selPyrNeurons2P(path,fileName,onlyRun,minFR)
% select neurons with high similarity based on single neuron correlation

    indNeuSel = [];
    fileNameInfo = [fileName '_Info.mat'];
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        disp('_info File does not exist.');
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
    indMFR = mFRStruct.mFR > minFR;
    indNeuSel = indPyr & indMFR;
    
end
