function spikeThreshlding(path,fileName)
% threshold spikes and dFF based on noise level

    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(isempty(indexFileName))
        fileNameMinLen = [fileName '_Info.mat'];
        fileName = [fileName '.mat'];
    else
        fileNameMinLen = [fileName(1:indexFileName(end)-1) '_Info.mat'];
    end 
    fullPath = [path fileName];
    if(exist(fullPath) == 0)
        disp('File does not exist.');
        return;
    end
    load(fullPath,'trials');
    
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo_smTr2P(path,fileName);
    end
    load(fullPath);
    
    GlobalConst2P;
    
    numNeurons = rec.numNeurons;
    numTr = beh.numTrials;
    
    for i = 1:numTr
        % smooth and threshold dFF (Following David Tank, nature, 2021)
        param.stdsm = stdSM;
        h = gaussFilter(6*param.stdsm,param.stdsm);
        lenGaussKernel = length(h);
        dFFSM = zero(size(dFF,1),size(dFF,2));
        spikesSM = trials{i}.spikes;
        stdDFF = std(dFF,0,2);
        meanDFF = mean(dFF,2);
        param.threStd = threStd;
        for j = 1:length(Clu.localClu)
            dFFSM(i,:) = gauss_filter(dFF(i,:),h,lenGaussKernel);
            indNoise = dFFSM(i,:) <= meanDFF + param.threStd * stdDFF;
            dFFSM(i,indNoise) = 0;
            spikesSM(i,indNoise) = 0;
        end 
        trials{i}.dFFSM = dFFSM;
        trials{i}.spikesSM = spikesSM;
    end
    
    fullPath = [path fileName];
    save(fullpath,'trials','-append');
end
