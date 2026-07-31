function RunSpeedOverDist2P(path, fileName, onlyRun, mazeSess)

    
    %%%%%%%%% initialize constants
    fileNameInfo = [fileName '_Info.mat'];     
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo_smTr2P(path,fileName);
    end
    load(fullPath);
    indTr = beh.indTrCtrl(beh.mazeSess(beh.indTrCtrl) == mazeSess);
    nTrials = length(indTr);
    
    fileName = [fileName '.mat'];
    fullPath = [path fileName];
    if(exist(fullPath) == 0)
        disp('Recording file does not exist.');
        return;
    end
    load(fullPath,'trials');
        
    GlobalConst2P;
    tracks = 1800;
    if(spaceMergeBin ~= 0)
        param.spaceSteps = [0:spaceMergeBin:tracks];
    else
        param.spaceSteps = [0:tracks];
    end
    
    numBins = length(param.spaceSteps);
    step = spaceMergeBin;
    timePerBin = zeros(nTrials,numBins);
    speedOverDist = zeros(nTrials,numBins);
    for tr = 1:nTrials
        indTrCur = indTr(tr);
        speedTmp = zeros(1,numBins);
        timePerBinTmp = zeros(1,numBins);
        if(onlyRun == 1)
            indSpeed = find(trials{indTrCur}.speedAll > minSpeed);
            distRun = trials{indTrCur}.xMM(indSpeed);
        else
            distRun = trials{indTrCur}.xMM;
            indSpeed = 1:length(trials{indTrCur}.xMM);
        end
        for i = 1:numBins
            ind = find(distRun >= param.spaceSteps(i)-step/2 & distRun < param.spaceSteps(i)+step/2);
            indOrig = indSpeed(ind);
            if(~isempty(ind))
                timePerBinTmp(i) = length(ind);
                speedTmp(i) = mean(trials{indTrCur}.speedAll(indOrig));
            else
                timePerBinTmp(i) = 1;
                if(i > 1)
                    speedTmp(i) = speedTmp(i-1);
                end
            end
        end
        speedOverDist(tr,:) = speedTmp;
        timePerBin(tr,:) = timePerBinTmp;
    end
    
    speedOverDistMean = mean(speedOverDist);
    speedOverDistStd = std(speedOverDist);
    speedOverDistSEM = std(speedOverDist)/sqrt(nTrials);
    
    save([path fileName(1:end-4) '_runSpeedDist_Run' num2str(onlyRun) '.mat'],'speedOverDist','timePerBin',...
        'speedOverDistMean','speedOverDistStd','speedOverDistSEM','param');
end