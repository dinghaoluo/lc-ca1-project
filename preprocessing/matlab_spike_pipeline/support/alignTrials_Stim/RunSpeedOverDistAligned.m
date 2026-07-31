function RunSpeedOverDistAligned(path, fileName, onlyRun, mazeSess)

    fullPath = [path fileName '_alignRun_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to run file does not exist');
        return;
    end
    load(fullPath,'trialsRun');
    
    fullPath = [path fileName '_alignRew_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to reward file does not exist');
        return;
    end
    load(fullPath,'trialsRew');
    
    fullPath = [path fileName '_alignCue_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to cue file does not exist');
        return;
    end
    load(fullPath,'trialsCue');
    
    fullPath = [path fileName '_PeakFRAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];    
    if(exist(fullPath) == 0)
        disp('The _PeakFRAligned file does not exist');
        return;
    end
    load(fullPath,'trialNoNonStimGood','trialNoNonStimBad','trialNoStim');
            
    GlobalConst;
    tracks = 1800;
    if(spaceMergeBin ~= 0)
        param.spaceSteps = [0:spaceMergeBin:tracks];
    else
        param.spaceSteps = [0:tracks];
    end
    
    nTrials = length(trialsRun.goodTrial);
    
    trInd = 1:nTrials;
    nTr = nTrials;
    
    numBins = length(param.spaceSteps);
    step = spaceMergeBin;
    timePerBin.Run = zeros(nTrials,numBins);
    speedOverDist.Run = zeros(nTrials,numBins);
    timePerBin.Rew = zeros(nTrials,numBins);
    speedOverDist.Rew = zeros(nTrials,numBins);
    timePerBin.Cue = zeros(nTrials,numBins);
    speedOverDist.Cue = zeros(nTrials,numBins);
    for tr = 1:nTrials
        %% aligned to run onset
        speedTmp = zeros(1,numBins);
        timePerBinTmp = zeros(1,numBins);
        if(onlyRun == 1)
            indSpeed = find(trialsRun.speed_MMsec{tr} > minSpeed);
            distRun = trialsRun.xMM{tr}(indSpeed);
        else
            distRun = trialsRun.xMM{tr};
            indSpeed = 1:length(trialsRun.xMM{tr});
        end
        for i = 1:numBins
            ind = find(distRun >= param.spaceSteps(i)-step/2 & distRun < param.spaceSteps(i)+step/2);
            indOrig = indSpeed(ind);
            if(~isempty(ind))
                timePerBinTmp(i) = length(ind);
                speedTmp(i) = mean(trialsRun.speed_MMsec{tr}(indOrig));
            else
                timePerBinTmp(i) = 1;
                if(i > 1)
                    speedTmp(i) = speedTmp(i-1);
                end
            end
        end
        speedOverDist.Run(tr,:) = speedTmp;
        timePerBin.Run(tr,:) = timePerBinTmp;
        
        %% aligned to reward
        speedTmp = zeros(1,numBins);
        timePerBinTmp = zeros(1,numBins);
        if(onlyRun == 1)
            indSpeed = find(trialsRew.speed_MMsec{tr} > minSpeed);
            distRun = trialsRew.xMM{tr}(indSpeed);
        else
            distRun = trialsRew.xMM{tr};
            indSpeed = 1:length(trialsRew.xMM{tr});
        end
        for i = 1:numBins
            ind = find(distRun >= param.spaceSteps(i)-step/2 & distRun < param.spaceSteps(i)+step/2);
            indOrig = indSpeed(ind);
            if(~isempty(ind))
                timePerBinTmp(i) = length(ind);
                speedTmp(i) = mean(trialsRew.speed_MMsec{tr}(indOrig));
            else
                timePerBinTmp(i) = 1;
                if(i > 1)
                    speedTmp(i) = speedTmp(i-1);
                end
            end
        end
        speedOverDist.Rew(tr,:) = speedTmp;
        timePerBin.Rew(tr,:) = timePerBinTmp;
        
        %% aligned to reward
        speedTmp = zeros(1,numBins);
        timePerBinTmp = zeros(1,numBins);
        if(onlyRun == 1)
            indSpeed = find(trialsCue.speed_MMsec{tr} > minSpeed);
            distRun = trialsCue.xMM{tr}(indSpeed);
        else
            distRun = trialsCue.xMM{tr};
            indSpeed = 1:length(trialsCue.xMM{tr});
        end
        for i = 1:numBins
            ind = find(distRun >= param.spaceSteps(i)-step/2 & distRun < param.spaceSteps(i)+step/2);
            indOrig = indSpeed(ind);
            if(~isempty(ind))
                timePerBinTmp(i) = length(ind);
                speedTmp(i) = mean(trialsCue.speed_MMsec{tr}(indOrig));
            else
                timePerBinTmp(i) = 1;
                if(i > 1)
                    speedTmp(i) = speedTmp(i-1);
                end
            end
        end
        speedOverDist.Cue(tr,:) = speedTmp;
        timePerBin.Cue(tr,:) = timePerBinTmp;
    end
    
    speedOverDist.meanRun = mean(speedOverDist.Run(trInd,:));
    speedOverDist.stdRun = std(speedOverDist.Run(trInd,:));
    speedOverDist.SEMRun = std(speedOverDist.Run(trInd,:))/sqrt(nTr);
    
    speedOverDist.meanRunNonStimGood = mean(speedOverDist.Run(trialNoNonStimGood,:));
    speedOverDist.stdRunNonStimGood = std(speedOverDist.Run(trialNoNonStimGood,:));
    speedOverDist.SEMRunNonStimGood = std(speedOverDist.Run(trialNoNonStimGood,:))...
        /sqrt(length(trialNoNonStimGood));
    
    speedOverDist.meanRunNonStimBad = mean(speedOverDist.Run(trialNoNonStimBad,:));
    speedOverDist.stdRunNonStimBad = std(speedOverDist.Run(trialNoNonStimBad,:));
    speedOverDist.SEMRunNonStimBad = std(speedOverDist.Run(trialNoNonStimBad,:))...
        /sqrt(length(trialNoNonStimBad));

    speedOverDist.meanRunStim = mean(speedOverDist.Run(trialNoStim,:));
    speedOverDist.stdRunStim = std(speedOverDist.Run(trialNoStim,:));
    speedOverDist.SEMRunStim = std(speedOverDist.Run(trialNoStim,:))...
        /sqrt(length(trialNoStim));
    
    speedOverDist.meanRew = mean(speedOverDist.Rew(trInd,:));
    speedOverDist.stdRew = std(speedOverDist.Rew(trInd,:));
    speedOverDist.SEMRew = std(speedOverDist.Rew(trInd,:))/sqrt(nTr);
    
    speedOverDist.meanCue = mean(speedOverDist.Cue(trInd,:));
    speedOverDist.stdCue = std(speedOverDist.Cue(trInd,:));
    speedOverDist.SEMCue = std(speedOverDist.Cue(trInd,:))/sqrt(nTr);
    
    save([path fileName '_runSpeedDist_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'],...
        'speedOverDist','timePerBin','param');
end