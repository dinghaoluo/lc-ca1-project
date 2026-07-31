function LickOverDistBeh(path, fileName, mazeSess)

    fullPath = [path fileName '_alignRun_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to run file does not exist');
        disp(fullPath);
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
    
    fullPath = [path fileName '_behPar.mat'];
    if(exist(fullPath) == 0)
        disp('The behPar file does not exist');
        return;
    end
    load(fullPath,'param');
        
    GlobalConst;
    tracks = 2200;
    spaceMergeBin = 10; %mm
    if(spaceMergeBin ~= 0)
        param.spaceSteps = [-spaceMergeBin/2:spaceMergeBin:tracks+spaceMergeBin/2];
    else
        param.spaceSteps = [0:tracks];
    end
    
    nTrials = length(trialsRun.goodTrial);
    
    trInd = param.startTr:param.endTr;
    nTr = length(trInd);

    numBins = length(param.spaceSteps);
    step = spaceMergeBin;
    lickOverDist.Run = zeros(nTrials,numBins);
    lickOverDist.Rew = zeros(nTrials,numBins);
    lickOverDist.Cue = zeros(nTrials,numBins);
    for tr = 1:nTrials
        %% aligned to run onset
        licks = trialsRun.lickLfpInd{tr};
        if(~isempty(licks))
            licks = trialsRun.xMM{tr}(licks);
        else
            continue;
        end
        lickTmp = hist(licks,param.spaceSteps);
        lickOverDist.Run(tr,:) = lickTmp;
        
        %% aligned to reward
        licks = trialsRew.lickLfpInd{tr};
        if(~isempty(licks))
            licks = trialsRew.xMM{tr}(licks);
        else
            continue;
        end
        lickTmp = hist(licks,param.spaceSteps);
        lickOverDist.Rew(tr,:) = lickTmp;
        
        %% aligned to cue
        licks = trialsCue.lickLfpInd{tr};
        if(~isempty(licks))
            licks = trialsCue.xMM{tr}(licks);
        else
            continue;
        end
        lickTmp = hist(licks,param.spaceSteps);
        lickOverDist.Cue(tr,:) = lickTmp;
    end
    
    %%
    lickOverDist.meanRun = mean(lickOverDist.Run(trInd,:));
    lickOverDist.stdRun = std(lickOverDist.Run(trInd,:));
    lickOverDist.SEMRun = std(lickOverDist.Run(trInd,:))/sqrt(nTr);
        
    lickOverDist.meanRew = mean(lickOverDist.Rew(trInd,:));
    lickOverDist.stdRew = std(lickOverDist.Rew(trInd,:));
    lickOverDist.SEMRew = std(lickOverDist.Rew(trInd,:))/sqrt(nTr);
    
    lickOverDist.meanCue = mean(lickOverDist.Cue(trInd,:));
    lickOverDist.stdCue = std(lickOverDist.Cue(trInd,:));
    lickOverDist.SEMCue = std(lickOverDist.Cue(trInd,:))/sqrt(nTr);
    
    save([path fileName '_lickDist.mat'],'lickOverDist','param');
end