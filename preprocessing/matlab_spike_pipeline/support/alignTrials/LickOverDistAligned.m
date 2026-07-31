function LickOverDistAligned(path, fileName, onlyRun, mazeSess)

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
    
    fullPath = [path fileName '_PeakFRAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];    
    if(exist(fullPath) == 0)
        disp('The _PeakFRAligned file does not exist');
        return;
    end
    load(fullPath,'trialNoNonStimGood','trialNoNonStimBad','trialNoStim','trialNoStimCtrl','pulseMeth');
        
    GlobalConst;
    tracks = 2200;
    spaceMergeBin = 10; %mm
    if(spaceMergeBin ~= 0)
        param.spaceSteps = [-spaceMergeBin/2:spaceMergeBin:tracks+spaceMergeBin/2];
    else
        param.spaceSteps = [0:tracks];
    end
    
    nTrials = length(trialsRun.goodTrial);
    
    trInd = 1:nTrials;
    nTr = nTrials;
    
    numBins = length(param.spaceSteps);
    step = spaceMergeBin;
    lickOverDist.Run = zeros(nTrials,numBins);
    lickOverDist.Rew = zeros(nTrials,numBins);
    lickOverDist.Cue = zeros(nTrials,numBins);
    for tr = 1:nTrials
        %% aligned to run onset
        licks = trialsRun.lickLfpInd{tr} - trialsRun.startLfpInd(tr);
        if(~isempty(licks))
            licks = trialsRun.xMM{tr}(licks);
        else
            continue;
        end
        lickTmp = hist(licks,param.spaceSteps);
        lickOverDist.Run(tr,:) = lickTmp;
        
        %% aligned to reward
        licks = trialsRew.lickLfpInd{tr} - trialsRew.startLfpInd(tr);
        if(~isempty(licks))
            licks = trialsRew.xMM{tr}(licks);
        else
            continue;
        end
        lickTmp = hist(licks,param.spaceSteps);
        lickOverDist.Rew(tr,:) = lickTmp;
        
        %% aligned to cue
        licks = trialsCue.lickLfpInd{tr} - trialsCue.startLfpInd(tr);
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
    
    lickOverDist.meanRunNonStimGood = mean(lickOverDist.Run(trialNoNonStimGood,:));
    lickOverDist.stdRunNonStimGood = std(lickOverDist.Run(trialNoNonStimGood,:));
    lickOverDist.SEMRunNonStimGood = std(lickOverDist.Run(trialNoNonStimGood,:))...
        /sqrt(length(trialNoNonStimGood));
    
    lickOverDist.meanRunNonStimBad = mean(lickOverDist.Run(trialNoNonStimBad,:));
    lickOverDist.stdRunNonStimBad = std(lickOverDist.Run(trialNoNonStimBad,:));
    lickOverDist.SEMRunNonStimBad = std(lickOverDist.Run(trialNoNonStimBad,:))...
        /sqrt(length(trialNoNonStimBad));
    
    %% added by Yingxue on 2/15/2021
    lickOverDist.meanRunStim = [];
    lickOverDist.stdRunStim = [];
    lickOverDist.SEMRunStim = [];
    
    lickOverDist.meanRunStimCtrl = [];
    lickOverDist.stdRunStimCtrl = [];
    lickOverDist.SEMRunStimCtrl = [];
    for i = 1:length(pulseMeth)
        lickOverDist.meanRunStim{i} = mean(lickOverDist.Run(trialNoStim{i},:));
        lickOverDist.stdRunStim{i} = std(lickOverDist.Run(trialNoStim{i},:));
        lickOverDist.SEMRunStim{i} = std(lickOverDist.Run(trialNoStim{i},:))...
            /sqrt(length(trialNoStim{i}));

        %% added by Yingxue on 2/14/2021
        lickOverDist.meanRunStimCtrl{i} = mean(lickOverDist.Run(trialNoStimCtrl{i},:));
        lickOverDist.stdRunStimCtrl{i} = std(lickOverDist.Run(trialNoStimCtrl{i},:));
        lickOverDist.SEMRunStimCtrl{i} = std(lickOverDist.Run(trialNoStimCtrl{i},:))...
            /sqrt(length(trialNoStimCtrl{i}));
        %%
    end

    lickOverDist.meanRew = mean(lickOverDist.Rew(trInd,:));
    lickOverDist.stdRew = std(lickOverDist.Rew(trInd,:));
    lickOverDist.SEMRew = std(lickOverDist.Rew(trInd,:))/sqrt(nTr);

    lickOverDist.meanCue = mean(lickOverDist.Cue(trInd,:));
    lickOverDist.stdCue = std(lickOverDist.Cue(trInd,:));
    lickOverDist.SEMCue = std(lickOverDist.Cue(trInd,:))/sqrt(nTr);
    
    save([path fileName '_lickDist_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'],'lickOverDist','param');
end