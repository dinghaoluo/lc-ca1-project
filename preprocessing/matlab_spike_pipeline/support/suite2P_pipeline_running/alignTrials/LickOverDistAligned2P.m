function LickOverDistAligned2P(path, fileName, onlyRun, mazeSess)

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
    
    fullPath = [path fileName '_PeakFRAlignedRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];    
    if(exist(fullPath) == 0)
        disp('The _PeakFRAligned file does not exist');
        return;
    end
    load(fullPath,'trialNoNonStimGood','trialNoNonStimBad','trialNoStim','trialNoStimCtrl','pulseMeth');
        
    GlobalConst2P;
    tracks = 2200;
    spaceMergeBin = 10; %mm  
    if(spaceMergeBin ~= 0)
        param.spaceSteps = [-spaceMergeBin/2:spaceMergeBin:tracks+spaceMergeBin/2];
    else
        param.spaceSteps = [0:tracks];
    end
    
    nTrials = length(trialsRun.goodTrial);
    
    trInd = startTrNo+1:nTrials;
    nTr = length(trInd);
    
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
    
    trIndTmp = setdiff(trialNoNonStimGood,1:startTrNo);
    lickOverDist.meanRunNonStimGood = mean(lickOverDist.Run(trIndTmp,:));
    lickOverDist.stdRunNonStimGood = std(lickOverDist.Run(trIndTmp,:));
    lickOverDist.SEMRunNonStimGood = std(lickOverDist.Run(trIndTmp,:))...
        /sqrt(length(trIndTmp));
    
    trIndTmp = setdiff(trialNoNonStimBad,1:startTrNo);
    lickOverDist.meanRunNonStimBad = mean(lickOverDist.Run(trIndTmp,:));
    lickOverDist.stdRunNonStimBad = std(lickOverDist.Run(trIndTmp,:));
    lickOverDist.SEMRunNonStimBad = std(lickOverDist.Run(trIndTmp,:))...
        /sqrt(length(trIndTmp));
    
    %% added by Yingxue on 2/15/2021
    lickOverDist.meanRunStim = [];
    lickOverDist.stdRunStim = [];
    lickOverDist.SEMRunStim = [];
    
    lickOverDist.meanRunStimCtrl = [];
    lickOverDist.stdRunStimCtrl = [];
    lickOverDist.SEMRunStimCtrl = [];
    for i = 1:length(pulseMeth)
        trIndTmp = setdiff(trialNoStim{i},1:startTrNo);
        lickOverDist.meanRunStim{i} = mean(lickOverDist.Run(trIndTmp,:));
        lickOverDist.stdRunStim{i} = std(lickOverDist.Run(trIndTmp,:));
        lickOverDist.SEMRunStim{i} = std(lickOverDist.Run(trIndTmp,:))...
            /sqrt(length(trIndTmp));

        %% added by Yingxue on 2/14/2021
        trIndTmp = setdiff(trialNoStimCtrl{i},1:startTrNo);
        lickOverDist.meanRunStimCtrl{i} = mean(lickOverDist.Run(trIndTmp,:));
        lickOverDist.stdRunStimCtrl{i} = std(lickOverDist.Run(trIndTmp,:));
        lickOverDist.SEMRunStimCtrl{i} = std(lickOverDist.Run(trIndTmp,:))...
            /sqrt(length(trIndTmp));
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