function LickOverDistAlignedCtrlOnly(path, fileName, onlyRun, mazeSess)

    fullPath = [path fileName '_lickDist_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _lickDist file does not exist');
        return;
    end
    load(fullPath,'lickOverDist','param');
    
    fullPath = [path fileName '_PeakFRAlignedCtrl_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];    
    if(exist(fullPath) == 0)
        disp('The _PeakFRAlignedCtrl file does not exist');
        return;
    end
    load(fullPath,'trialNoNonStim');
        
    GlobalConst;
    
    lickOverDistCtrl.meanRunNonStim = mean(lickOverDist.Run(trialNoNonStim,:));
    lickOverDistCtrl.stdRunNonStim = std(lickOverDist.Run(trialNoNonStim,:));
    lickOverDistCtrl.SEMRunNonStim = std(lickOverDist.Run(trialNoNonStim,:))...
        /sqrt(length(trialNoNonStim)); 

    lickOverDistCtrl.meanRewNoStim = mean(lickOverDist.Rew(trialNoNonStim,:));
    lickOverDistCtrl.stdRewNoStim = std(lickOverDist.Rew(trialNoNonStim,:));
    lickOverDistCtrl.SEMRewNoStim = std(lickOverDist.Rew(trialNoNonStim,:))...
        /sqrt(length(trialNoNonStim));

    lickOverDistCtrl.meanCueStim = mean(lickOverDist.Cue(trialNoNonStim,:));
    lickOverDistCtrl.stdCueStim = std(lickOverDist.Cue(trialNoNonStim,:));
    lickOverDistCtrl.SEMCueStim = std(lickOverDist.Cue(trialNoNonStim,:))...
        /sqrt(length(trialNoNonStim));
    
    save([path fileName '_lickDistCtrl_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'],'lickOverDistCtrl','param');
end