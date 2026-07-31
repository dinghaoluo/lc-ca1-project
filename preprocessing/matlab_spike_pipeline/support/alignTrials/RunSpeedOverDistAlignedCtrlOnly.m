function RunSpeedOverDistAlignedCtrlOnly(path, fileName, onlyRun, mazeSess)

    fullPath = [path fileName '_runSpeedDist_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _runSpeedDist file does not exist');
        return;
    end
    load(fullPath,'speedOverDist','timePerBin','param');
    
    fullPath = [path fileName '_PeakFRAlignedCtrl_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];    
    if(exist(fullPath) == 0)
        disp('The _PeakFRAlignedCtrl file does not exist');
        return;
    end
    load(fullPath,'trialNoNonStim');
            
    GlobalConst;
    
    speedOverDistCtrl.meanRunNonStim = mean(speedOverDist.Run(trialNoNonStim,:));
    speedOverDistCtrl.stdRunNonStim = std(speedOverDist.Run(trialNoNonStim,:));
    speedOverDistCtrl.SEMRunNonStim = std(speedOverDist.Run(trialNoNonStim,:))...
        /sqrt(length(trialNoNonStim));
    
    speedOverDistCtrl.meanRewNonStim = mean(speedOverDist.Rew(trialNoNonStim,:));
    speedOverDistCtrl.stdRewNonStim = std(speedOverDist.Rew(trialNoNonStim,:));
    speedOverDistCtrl.SEMRewNonStim = std(speedOverDist.Rew(trialNoNonStim,:))...
        /sqrt(length(trialNoNonStim));
    
    speedOverDistCtrl.meanCueNonStim = mean(speedOverDist.Cue(trialNoNonStim,:));
    speedOverDistCtrl.stdCueNonStim = std(speedOverDist.Cue(trialNoNonStim,:));
    speedOverDistCtrl.SEMCueNonStim = std(speedOverDist.Cue(trialNoNonStim,:))...
        /sqrt(length(trialNoNonStim));
    
    save([path fileName '_runSpeedDistCtrl_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'],...
        'speedOverDistCtrl','timePerBin','param');
end