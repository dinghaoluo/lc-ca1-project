function thetaPhaseOverTimeAlignedCtrlOnly(path, fileName, mazeSess, onlyRun)
% convolve spike train with gaussian filter in time
% e.g. ConvSpikeTrain_Aligned('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',0)

    fullPath = [path fileName '_alignRun_msess' num2str(mazeSess) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _alignedRun file does not exist');
        return;
    end
    load(fullPath,'trialsRun');
    
    fullPath = [path fileName '_alignCue_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to cue file does not exist');
        return;
    end
    load(fullPath,'trialsCue');    
    
    fileNameThetaPhCtrl = [fileName '_thetaPhaseOverTimelignedCtrl_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    
    fullPath = [path fileName '_thetaPhaseOverTimeligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _thetaPhaseOverTimeligned file does not exist');
        return;
    end
    load(fullPath,'paramC');
    
    fullPath = [path fileName '_PeakFRAlignedCtrl_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];    
    if(exist(fullPath) == 0)
        disp('The _PeakFRAlignedCtrl file does not exist');
        return;
    end
    load(fullPath,'trialNoNonStim');
        
    GlobalConst;
    
    runSpeedNonStim = ...
        speedAlignedRun(trialsRun,trialNoNonStim,paramC.trialLenT*sampleFq);
        
    [thetaPhaseRunNonStimL,thetaPhaseRunNonStimH] = ...
        ThetaPhaseAlignedRun(trialsRun,trialNoNonStim,paramC.trialLenT*sampleFq);
    
    [thetaPhaseCueNonStimL,thetaPhaseCueNonStimH] = ...
        ThetaPhaseAlignedCue(trialsCue,trialNoNonStim,paramC.trialLenT*sampleFq-nSampBef);
   
    fullPath = [path fileNameThetaPhCtrl];
    save(fullPath, 'runSpeedNonStim',...
        'thetaPhaseRunNonStimL','thetaPhaseRunNonStimH',...
        'thetaPhaseCueNonStimL','thetaPhaseCueNonStimH',...
        'paramC','-v7.3'); 
    
end

function runSpeed = speedAlignedRun(trialsRun,trialNo,trialLen)
    runSpeed = zeros(length(trialNo),trialLen);
    for i = 1:length(trialNo)
        tmpSpeed = [trialsRun.speed_MMsecBef{trialNo(i)}',...
            trialsRun.speed_MMsec{trialNo(i)}'];
        tmpSpeed(tmpSpeed < 0) = 0;
        if(length(tmpSpeed) > trialLen)
            runSpeed(i,:) = tmpSpeed(1:trialLen);
        else
            runSpeed(i,1:length(tmpSpeed)) = tmpSpeed;
        end
    end
end

function runAcc = accAlignedRun(trialsRun,trialNo,trialLen)
    runAcc = zeros(length(trialNo),trialLen);
    for i = 1:length(trialNo)
        tmpAcc = trialsRun.accel_MMsecSq{trialNo(i)}';
        tmpSpeed = trialsRun.speed_MMsec{trialNo(i)}';
        tmpAcc(tmpSpeed < 0) = 0;
        if(length(tmpAcc) > trialLen)
            runAcc(i,:) = tmpAcc(1:trialLen);
        else
            runAcc(i,1:length(tmpAcc)) = tmpAcc;
        end
    end
end

function [thetaPhaseL,thetaPhaseH] = ThetaPhaseAlignedRun(trialsRun,trialNo,trialLen)
    thetaPhaseL = zeros(length(trialNo),trialLen);
    thetaPhaseH = zeros(length(trialNo),trialLen);
    for i = 1:length(trialNo)
        tmpThetaPhase = [trialsRun.thetaPhLinInterpBef{trialNo(i)}',...
            trialsRun.thetaPhLinInterp{trialNo(i)}'];
        tmpThetaPhaseH = [trialsRun.thetaPhHilbBef{trialNo(i)}',...
            trialsRun.thetaPhHilb{trialNo(i)}'];
        if(length(tmpThetaPhase) > trialLen)
            thetaPhaseL(i,:) = tmpThetaPhase(1:trialLen);
            thetaPhaseH(i,:) = tmpThetaPhaseH(1:trialLen);
        else
            thetaPhaseL(i,1:length(tmpThetaPhase)) = tmpThetaPhase;
            thetaPhaseH(i,1:length(tmpThetaPhase)) = tmpThetaPhaseH;
        end
    end
end

function [thetaPhaseL,thetaPhaseH] = ThetaPhaseAlignedCue(trialsCue,trialNo,trialLen)
    thetaPhaseL = zeros(length(trialNo),trialLen);
    thetaPhaseH = zeros(length(trialNo),trialLen);
    for i = 1:length(trialNo)
        tmpThetaPhase = trialsCue.thetaPhLinInterp{trialNo(i)}';
        tmpThetaPhaseH = trialsCue.thetaPhHilb{trialNo(i)}';
        if(length(tmpThetaPhase) > trialLen)
            thetaPhaseL(i,:) = tmpThetaPhase(1:trialLen);
            thetaPhaseH(i,:) = tmpThetaPhaseH(1:trialLen);
        else
            thetaPhaseL(i,1:length(tmpThetaPhase)) = tmpThetaPhase;
            thetaPhaseH(i,1:length(tmpThetaPhase)) = tmpThetaPhaseH;
        end
    end
end