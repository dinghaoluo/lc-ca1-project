function RunSpeedAlignedNonRun(path, fileName, onlyRun, mazeSess)

    fullPath = [path fileName '_alignRun_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to run file does not exist');
        return;
    end
    load(fullPath,'trialsRun');
        
    fullPath = [path fileName '_PeakFRAligned_msess' num2str(mazeSess) '_Run1.mat'];    
    if(exist(fullPath) == 0)
        disp('The _PeakFRAligned file does not exist');
        return;
    end
    load(fullPath,'trialNoNonStimGood','trialNoNonStimBad','trialNoStim');
            
    GlobalConst;
    tracks = 1800;
    trialLenT = 20; % sec
    trialLenTBef = 3; % sec
    
    nTrials = length(trialsRun.goodTrial);
    
    trInd = 2:nTrials;
    nTr = nTrials-1;
    speedOverTime.Run = zeros(nTrials,(trialLenTBef+trialLenT)*sampleFq);
    
    numSampleAll = (trialLenTBef+trialLenT) * sampleFq;
    for tr = 2:nTrials
        %% aligned to run onset
        speed = [trialsRun.speed_MMsecBef{tr}' trialsRun.speed_MMsec{tr}'];
        speed(speed < 0) = 0;
        trLenCur = length(speed);
        if(trLenCur > numSampleAll)
            speedOverTime.Run(tr,:) = speed(1:numSampleAll);
        else
            speedOverTime.Run(tr,1:trLenCur) = speed;
        end
        
    end
    
    speedOverTime.meanRun = mean(speedOverTime.Run(trInd,:));
    speedOverTime.stdRun = std(speedOverTime.Run(trInd,:));
    speedOverTime.SEMRun = std(speedOverTime.Run(trInd,:))/sqrt(nTr);
    
    trialNoNonStimGood = trialNoNonStimGood(trialNoNonStimGood ~= 1);
    speedOverTime.meanRunNonStimGood = mean(speedOverTime.Run(trialNoNonStimGood,:));
    speedOverTime.stdRunNonStimGood = std(speedOverTime.Run(trialNoNonStimGood,:));
    speedOverTime.SEMRunNonStimGood = std(speedOverTime.Run(trialNoNonStimGood,:))...
        /sqrt(length(trialNoNonStimGood));
    
    trialNoNonStimBad = trialNoNonStimBad(trialNoNonStimBad ~= 1);
    speedOverTime.meanRunNonStimBad = mean(speedOverTime.Run(trialNoNonStimBad,:));
    speedOverTime.stdRunNonStimBad = std(speedOverTime.Run(trialNoNonStimBad,:));
    speedOverTime.SEMRunNonStimBad = std(speedOverTime.Run(trialNoNonStimBad,:))...
        /sqrt(length(trialNoNonStimBad));
    
    trialNoStim = trialNoStim(trialNoStim ~= 1);
    speedOverTime.meanRunStim = mean(speedOverTime.Run(trialNoStim,:));
    speedOverTime.stdRunStim = std(speedOverTime.Run(trialNoStim,:));
    speedOverTime.SEMRunStim = std(speedOverTime.Run(trialNoStim,:))...
        /sqrt(length(trialNoStim));
    
    speedOverTime.timeStep = (-trialLenTBef*sampleFq+1:trialLenT*sampleFq)/sampleFq;
    
    save([path fileName '_runSpeedDist_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'],...
        'speedOverTime');
end