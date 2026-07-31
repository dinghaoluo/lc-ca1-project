function thetaPhaseOverTimeAligned_Stim(path, fileName, mazeSess)
% convolve spike train with gaussian filter in time, stimulation session
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
    
    fileNameThetaPh = [fileName '_thetaPhaseOverTimeligned_msess' num2str(mazeSess) '.mat'];
    
    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    fullPath = [path fileName '_Info.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to run file does not exist');
        return;
    end
    load(fullPath,'beh');
    
    GlobalConst;
    if(isfield(beh,'indStimLap'))
        ind = beh.mazeSess == mazeSess;
        trialNoNonStimGood = find(beh.indStimLap(ind) == 0 & behPar.indTrBadBeh == 0);
        trialNoNonStimBad = find(beh.indStimLap(ind) == 0 & behPar.indTrBadBeh == 1);
        trialNoStim = find(beh.indStimLap(ind) == 1);
        %% added by Yingxue on 2/3/2021
        trialNoStimGood = find(beh.indStimLap(ind) == 1 & behPar.indTrBadBeh == 0);
        trialNoStimBad = find(beh.indStimLap(ind) == 1 & behPar.indTrBadBeh == 1);
    else
        trialNoNonStimGood = find(behPar.indTrBadBeh == 0);
        trialNoNonStimBad = find(behPar.indTrBadBeh == 1);
        trialNoStim = [];
        %% added by Yingxue on 2/3/2021
        trialNoStimGood = [];
        trialNoStimBad = [];
    end
    
    paramC.trialLenT = 20; %sec
    paramC.timeIncBef = ((0:trialLenT*sampleFq-1) - nSampBef)/sampleFq;
    paramC.time = (1:(trialLenT*sampleFq - nSampBef))/sampleFq;
    
    runSpeedNonStimGood = ...
        speedAlignedRun(trialsRun,trialNoNonStimGood,paramC.trialLenT*sampleFq);
    runSpeedNonStimBad = ...
        speedAlignedRun(trialsRun,trialNoNonStimBad,paramC.trialLenT*sampleFq);
    %% added by Yingxue on 2/3/2021
    runSpeedStimGood = ...
        speedAlignedRun(trialsRun,trialNoStimGood,paramC.trialLenT*sampleFq);
    runSpeedStimBad = ...
        speedAlignedRun(trialsRun,trialNoStimBad,paramC.trialLenT*sampleFq);
    %%
    
%     runAccNonStimGood = ...
%         accAlignedRun(trialsRun,trialNoNonStimGood,paramC.trialLenT*sampleFq);
    
    [thetaPhaseRunNonStimGoodL,thetaPhaseRunNonStimGoodH] = ...
        ThetaPhaseAlignedRun(trialsRun,trialNoNonStimGood,paramC.trialLenT*sampleFq);
    [thetaPhaseRunNonStimBadL,thetaPhaseRunNonStimBadH] = ...
        ThetaPhaseAlignedRun(trialsRun,trialNoNonStimBad,paramC.trialLenT*sampleFq);
    [thetaPhaseRunStimL,thetaPhaseRunStimH] = ...
        ThetaPhaseAlignedRun(trialsRun,trialNoStim,paramC.trialLenT*sampleFq);
    %% added by Yingxue on 2/3/2021
    [thetaPhaseRunStimGoodL,thetaPhaseRunStimGoodH] = ...
        ThetaPhaseAlignedRun(trialsRun,trialNoStimGood,paramC.trialLenT*sampleFq);
    [thetaPhaseRunStimBadL,thetaPhaseRunStimBadH] = ...
        ThetaPhaseAlignedRun(trialsRun,trialNoStimBad,paramC.trialLenT*sampleFq);
    %%
    
    [thetaPhaseCueNonStimGoodL,thetaPhaseCueNonStimGoodH] = ...
        ThetaPhaseAlignedCue(trialsCue,trialNoNonStimGood,paramC.trialLenT*sampleFq-nSampBef);
    [thetaPhaseCueNonStimBadL,thetaPhaseCueNonStimBadH] = ...
        ThetaPhaseAlignedCue(trialsCue,trialNoNonStimBad,paramC.trialLenT*sampleFq-nSampBef);
    [thetaPhaseCueStimL,thetaPhaseCueStimH] = ...
        ThetaPhaseAlignedCue(trialsCue,trialNoStim,paramC.trialLenT*sampleFq-nSampBef);
    %% added by Yingxue on 2/3/2021
    [thetaPhaseCueStimGoodL,thetaPhaseCueStimGoodH] = ...
        ThetaPhaseAlignedCue(trialsCue,trialNoStimGood,paramC.trialLenT*sampleFq-nSampBef);
    [thetaPhaseCueStimBadL,thetaPhaseCueStimBadH] = ...
        ThetaPhaseAlignedCue(trialsCue,trialNoStimBad,paramC.trialLenT*sampleFq-nSampBef);
    %%
    
    fullPath = [path fileNameThetaPh];
    save(fullPath, 'runSpeedNonStimGood','runSpeedNonStimBad',...
        'runSpeedStimGood','runSpeedStimBad',...
        'thetaPhaseRunNonStimGoodL','thetaPhaseRunNonStimGoodH',...
        'thetaPhaseRunNonStimBadL','thetaPhaseRunNonStimBadH',...
        'thetaPhaseRunStimL','thetaPhaseRunStimH',...
        'thetaPhaseRunStimGoodL','thetaPhaseRunStimGoodH',...
        'thetaPhaseRunStimBadL','thetaPhaseRunStimBadH',...
        'thetaPhaseCueNonStimGoodL','thetaPhaseCueNonStimGoodH',...
        'thetaPhaseCueNonStimBadL','thetaPhaseCueNonStimBadH',...
        'thetaPhaseCueStimL','thetaPhaseCueStimH',...
        'thetaPhaseCueStimGoodL','thetaPhaseCueStimGoodH',...
        'thetaPhaseCueStimBadL','thetaPhaseCueStimBadH','paramC','-v7.3'); 
    
%     [~,indSpeed] = sort(behPar.meanSpeed(trialNoNonStimGood));
% %     [~,indSpeed] = sort(behPar.maxSpeed(trialNoNonStimGood));
% %     [~,indSpeed] = sort(behPar.numSamplesRun(trialNoNonStimGood));
%     [figNew,pos] = CreateFig();
%     set(0,'Units','pixels') 
%     set(figure(figNew),'OuterPosition',[pos(1) pos(2) pos(3) pos(4)])
%     subplot(3,1,1)
%     imagesc(paramC.timeIncBef,1:length(trialNoNonStimGood),runSpeedNonStimGood(indSpeed,:));
%     set(gca,'XLim',[-0.5 3])
%     ylabel('Trial no. (RunOnset speed)')
%     subplot(3,1,2)
%     imagesc(paramC.timeIncBef,1:length(trialNoNonStimGood),thetaPhaseRunNonStimGoodL(indSpeed,:));
%     set(gca,'XLim',[-0.5 3])
%     ylabel('Trial no. (RunOnset)')
%     title(fileName(1:16));
%     subplot(3,1,3)
%     imagesc(paramC.time,1:length(trialNoNonStimGood),thetaPhaseCueNonStimGoodL(indSpeed,:));
%     set(gca,'XLim',[0 3])
%     ylabel('Trial no. (Cue)')
%     xlabel('Time (s)')
    
%     [~,indStartPhase] = sort(thetaPhaseRunNonStimGoodL(:,nSampBef+1));
%     [figNew,pos] = CreateFig();
%     set(0,'Units','pixels') 
%     set(figure(figNew),'OuterPosition',[pos(1) pos(2) pos(3) pos(4)])
%     subplot(3,1,1)
%     imagesc(paramC.timeIncBef,1:length(trialNoNonStimGood),runSpeedNonStimGood(indStartPhase,:));
%     set(gca,'XLim',[-0.5 3])
%     ylabel('Trial no. (RunOnset speed)')
%     subplot(3,1,2)
%     imagesc(paramC.timeIncBef,1:length(trialNoNonStimGood),thetaPhaseRunNonStimGoodL(indStartPhase,:));
%     set(gca,'XLim',[-0.5 3])
%     ylabel('Trial no. (RunOnset)')
%     title(fileName(1:16));
%     subplot(3,1,3)
%     imagesc(paramC.time,1:length(trialNoNonStimGood),thetaPhaseCueNonStimGoodL(indStartPhase,:));
%     set(gca,'XLim',[0 3])
%     ylabel('Trial no. (Cue)')
%     xlabel('Time (s)')
    
%     [figNew,pos] = CreateFig();
%     set(0,'Units','pixels') 
%     set(figure(figNew),'OuterPosition',[pos(1) pos(2) pos(3) pos(4)])
%     hist(thetaPhaseRunNonStimGoodL(:,nSampBef+1),-pi:pi/20:pi)
end

function runSpeed = speedAlignedRun(trialsRun,trialNo,trialLen)
    runSpeed = zeros(length(trialNo),trialLen);
    for i = 1:length(trialNo)
        tmpSpeed = [trialsRun.speed_MMsecBef{trialNo(i)}',...
            trialsRun.speed_MMsec{trialNo(i)}'];
        tmpSpeed = tmpSpeed(tmpSpeed >= 0);
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
        tmpAcc = tmpAcc(tmpSpeed >= 0);
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