function meanPopSpikeSimTStim(path,fileName,onlyRun,mazeSess,intervalT)
% population level mean spike correlation across trials

    if(nargin == 4)
        intervalT = 0;
    end
    fullPath = [path fileName '_popSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    fileNameCorr = [fileName '_meanPopCorrTStim_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _popSimTAligned_Run file does not exist');
        return;
    end
    load(fullPath,'popSimTRun','popSimTRew','popSimTCue');
    
    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    fullPath = [path fileName '_PeakFRAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];    
    if(exist(fullPath) == 0)
        disp('The _PeakFRAligned file does not exist');
        return;
    end
    load(fullPath,'trialNoNonStimGood','trialNoNonStimBad','trialNoStim','trialNoStimCtrl','pulseMeth');
        
    nPulseMeth = length(pulseMeth);
    trialNoNonStim = [trialNoNonStimGood; trialNoNonStimBad];
    nTrialNonStim = length(trialNoNonStim);
    nElemNonStim = (nTrialNonStim*nTrialNonStim-nTrialNonStim)/2;
    
    nTrialNonStimGood = length(trialNoNonStimGood);
    nElemNonStimGood = (nTrialNonStimGood*nTrialNonStimGood-nTrialNonStimGood)/2;
    
    nTrialNonStimBad = length(trialNoNonStimBad);
    nElemNonStimBad = (nTrialNonStimBad*nTrialNonStimBad-nTrialNonStimBad)/2;
    
    nTrialStim = zeros(1,length(pulseMeth));
    nTrialStimCtrl = zeros(1,length(pulseMeth));
    nElemStim = zeros(1,length(pulseMeth));
    nElemStimCtrl = zeros(1,length(pulseMeth));
    for i = 1:length(pulseMeth)
        nTrialStim(i) = length(trialNoStim{i});
        nElemStim(i) = (nTrialStim(i)*nTrialStim(i)-nTrialStim(i))/2;
        
        nTrialStimCtrl(i) = length(trialNoStimCtrl{i});
        nElemStimCtrl(i) = (nTrialStimCtrl(i)*nTrialStimCtrl(i)-nTrialStimCtrl(i))/2;
    end
    
    meanPopSimT = struct('meanNonStim',0,...
                       ...
                       'meanNonStimGood',0,...
                       ...
                       'meanNonStimBad',0,...
                       ...
                       'meanStim',zeros(nPulseMeth,1),...
                       ...
                       'meanStimCtrl',zeros(nPulseMeth,1));
    
    meanPopSimTRun = meanPopSimT;
    meanPopSimTRew = meanPopSimT;
    meanPopSimTCue = meanPopSimT;
    
    %% non-stim trials
    corrArr = triu(popSimTRun(trialNoNonStim,trialNoNonStim),1);
    corrArr = corrArr(:);        
    meanPopSimTRun.meanNonStim = sum(corrArr(isnan(corrArr) == 0))/nElemNonStim;

    corrArr = triu(popSimTRew(trialNoNonStim,trialNoNonStim),1);
    corrArr = corrArr(:);        
    meanPopSimTRew.meanNonStim = sum(corrArr(isnan(corrArr) == 0))/nElemNonStim;

    corrArr = triu(popSimTCue(trialNoNonStim,trialNoNonStim),1);
    corrArr = corrArr(:);        
    meanPopSimTCue.meanNonStim = sum(corrArr(isnan(corrArr) == 0))/nElemNonStim;

    %% non-stim good trials
    corrArr = triu(popSimTRun(trialNoNonStimGood,trialNoNonStimGood),1);
    corrArr = corrArr(:);        
    meanPopSimTRun.meanNonStimGood = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimGood;

    corrArr = triu(popSimTRew(trialNoNonStimGood,trialNoNonStimGood),1);
    corrArr = corrArr(:);        
    meanPopSimTRew.meanNonStimGood = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimGood;

    corrArr = triu(popSimTCue(trialNoNonStimGood,trialNoNonStimGood),1);
    corrArr = corrArr(:);        
    meanPopSimTCue.meanNonStimGood = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimGood;

    %% non-stim bad trials
    corrArr = triu(popSimTRun(trialNoNonStimBad,trialNoNonStimBad),1);
    corrArr = corrArr(:);        
    meanPopSimTRun.meanNonStimBad = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimBad;

    corrArr = triu(popSimTRew(trialNoNonStimBad,trialNoNonStimBad),1);
    corrArr = corrArr(:);        
    meanPopSimTRew.meanNonStimBad = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimBad;

    corrArr = triu(popSimTCue(trialNoNonStimBad,trialNoNonStimBad),1);
    corrArr = corrArr(:);        
    meanPopSimTCue.meanNonStimBad = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimBad;

    for i = 1:length(pulseMeth)
        %% stim trials
        corrArr = triu(popSimTRun(trialNoStim{i},trialNoStim{i}),1);
        corrArr = corrArr(:);        
        meanPopSimTRun.meanStim(i) = sum(corrArr(isnan(corrArr) == 0))/nElemStim(i);

        corrArr = triu(popSimTRew(trialNoStim{i},trialNoStim{i}),1);
        corrArr = corrArr(:);        
        meanPopSimTRew.meanStim(i) = sum(corrArr(isnan(corrArr) == 0))/nElemStim(i);

        corrArr = triu(popSimTCue(trialNoStim{i},trialNoStim{i}),1);
        corrArr = corrArr(:);        
        meanPopSimTCue.meanStim(i) = sum(corrArr(isnan(corrArr) == 0))/nElemStim(i);

        %% stim ctrl trials
        corrArr = triu(popSimTRun(trialNoStimCtrl{i},trialNoStimCtrl{i}),1);
        corrArr = corrArr(:);        
        meanPopSimTRun.meanStimCtrl(i) = sum(corrArr(isnan(corrArr) == 0))/nElemStimCtrl(i);

        corrArr = triu(popSimTRew(trialNoStimCtrl{i},trialNoStimCtrl{i}),1);
        corrArr = corrArr(:);        
        meanPopSimTRew.meanStimCtrl(i) = sum(corrArr(isnan(corrArr) == 0))/nElemStimCtrl(i);

        corrArr = triu(popSimTCue(trialNoStimCtrl{i},trialNoStimCtrl{i}),1);
        corrArr = corrArr(:);        
        meanPopSimTCue.meanStimCtrl(i) = sum(corrArr(isnan(corrArr) == 0))/nElemStimCtrl(i);

    end
    
    save([path fileNameCorr],'meanPopSimTRun','meanPopSimTRew','meanPopSimTCue');
end
