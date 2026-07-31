function meanPopSpikeCorrTStim(path,fileName,onlyRun,mazeSess,intervalT)
% population level mean spike correlation across trials

    if(nargin == 4)
        intervalT = 0;
    end
    fullPath = [path fileName '_popCorrT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    fileNameCorr = [fileName '_meanPopCorrTStim_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _popCorrTAligned_Run file does not exist');
        return;
    end
    load(fullPath,'popCorrTRun','popCorrTRew','popCorrTCue');
    
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
    
    meanPopCorrT = struct('meanNonStim',0,...
                       ...
                       'meanNonStimGood',0,...
                       ...
                       'meanNonStimBad',0,...
                       ...
                       'meanStim',zeros(nPulseMeth,1),...
                       ...
                       'meanStimCtrl',zeros(nPulseMeth,1));
    
    meanPopCorrTRun = meanPopCorrT;
    meanPopCorrTRew = meanPopCorrT;
    meanPopCorrTCue = meanPopCorrT;
    
    %% non-stim trials
    corrArr = triu(popCorrTRun(trialNoNonStim,trialNoNonStim),1);
    corrArr = corrArr(:);        
    meanPopCorrTRun.meanNonStim = sum(corrArr(isnan(corrArr) == 0))/nElemNonStim;

    corrArr = triu(popCorrTRew(trialNoNonStim,trialNoNonStim),1);
    corrArr = corrArr(:);        
    meanPopCorrTRew.meanNonStim = sum(corrArr(isnan(corrArr) == 0))/nElemNonStim;

    corrArr = triu(popCorrTCue(trialNoNonStim,trialNoNonStim),1);
    corrArr = corrArr(:);        
    meanPopCorrTCue.meanNonStim = sum(corrArr(isnan(corrArr) == 0))/nElemNonStim;

    %% non-stim good trials
    corrArr = triu(popCorrTRun(trialNoNonStimGood,trialNoNonStimGood),1);
    corrArr = corrArr(:);        
    meanPopCorrTRun.meanNonStimGood = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimGood;

    corrArr = triu(popCorrTRew(trialNoNonStimGood,trialNoNonStimGood),1);
    corrArr = corrArr(:);        
    meanPopCorrTRew.meanNonStimGood = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimGood;

    corrArr = triu(popCorrTCue(trialNoNonStimGood,trialNoNonStimGood),1);
    corrArr = corrArr(:);        
    meanPopCorrTCue.meanNonStimGood = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimGood;

    %% non-stim bad trials
    corrArr = triu(popCorrTRun(trialNoNonStimBad,trialNoNonStimBad),1);
    corrArr = corrArr(:);        
    meanPopCorrTRun.meanNonStimBad = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimBad;

    corrArr = triu(popCorrTRew(trialNoNonStimBad,trialNoNonStimBad),1);
    corrArr = corrArr(:);        
    meanPopCorrTRew.meanNonStimBad = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimBad;

    corrArr = triu(popCorrTCue(trialNoNonStimBad,trialNoNonStimBad),1);
    corrArr = corrArr(:);        
    meanPopCorrTCue.meanNonStimBad = sum(corrArr(isnan(corrArr) == 0))/nElemNonStimBad;

    for i = 1:length(pulseMeth)
        %% stim trials
        corrArr = triu(popCorrTRun(trialNoStim{i},trialNoStim{i}),1);
        corrArr = corrArr(:);        
        meanPopCorrTRun.meanStim(i) = sum(corrArr(isnan(corrArr) == 0))/nElemStim(i);

        corrArr = triu(popCorrTRew(trialNoStim{i},trialNoStim{i}),1);
        corrArr = corrArr(:);        
        meanPopCorrTRew.meanStim(i) = sum(corrArr(isnan(corrArr) == 0))/nElemStim(i);

        corrArr = triu(popCorrTCue(trialNoStim{i},trialNoStim{i}),1);
        corrArr = corrArr(:);        
        meanPopCorrTCue.meanStim(i) = sum(corrArr(isnan(corrArr) == 0))/nElemStim(i);

        %% stim ctrl trials
        corrArr = triu(popCorrTRun(trialNoStimCtrl{i},trialNoStimCtrl{i}),1);
        corrArr = corrArr(:);        
        meanPopCorrTRun.meanStimCtrl(i) = sum(corrArr(isnan(corrArr) == 0))/nElemStimCtrl(i);

        corrArr = triu(popCorrTRew(trialNoStimCtrl{i},trialNoStimCtrl{i}),1);
        corrArr = corrArr(:);        
        meanPopCorrTRew.meanStimCtrl(i) = sum(corrArr(isnan(corrArr) == 0))/nElemStimCtrl(i);

        corrArr = triu(popCorrTCue(trialNoStimCtrl{i},trialNoStimCtrl{i}),1);
        corrArr = corrArr(:);        
        meanPopCorrTCue.meanStimCtrl(i) = sum(corrArr(isnan(corrArr) == 0))/nElemStimCtrl(i);

    end
    
    save([path fileNameCorr],'meanPopCorrTRun','meanPopCorrTRew','meanPopCorrTCue');
end
