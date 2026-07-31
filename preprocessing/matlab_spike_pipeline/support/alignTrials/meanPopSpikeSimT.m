function meanPopSpikeSimT(path,fileName,onlyRun,mazeSess,intervalT)
% single neuron level mean spike correlation across trials
    
    if(nargin == 4)
        intervalT = 0;
    end
    fullPath = [path fileName '_popSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    fileNameCorr = [fileName '_meanPopSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _popSimT_Run file does not exist');
        return;
    end
    load(fullPath,'popSimTRun','popSimTRew','popSimTCue');
    
    fullPath = [path fileName '_PeakFRAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];    
    if(exist(fullPath) == 0)
        disp('The _PeakFRAligned file does not exist');
        return;
    end
    load(fullPath,'trialNoNonStimGood','trialNoNonStimBad','trialNoStim','trialNoStimCtrl','pulseMeth');
    
    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    trialNo = size(popSimTRun,1);
    nElem = (trialNo*trialNo-trialNo)/2;
    
    indBadBeh = behPar.indTrBadBeh;
    indGoodTr = indBadBeh == 0;
    indBadTr = indBadBeh == 1;
    nGoodTr = sum(indGoodTr);
    nBadTr = sum(indBadTr);
    nElemGood = (nGoodTr*nGoodTr-nGoodTr)/2;
    nElemBad = (nBadTr*nBadTr-nBadTr)/2;
       
    corrArr = triu(popSimTRun,1);
    corrArr = corrArr(:);        
    meanPopSimTRun.mean = sum(corrArr(isnan(corrArr) == 0))/nElem;
      
    corrArr = triu(popSimTRew,1);
    corrArr = corrArr(:);  
    meanPopSimTRew.mean = sum(corrArr(isnan(corrArr) == 0))/nElem;
        
    corrArr = triu(popSimTCue,1);
    corrArr = corrArr(:);
    meanPopSimTCue.mean = sum(corrArr(isnan(corrArr) == 0))/nElem;
        
    goodTrCorr = popSimTRun(indGoodTr,indGoodTr);
    goodTrCorr = triu(goodTrCorr,1);
    goodTrCorr = goodTrCorr(:);
    meanPopSimTRun.meanGood = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGood;
        
    goodTrCorr = triu(popSimTRew(indGoodTr,indGoodTr),1);
    goodTrCorr = goodTrCorr(:);
    meanPopSimTRew.meanGood = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGood;
        
    goodTrCorr = triu(popSimTCue(indGoodTr,indGoodTr),1);
    goodTrCorr = goodTrCorr(:);
    meanPopSimTCue.meanGood = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGood;
        
    badTrCorr = triu(popSimTRun(indBadTr,indBadTr),1);
    badTrCorr = badTrCorr(:);
    meanPopSimTRun.meanBad = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBad;
        
    badTrCorr = triu(popSimTRew(indBadTr,indBadTr),1);
    badTrCorr = badTrCorr(:);
    meanPopSimTRew.meanBad = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBad;
        
    badTrCorr = triu(popSimTCue(indBadTr,indBadTr),1);
    badTrCorr = badTrCorr(:);
    meanPopSimTCue.meanBad = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBad;
    
    %% added by YIngxue on 2/14/2021
    goodNoStimTrCorr = popSimTRun(trialNoNonStimGood,trialNoNonStimGood);
    goodNoStimTrCorr = triu(goodNoStimTrCorr,1);
    goodNoStimTrCorr = goodNoStimTrCorr(:);
    nNoStimGoodTr = length(trialNoNonStimGood);
    nElemNoStimGood = (nNoStimGoodTr*nNoStimGoodTr-nNoStimGoodTr)/2;
    meanPopSimTRun.meanNoStimGood = sum(goodNoStimTrCorr(isnan(goodNoStimTrCorr) == 0))/nElemNoStimGood;
    
    badNoStimTrCorr = popSimTRun(trialNoNonStimBad,trialNoNonStimBad);
    badNoStimTrCorr = triu(badNoStimTrCorr,1);
    badNoStimTrCorr = badNoStimTrCorr(:);
    nNoStimBadTr = length(trialNoNonStimBad);
    nElemNoStimBad = (nNoStimBadTr*nNoStimBadTr-nNoStimBadTr)/2;
    meanPopSimTRun.meanNoStimBad = sum(badNoStimTrCorr(isnan(badNoStimTrCorr) == 0))/nElemNoStimBad;
    
    meanPopSimTRun.meanStim = [];
    meanPopSimTRun.meanStimCtrl = [];
    for i = 1:length(pulseMeth) 
        stimTrCorr = popSimTRun(trialNoStim{i},trialNoStim{i});
        stimTrCorr = triu(stimTrCorr,1);
        stimTrCorr = stimTrCorr(:);
        nStimTr = length(trialNoStim{i});
        nElemStim = (nStimTr*nStimTr-nStimTr)/2;
        meanPopSimTRun.meanStim{i} = sum(stimTrCorr(isnan(stimTrCorr) == 0))/nElemStim;

        stimCtrlTrCorr = popSimTRun(trialNoStimCtrl{i},trialNoStimCtrl{i});
        stimCtrlTrCorr = triu(stimCtrlTrCorr,1);
        stimCtrlTrCorr = stimCtrlTrCorr(:);
        nStimCtrlTr = length(trialNoStimCtrl{i});
        nElemStimCtrl = (nStimCtrlTr*nStimCtrlTr-nStimCtrlTr)/2;
        meanPopSimTRun.meanStimCtrl{i} = sum(stimCtrlTrCorr(isnan(stimCtrlTrCorr) == 0))/nElemStimCtrl;
    end

    save([path fileNameCorr],'meanPopSimTRun','meanPopSimTRew','meanPopSimTCue');
end

function plotCompCorr(x,y,xlab,ylab)
    figure
    plot(x,y,'ro');
    hold on;
    corrAll = [x y];
    maxCorr = max(corrAll);
    minCorr = min(corrAll);
    plot([minCorr maxCorr],[minCorr maxCorr],'k:');
    xlim(gca,[minCorr maxCorr]);
    ylim(gca,[minCorr maxCorr]);
    xlabel(xlab);
    ylabel(ylab);
end
