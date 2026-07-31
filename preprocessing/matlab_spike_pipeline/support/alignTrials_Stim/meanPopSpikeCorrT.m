function meanPopSpikeCorrT(path,fileName,onlyRun,mazeSess,intervalT)
% single neuron level mean spike correlation across trials
    
    if(nargin == 4)
        intervalT = 0;
    end
    fullPath = [path fileName '_popCorrT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    fileNameCorr = [fileName '_meanPopCorrT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
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
    
    trialNo = size(popCorrTRun,1);
    nElem = (trialNo*trialNo-trialNo)/2;
    
    indBadBeh = behPar.indTrBadBeh;
    indGoodTr = indBadBeh == 0;
    indBadTr = indBadBeh == 1;
    nGoodTr = sum(indGoodTr);
    nBadTr = sum(indBadTr);
    nElemGood = (nGoodTr*nGoodTr-nGoodTr)/2;
    nElemBad = (nBadTr*nBadTr-nBadTr)/2;
       
    corrArr = triu(popCorrTRun,1);
    corrArr = corrArr(:);        
    meanPopCorrTRun.mean = sum(corrArr(isnan(corrArr) == 0))/nElem;
      
    corrArr = triu(popCorrTRew,1);
    corrArr = corrArr(:);  
    meanPopCorrTRew.mean = sum(corrArr(isnan(corrArr) == 0))/nElem;
        
    corrArr = triu(popCorrTCue,1);
    corrArr = corrArr(:);
    meanPopCorrTCue.mean = sum(corrArr(isnan(corrArr) == 0))/nElem;
        
    goodTrCorr = popCorrTRun(indGoodTr,indGoodTr);
    goodTrCorr = triu(goodTrCorr,1);
    goodTrCorr = goodTrCorr(:);
    meanPopCorrTRun.meanGood = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGood;
        
    goodTrCorr = triu(popCorrTRew(indGoodTr,indGoodTr),1);
    goodTrCorr = goodTrCorr(:);
    meanPopCorrTRew.meanGood = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGood;
        
    goodTrCorr = triu(popCorrTCue(indGoodTr,indGoodTr),1);
    goodTrCorr = goodTrCorr(:);
    meanPopCorrTCue.meanGood = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGood;
        
    badTrCorr = triu(popCorrTRun(indBadTr,indBadTr),1);
    badTrCorr = badTrCorr(:);
    meanPopCorrTRun.meanBad = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBad;
        
    badTrCorr = triu(popCorrTRew(indBadTr,indBadTr),1);
    badTrCorr = badTrCorr(:);
    meanPopCorrTRew.meanBad = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBad;
        
    badTrCorr = triu(popCorrTCue(indBadTr,indBadTr),1);
    badTrCorr = badTrCorr(:);
    meanPopCorrTCue.meanBad = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBad;

    save([path fileNameCorr],'meanPopCorrTRun','meanPopCorrTRew','meanPopCorrTCue');
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
