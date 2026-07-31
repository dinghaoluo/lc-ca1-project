function meanPopSpikeCorrTLastToCurTr(path,fileName,onlyRun,mazeSess,intervalT,intervalTMin)
% single neuron level mean spike correlation across trials
    
    if(nargin == 4)
        intervalT = 0;
    end
    fullPath = [path fileName '_popCorrT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    fileNameCorr = [fileName '_meanPopCorrT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
             num2str(intervalTMin) '-' num2str(intervalT) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _popCorrTAligned_Run file does not exist');
        return;
    end
    load(fullPath,'popCorrTRun_LasttoCurTr','popCorrTRew_LasttoCurTr','popCorrTCue_LasttoCurTr');
    
    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    trialNo = size(popCorrTRun_LasttoCurTr,1);
    nElem = (trialNo*trialNo-trialNo)/2;
    
    indBadBeh = behPar.indTrBadBeh;
    indGoodTr = indBadBeh == 0;
    indBadTr = indBadBeh == 1;
    nGoodTr = sum(indGoodTr);
    nBadTr = sum(indBadTr);
    nElemGood = (nGoodTr*nGoodTr-nGoodTr)/2;
    nElemBad = (nBadTr*nBadTr-nBadTr)/2;
       
    corrArr = triu(popCorrTRun_LasttoCurTr,1);
    corrArr = corrArr(:);        
    meanPopCorrTRun_LasttoCurTr.mean = sum(corrArr(isnan(corrArr) == 0))/nElem;
      
    corrArr = triu(popCorrTRew_LasttoCurTr,1);
    corrArr = corrArr(:);  
    meanPopCorrTRew_LasttoCurTr.mean = sum(corrArr(isnan(corrArr) == 0))/nElem;
        
    corrArr = triu(popCorrTCue_LasttoCurTr,1);
    corrArr = corrArr(:);
    meanPopCorrTCue_LasttoCurTr.mean = sum(corrArr(isnan(corrArr) == 0))/nElem;
        
    goodTrCorr = popCorrTRun_LasttoCurTr(indGoodTr,indGoodTr);
    goodTrCorr = triu(goodTrCorr,1);
    goodTrCorr = goodTrCorr(:);
    meanPopCorrTRun_LasttoCurTr.meanGood = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGood;
        
    goodTrCorr = triu(popCorrTRew_LasttoCurTr(indGoodTr,indGoodTr),1);
    goodTrCorr = goodTrCorr(:);
    meanPopCorrTRew_LasttoCurTr.meanGood = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGood;
        
    goodTrCorr = triu(popCorrTCue_LasttoCurTr(indGoodTr,indGoodTr),1);
    goodTrCorr = goodTrCorr(:);
    meanPopCorrTCue_LasttoCurTr.meanGood = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGood;
        
    badTrCorr = triu(popCorrTRun_LasttoCurTr(indBadTr,indBadTr),1);
    badTrCorr = badTrCorr(:);
    meanPopCorrTRun_LasttoCurTr.meanBad = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBad;
        
    badTrCorr = triu(popCorrTRew_LasttoCurTr(indBadTr,indBadTr),1);
    badTrCorr = badTrCorr(:);
    meanPopCorrTRew_LasttoCurTr.meanBad = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBad;
        
    badTrCorr = triu(popCorrTCue_LasttoCurTr(indBadTr,indBadTr),1);
    badTrCorr = badTrCorr(:);
    meanPopCorrTCue_LasttoCurTr.meanBad = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBad;

    save([path fileNameCorr],'meanPopCorrTRun_LasttoCurTr','meanPopCorrTRew_LasttoCurTr','meanPopCorrTCue_LasttoCurTr');
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
