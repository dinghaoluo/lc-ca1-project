function meanPopSpikeSimTLastToCurTr(path,fileName,onlyRun,mazeSess,intervalT,intervalTMin)
% single neuron level mean spike correlation across trials
    
    if(nargin == 4)
        intervalT = 0;
    end
    fullPath = [path fileName '_popSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    fileNameCorr = [fileName '_meanPopSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalTMin) '-' num2str(intervalT) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _popSimT_Run file does not exist');
        return;
    end
    load(fullPath,'popSimTRun_LasttoCurTr','popSimTRew_LasttoCurTr','popSimTCue_LasttoCurTr');
    
    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    trialNo = size(popSimTRun_LasttoCurTr,1);
    nElem = (trialNo*trialNo-trialNo)/2;
    
    indBadBeh = behPar.indTrBadBeh;
    indGoodTr = indBadBeh == 0;
    indBadTr = indBadBeh == 1;
    nGoodTr = sum(indGoodTr);
    nBadTr = sum(indBadTr);
    nElemGood = (nGoodTr*nGoodTr-nGoodTr)/2;
    nElemBad = (nBadTr*nBadTr-nBadTr)/2;
       
    corrArr = triu(popSimTRun_LasttoCurTr,1);
    corrArr = corrArr(:);        
    meanPopSimTRun_LasttoCurTr.mean = sum(corrArr(isnan(corrArr) == 0))/nElem;
      
    corrArr = triu(popSimTRew_LasttoCurTr,1);
    corrArr = corrArr(:);  
    meanPopSimTRew_LasttoCurTr.mean = sum(corrArr(isnan(corrArr) == 0))/nElem;
        
    corrArr = triu(popSimTCue_LasttoCurTr,1);
    corrArr = corrArr(:);
    meanPopSimTCue_LasttoCurTr.mean = sum(corrArr(isnan(corrArr) == 0))/nElem;
        
    goodTrCorr = popSimTRun_LasttoCurTr(indGoodTr,indGoodTr);
    goodTrCorr = triu(goodTrCorr,1);
    goodTrCorr = goodTrCorr(:);
    meanPopSimTRun_LasttoCurTr.meanGood = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGood;
        
    goodTrCorr = triu(popSimTRew_LasttoCurTr(indGoodTr,indGoodTr),1);
    goodTrCorr = goodTrCorr(:);
    meanPopSimTRew_LasttoCurTr.meanGood = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGood;
        
    goodTrCorr = triu(popSimTCue_LasttoCurTr(indGoodTr,indGoodTr),1);
    goodTrCorr = goodTrCorr(:);
    meanPopSimTCue_LasttoCurTr.meanGood = sum(goodTrCorr(isnan(goodTrCorr) == 0))/nElemGood;
        
    badTrCorr = triu(popSimTRun_LasttoCurTr(indBadTr,indBadTr),1);
    badTrCorr = badTrCorr(:);
    meanPopSimTRun_LasttoCurTr.meanBad = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBad;
        
    badTrCorr = triu(popSimTRew_LasttoCurTr(indBadTr,indBadTr),1);
    badTrCorr = badTrCorr(:);
    meanPopSimTRew_LasttoCurTr.meanBad = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBad;
        
    badTrCorr = triu(popSimTCue_LasttoCurTr(indBadTr,indBadTr),1);
    badTrCorr = badTrCorr(:);
    meanPopSimTCue_LasttoCurTr.meanBad = sum(badTrCorr(isnan(badTrCorr) == 0))/nElemBad;

    save([path fileNameCorr],'meanPopSimTRun_LasttoCurTr','meanPopSimTRew_LasttoCurTr','meanPopSimTCue_LasttoCurTr');
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
