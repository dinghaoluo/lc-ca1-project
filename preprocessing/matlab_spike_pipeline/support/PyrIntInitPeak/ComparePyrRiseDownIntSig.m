function ComparePyrRiseDownIntSig(pathAnal,timeStepRun,InitAllPyr,InitAllInt,PyrRise,PyrDown,IntRise)

    %% compare pyramidal cells profile with each type of interneuron
    for i = 1:max(InitAllInt.idxC2)
        
        meanAvgFRProfile = mean(InitAllPyr.avgFRProfileNorm(PyrRise.idxRise,:));
        meanAvgFRProfile1 = mean(InitAllPyr.avgFRProfileNorm(PyrDown.idxDown,:));
        avgFRProfileIntC = InitAllInt.avgFRProfileNorm(IntRise.idxRise{i},:);
        meanAvgFRProfileInt = mean(avgFRProfileIntC);
        
        indTmp = timeStepRun >=-1 & timeStepRun <= 4 ;
        plotMeanFRProfilePyrVsInt(timeStepRun(indTmp),...
            (meanAvgFRProfile(indTmp)-min(meanAvgFRProfile(indTmp)))...
            /(max(meanAvgFRProfile(indTmp))-min(meanAvgFRProfile(indTmp))),...
            (meanAvgFRProfileInt(indTmp)-min(meanAvgFRProfileInt(indTmp)))...
            /(max(meanAvgFRProfileInt(indTmp))-min(meanAvgFRProfileInt(indTmp))),...
            ['Norm. FR PyrInc Vs Int C' num2str(i)],...
            ['meanFRProfileNormPyrIncVsIntC' num2str(i)],pathAnal,[]);
        
        plotMeanFRProfilePyrVsInt(timeStepRun(indTmp),...
            (meanAvgFRProfile1(indTmp)-min(meanAvgFRProfile1(indTmp)))...
            /(max(meanAvgFRProfile1(indTmp))-min(meanAvgFRProfile1(indTmp))),...
            (meanAvgFRProfileInt(indTmp)-min(meanAvgFRProfileInt(indTmp)))...
            /(max(meanAvgFRProfileInt(indTmp))-min(meanAvgFRProfileInt(indTmp))),...
            ['Norm. FR PyrDec Vs Int C' num2str(i)],...
            ['meanFRProfileNormPyrDecVsIntC' num2str(i)],pathAnal,[]);
        
        plotAvgFRProfileCmp(timeStepRun,...
            InitAllPyr.avgFRProfileNorm(PyrRise.idxRise,:),avgFRProfileIntC,...
            ['Norm FR PyrInc vs. Int C' num2str(i)],...
            ['FRProfileNormPyrIncVsIntC' num2str(i)],...
            pathAnal,[0 1],[{'Pyr'} {'Int'}])
    end
    
end

function plotMeanFRProfilePyrVsInt(timeStepRun,avgFRProfilex,avgFRProfiley,yl,fileName,pathAnal,ylimit)
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 400 200]);
    
    h = plot(timeStepRun,avgFRProfilex);
    set(h,'LineWidth',1, 'Color',[ 39 169 225]./255);
    hold on;
    h = plot(timeStepRun,avgFRProfiley);
    set(h,'LineWidth',1, 'Color',[167 169  171]./255);
    
    if(~isempty(ylimit))
        h = plot([0 0],ylimit,'r-');
    else
        h = plot([0 0],[0 1],'r-');
    end
    set(h,'LineWidth',1)
    set(gca,'XLim',[-1 4]);
    
    xlabel('Time (s)')
    ylabel(yl)
    
    fileName1 = [pathAnal fileName];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end