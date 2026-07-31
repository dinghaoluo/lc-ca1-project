function meanPopActivityCorrDiffBeh(pathAnal,onlyRun,type1,type2)
% compare the passive no cue recordings with the active licking no cue
% behavior, and look at the population activity correlation over trials
% default analysis path ['Z:\Raphael_tests\mice_expdata\Analysis\'];
% e.g.  meanActivityCorrDiffBeh('Z:\Raphael_tests\mice_expdata\Analysis\',3,2)
% type1/type2 = 1, passive licking; 2, active licking; 3, passive no cue

    load([pathAnal 'popActivityCorr_RecList' num2str(type1) '_Run1.mat']);
    recList1 = recListPopActivityCorr;
    
    load([pathAnal 'popActivityCorr_RecList' num2str(type2) '_Run1.mat']);
    recList2 = recListPopActivityCorr;
    
    trialTypeDiffPopActCorr.meanMeanCorr = [mean(recList1.meanCorr),...
        mean(recList2.meanCorr)];
    trialTypeDiffPopActCorr.stdMeanCorr = ...
        [std(recList1.meanCorr)/sqrt(length(recList1.meanCorr)),...
        std(recList2.meanCorr)/sqrt(length(recList2.meanCorr))];
    trialTypeDiffPopActCorr.pMeanCorr = ...
        ranksum(recList1.meanCorr,recList2.meanCorr);
   
    fullPath = ['Z:\Raphael_tests\mice_expdata\Analysis\'];
    fullPath = [fullPath 'compRecList' num2str(type1) '-' num2str(type2) ...
        '_popActivityCorr_Run' num2str(onlyRun) '.mat'];
    save(fullPath,'trialTypeDiffPopActCorr');
    
    plotBars(recList1.meanCorr,recList2.meanCorr,...
        trialTypeDiffPopActCorr.meanMeanCorr,...
        trialTypeDiffPopActCorr.stdMeanCorr,'','mean Pop. Corr.', ...
        ['Mean population correlation: p=' ...
        num2str(trialTypeDiffPopActCorr.pMeanCorr)]);
    print('-painters','-dpdf',[pathAnal 'TrialTypeCompMeanPopActivityCorr' num2str(type1) '-' num2str(type2)],'-r600');
   
end

function plotBars(data1,data2,mean,std,x1,y1,t)
    fig = figure;
    set(0,'Units','pixels') 
    set(figure(fig),'OuterPosition',...
        [500 500 210 280])
    fig.Renderer = 'Painters';
        
    h = bar([1,2],mean,0.5);
    set(h,'EdgeColor',[0.3 0.3 0.3],'FaceColor',[187 189 192]/255);
    hold
    
    h = errorbar([1,2],mean,std);
    set(h,'Marker','.','MarkerSize',0.1,'Color',[0 0 0],'LineStyle','none')
    
    h = plot(1.05+0.1*rand(1,length(data1)),data1,'o');
    set(h,'MarkerSize',5,'Color',[167 169 171]/255);
    
    h = plot(2.05+0.1*rand(1,length(data2)),data2,'o');
    set(h,'MarkerSize',5,'Color',[27 117 187]/255);
    set(gca,'XLim',[0.5 2.5],'YLim',[0 max([data1 data2])+0.01]);
    
    xlabel(x1);
    ylabel(y1);
    title(t);
end
