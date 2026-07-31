function popSimilarityAllRecCompExp(onlyRun)
% compare the population corrT for different experimental paradigms (passive vs active) 

    load('Z:\Yingxue\DataAnalysisRaphi\PopCorrTAlignedActive_AllRec_Stat.mat',...
        'popSimRunAll','popSimRewAll','popSimCueAll');
    popSimRunAllAct = popSimRunAll;
    popSimRewAllAct = popSimRewAll;
    popSimCueAllAct = popSimCueAll;
    
    load('Z:\Yingxue\DataAnalysisRaphi\PopCorrTAlignedNoCue_AllRec_Stat.mat',...
        'popSimRunAll','popSimRewAll','popSimCueAll');
    popSimRunAllNoCue = popSimRunAll;
    popSimRewAllNoCue = popSimRewAll;
    popSimCueAllNoCue = popSimCueAll;
    
    popSimRunAllAL_NoCue.pRS_CorrTGood = ranksum(popSimRunAllAct.CorrTGood,...
        popSimRunAllNoCue.CorrTGood);
    popSimRunAllAL_NoCue.pRS_CorrTBad = ranksum(popSimRunAllAct.CorrTBad,...
        popSimRunAllNoCue.CorrTBad);
    popSimRunAllAL_NoCue.pRS_CorrT = ranksum(popSimRunAllAct.CorrT,...
        popSimRunAllNoCue.CorrT);

    plotSimCompxy(popSimRunAllNoCue.CorrTGood,...
        popSimRunAllAct.CorrTGood,...
        popSimRunAllNoCue.meanCorrTGood,...
        popSimRunAllAct.meanCorrTGood,...
        popSimRunAllAL_NoCue.pRS_CorrTGood,...
        'Pop. corr. good tr.',0.6);
    fileName1 = ['PopCorrTAligned_Good_AL_NoCue_AllRec'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');

    plotSimCompxy(popSimRunAllNoCue.CorrTBad,...
        popSimRunAllAct.CorrTBad,...
        popSimRunAllNoCue.meanCorrTBad,...
        popSimRunAllAct.meanCorrTBad,...
        popSimRunAllAL_NoCue.pRS_CorrTBad,...
        'Pop. corr. bad tr.',0.6);
    fileName1 = ['PopCorrTAligned_Bad_AL_NoCue_AllRec'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');

    plotSimCompxy(popSimRunAllNoCue.CorrT,...
        popSimRunAllAct.CorrT,...
        popSimRunAllNoCue.meanCorrT,...
        popSimRunAllAct.meanCorrT,...
        popSimRunAllAL_NoCue.pRS_CorrT,...
        'Pop. corr.',0.6);
    fileName1 = ['PopCorrTAligned_AL_NoCue_AllRec'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');
end

function plotSimCompxy(x,y,xm,ym,pxy,yl,maxXY)
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 280 280])

    h = bar([1 2],[xm,ym]);
    set(h,'FaceColor',[0.7 0.7 0.9],'EdgeColor',[0.5 0.5 0.5]);
    if(isempty(maxXY))
        maxXY = max([x,y]);
    end
    hold on;
    h = plot(ones(1,length(x)),x,'k.');
    set(h,'MarkerSize',11);
    h = plot(2*ones(1,length(y)),y,'k.');
    set(h,'MarkerSize',11);
    set(gca,'YLim',[0 maxXY]);
    title(['p12 = ' num2str(pxy) ]);
    ylabel(yl)
end