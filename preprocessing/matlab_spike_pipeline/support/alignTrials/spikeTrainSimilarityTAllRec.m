function spikeTrainSimilarityTAllRec(onlyRun)
% compare the spike train similarity for different alignment conditions over all the recordings 

    GlobalConst;
    intervalT = 20;
    intervalD = 1800;
    tc = 500; % ms
    cost = 2/(tc/1000*sampleFq);
    
    thrCorrT = -1; %0.03;
    
    RecordingList;
    
    spikeSimRunAll.SimTNonZeroGood = [];
    spikeSimRunAll.SimTGood = [];
    spikeSimRunAll.SimTNeuSelGood = [];
    spikeSimRunAll.SimTRecNoGood = [];
    spikeSimRunAll.SimTNonZero = [];
    spikeSimRunAll.SimT = [];
    spikeSimRunAll.SimTNeuSel = [];
    spikeSimRunAll.SimTRecNo = [];
        
    spikeSimRewAll.SimTNonZeroGood = [];
    spikeSimRewAll.SimTGood = [];
    spikeSimRewAll.SimTNonZero = [];
    spikeSimRewAll.SimT = [];
    
    spikeSimCueAll.SimTNonZeroGood = [];
    spikeSimCueAll.SimTGood = [];
    spikeSimCueAll.SimTNonZero = [];
    spikeSimCueAll.SimT = [];
    for i = 1:size(listRecordingsActiveLickPath,1)
        path = listRecordingsActiveLickPath(i,:);
        fileName = listRecordingsActiveLickFileName(i,:);
        mazeSess = mazeSessionActiveLick(i);
        fullPath = [path,fileName, '_meanSpikeTrainSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
        load(fullPath,'meanSimTRun','meanSimTRew','meanSimTCue');
        
        fileNameInfo = [fileName '_Info.mat'];
        fullPath = [path fileNameInfo];
        if(exist(fullPath) == 0)
            disp('_Info.mat file does not exist.');
            return;
        end
        load(fullPath,'autoCorr','beh'); 
        
        fullPathFR = [fileName '_FR_Run' num2str(onlyRun) '.mat'];
        fullPath = [path fullPathFR];
        if(exist(fullPath) == 0)
            disp('_FR_Run.mat file does not exist.');
            return;
        end
        load(fullPath,'mFRStruct','mFRStructSess'); 
        if(length(beh.mazeSessAll) > 1)
            mFR = mFRStructSess{mazeSess};
        else
            mFR = mFRStruct;
        end
        
        indNeuSel = mFR.mFR > minFR & mFR.mFR < maxFR &...
                    autoCorr.isPyrneuron == 1;
                
        indSelSimT = meanSimTRun.meanGoodNZ > thrCorrT;
        indSelSimT = indNeuSel & indSelSimT;
        spikeSimRunAll.SimTNeuSelGood = [spikeSimRunAll.SimTNeuSelGood ...
            find(indSelSimT == 1)];
        spikeSimRunAll.SimTRecNoGood = [spikeSimRunAll.SimTRecNoGood ...
            i*ones(1,sum(indSelSimT))];
        spikeSimRunAll.SimTNonZeroGood = [spikeSimRunAll.SimTNonZeroGood ...
            meanSimTRun.meanGoodNZ(indSelSimT)];
        spikeSimRewAll.SimTNonZeroGood = [spikeSimRewAll.SimTNonZeroGood ...
            meanSimTRew.meanGoodNZ(indSelSimT)];
        spikeSimCueAll.SimTNonZeroGood = [spikeSimCueAll.SimTNonZeroGood ...
            meanSimTCue.meanGoodNZ(indSelSimT)];
        
        spikeSimRunAll.SimTGood = [spikeSimRunAll.SimTGood ...
            meanSimTRun.meanGood(indSelSimT)];
        spikeSimRewAll.SimTGood = [spikeSimRewAll.SimTGood ...
            meanSimTRew.meanGood(indSelSimT)];
        spikeSimCueAll.SimTGood = [spikeSimCueAll.SimTGood ...
            meanSimTCue.meanGood(indSelSimT)];
               
        indSelSimT = meanSimTRun.meanNZ > thrCorrT;
        indSelSimT = indNeuSel & indSelSimT;
        spikeSimRunAll.SimTNeuSel = [spikeSimRunAll.SimTNeuSel ...
            find(indSelSimT == 1)];
        spikeSimRunAll.SimTRecNo = [spikeSimRunAll.SimTRecNo ...
            i*ones(1,sum(indSelSimT))];
        spikeSimRunAll.SimTNonZero = [spikeSimRunAll.SimTNonZero ...
            meanSimTRun.meanNZ(indSelSimT)];
        spikeSimRewAll.SimTNonZero = [spikeSimRewAll.SimTNonZero ...
            meanSimTRew.meanNZ(indSelSimT)];
        spikeSimCueAll.SimTNonZero = [spikeSimCueAll.SimTNonZero ...
            meanSimTCue.meanNZ(indSelSimT)];
        
        spikeSimRunAll.SimT = [spikeSimRunAll.SimT ...
            meanSimTRun.mean(indSelSimT)];
        spikeSimRewAll.SimT = [spikeSimRewAll.SimT ...
            meanSimTRew.mean(indSelSimT)];
        spikeSimCueAll.SimT = [spikeSimCueAll.SimT ...
            meanSimTCue.mean(indSelSimT)];
    end
    
    spikeSimRunAll.pRS_SimTNonZeroGood_RR = ranksum(spikeSimRunAll.SimTNonZeroGood,...
        spikeSimRewAll.SimTNonZeroGood);
    spikeSimRunAll.pRS_SimTNonZeroGood_RC = ranksum(spikeSimRunAll.SimTNonZeroGood,...
        spikeSimCueAll.SimTNonZeroGood);
    
    spikeSimRunAll.pRS_SimTGood_RR = ranksum(spikeSimRunAll.SimTGood,...
        spikeSimRewAll.SimTGood);
    spikeSimRunAll.pRS_SimTGood_RC = ranksum(spikeSimRunAll.SimTGood,...
        spikeSimCueAll.SimTGood);
    
    spikeSimRunAll.pRS_SimTNonZero_RR = ranksum(spikeSimRunAll.SimTNonZero,...
        spikeSimRewAll.SimTNonZero);
    spikeSimRunAll.pRS_SimTNonZero_RC = ranksum(spikeSimRunAll.SimTNonZero,...
        spikeSimCueAll.SimTNonZero);
    
    spikeSimRunAll.pRS_SimT_RR = ranksum(spikeSimRunAll.SimT,...
        spikeSimRewAll.SimT);
    spikeSimRunAll.pRS_SimT_RC = ranksum(spikeSimRunAll.SimT,...
        spikeSimCueAll.SimT);
    
    save('Z:\Yingxue\DataAnalysisRaphi\SpikeSimTAlignedAL.mat','spikeSimRunAll',...
        'spikeSimRewAll','spikeSimCueAll');
    
    plotSimComp(spikeSimRunAll.SimTNonZeroGood,...
        spikeSimRewAll.SimTNonZeroGood,...
        spikeSimRunAll.pRS_SimTNonZeroGood_RR,...
        'Spike simT NZgood aligned to run onset',...
        'Spike simT NZgood aligned to reward onset');
    fileName1 = ['SpikeSimTAligned_NZGood_RunVsRew_ALRec'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');

    plotSimComp(spikeSimRunAll.SimTNonZeroGood,...
        spikeSimCueAll.SimTNonZeroGood,...
        spikeSimRunAll.pRS_SimTNonZeroGood_RC,...
        'Spike simT NZgood aligned to run onset',...
        'Spike simT NZgood aligned to cue onset');
    fileName1 = ['SpikeSimTAligned_NZGood_RunVsCue_ALRec'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');
    
    plotSimComp(spikeSimRunAll.SimTGood,...
        spikeSimRewAll.SimTGood,...
        spikeSimRunAll.pRS_SimTGood_RR,...
        'Spike simT good aligned to run onset',...
        'Spike simT good aligned to reward onset');
    fileName1 = ['SpikeSimTAligned_Good_RunVsRew_ALRec'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');

    plotSimComp(spikeSimRunAll.SimTGood,...
        spikeSimCueAll.SimTGood,...
        spikeSimRunAll.pRS_SimTGood_RC,...
        'Spike simT good aligned to run onset',...
        'Spike simT good aligned to cue onset');
    fileName1 = ['SpikeSimTAligned_Good_RunVsCue_ALRec'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');
    
    plotSimComp(spikeSimRunAll.SimTNonZero,...
        spikeSimRewAll.SimTNonZero,...
        spikeSimRunAll.pRS_SimTNonZero_RR,...
        'Spike simT NZ aligned to run onset',...
        'Spike simT NZ aligned to reward onset');
    fileName1 = ['SpikeSimTAligned_NZ_RunVsRew_ALRec'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');

    plotSimComp(spikeSimRunAll.SimTNonZero,...
        spikeSimCueAll.SimTNonZero,...
        spikeSimRunAll.pRS_SimTNonZero_RC,...
        'Spike simT NZ aligned to run onset',...
        'Spike simT NZ aligned to cue onset');
    fileName1 = ['SpikeSimTAligned_NZ_RunVsCue_ALRec'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');
    
    plotSimComp(spikeSimRunAll.SimT,...
        spikeSimRewAll.SimT,...
        spikeSimRunAll.pRS_SimT_RR,...
        'Spike simT aligned to run onset',...
        'Spike simT aligned to reward onset');
    fileName1 = ['SpikeSimTAligned_RunVsRew_ALRec'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');

    plotSimComp(spikeSimRunAll.SimT,...
        spikeSimCueAll.SimT,...
        spikeSimRunAll.pRS_SimT_RC,...
        'Spike simT aligned to run onset',...
        'Spike simT aligned to cue onset');
    fileName1 = ['SpikeSimTAligned_RunVsCue_ALRec'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');
    
    colorSel= 0;
    plotBoxPlot(spikeSimRunAll.SimTGood,...
        spikeSimCueAll.SimTGood,'Spike corr. good',...
        'SpikeSimTAligned_Good_RunVsCue_ALRecBox',...
        'Z:\Yingxue\DataAnalysisRaphi\',[],spikeSimRunAll.pRS_SimTGood_RC,colorSel);
    
    plotBoxPlot(spikeSimRunAll.SimTNonZeroGood,...
        spikeSimCueAll.SimTNonZeroGood,'Spike corr. good(nonzero)',...
        'SpikeSimTAligned_NZGood_RunVsCue_ALRecBox',...
        'Z:\Yingxue\DataAnalysisRaphi\',[],spikeSimRunAll.pRS_SimTNonZeroGood_RC,colorSel);
    
    plotBoxPlot(spikeSimRunAll.SimT,...
        spikeSimCueAll.SimT,'Spike corr.',...
        'SpikeSimTAligned_RunVsCue_ALRecBox',...
        'Z:\Yingxue\DataAnalysisRaphi\',[],spikeSimRunAll.pRS_SimT_RC,colorSel);
    
    plotBoxPlot(spikeSimRunAll.SimTNonZero,...
        spikeSimCueAll.SimTNonZero,'Spike corr.(nonzero)',...
        'SpikeSimTAligned_NZ_RunVsCue_ALRecBox',...
        'Z:\Yingxue\DataAnalysisRaphi\',[],spikeSimRunAll.pRS_SimTNonZero_RC,colorSel);
end

function plotSimComp(x,y,p,xl,yl)
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 280 280])
    h = plot(x,y,'.');
    set(h,'MarkerSize',4,'Color',[0.5 0.5 0.9]);
    maxXY = max([x,y]);
    hold on;
    h = plot([0 maxXY],[0 maxXY],'r:');
    set(h,'LineWidth',1);
    set(gca,'XLim',[0 maxXY],'YLim',[0 maxXY]);
    title(['p = ' num2str(p)]);
    xlabel(xl)
    ylabel(yl)
end

function plotBoxPlot(x1,x2,yl,fn,pathAnal,ylimit,p,colorSel)
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 200 400])
    if(colorSel == 0)
        colorArr = [163 207 98;...
                234 131 114]/255;
    elseif(colorSel == 1)            
        colorArr = [234 131 114;...
                116 53 61]/255;
    else        
        colorArr = [163 207 98;... 
            63 79 37]/255;
    end
    x = [x1';x2'];
    g = [repmat({'C1'},length(x1),1);...
        repmat({'C2'},length(x2),1)];
    boxplot(x,g,'Notch','on','Widths',0.3,'Symbol','');
    h = findobj(gca,'Tag','Box');
    for j = 1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),colorArr(j,:),'FaceAlpha',0.5);
    end
    if(~isempty(ylimit))
        set(gca,'YLim',ylimit);
    end
    ylabel(yl);
    title(['p = ' num2str(p)]);
    
    fileName1 = [pathAnal fn];
    saveas(gcf,fileName1);
    print('-painters', '-dpdf', fileName1, '-r600')
end