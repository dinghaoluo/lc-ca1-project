function spikeTrainSimilarityTAllRec(onlyRun)
% compare the spike train similarity for different alignment conditions over all the recordings 

    GlobalConst;
    intervalT = 20;
    intervalD = 1800;
    tc = 500; % ms
    cost = 2/(tc/1000*sampleFq);
    
    thrCorrT = 0; %0.03;
    
    RecordingList;
    
    spikeSimRunAll.SimTNonZeroGood = [];
    spikeSimRunAll.SimTGood = [];
    spikeSimRunAll.SimTNeuSel = [];
    spikeSimRunAll.SimTRecNo = [];
        
    spikeSimRewAll.SimTNonZeroGood = [];
    spikeSimRewAll.SimTGood = [];
    
    spikeSimCueAll.SimTNonZeroGood = [];
    spikeSimCueAll.SimTGood = [];
    for i = 1:size(listRecordingsActiveLickPath,1)
        path = listRecordingsActiveLickPath(i,:);
        fileName = listRecordingsActiveLickFileName(i,:);
        mazeSess = mazeSessionActiveLick(i);
        fullPath = [path,fileName, '_meanSpikeTrainSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
        load(fullPath,'meanSimTRun','meanSimTRew','meanSimTCue');
        
        indNeuSel = selPyrNeurons(path,fileName,onlyRun,minFR);
        indSelSimT = meanSimTRun.meanGoodNZ > thrCorrT;
        indSelSimT = indNeuSel & indSelSimT;
        spikeSimRunAll.SimTNeuSel = [spikeSimRunAll.SimTNeuSel ...
            find(indSelSimT == 1)];
        spikeSimRunAll.SimTRecNo = [spikeSimRunAll.SimTRecNo ...
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
        
%         path = listRecordingsActiveLickPath(i,:);
%         fileName = listRecordingsActiveLickFileName(i,:);
%         fullPath = [path,fileName, '_meanSpikesCorrDistAligned_Run' num2str(onlyRun) '_intD' ...
%             num2str(intervalT) '.mat'];
%         load(fullPath,'meanCorrDistRun','meanCorrDistRew','meanCorrDistCue');
%         
%         indSelSimT = selNeuronsCorrD(path,fileName,onlyRun,intervalT,thrCorrD,minFR);
%         spikeSimRunAll.SimTNeuSel = [spikeSimRunAll.SimTNeuSel ...
%             find(indSelSimT == 1)];
%         spikeSimRunAll.SimTRecNo = [spikeSimRunAll.SimTRecNo ...
%             i*ones(1,sum(indSelSimT))];
%         spikeSimRunAll.SimTNonZeroGood = [spikeSimRunAll.SimTNonZeroGood ...
%             meanSimTRun.meanGoodNZ(indSelSimT)];
%         spikeSimRewAll.SimTNonZeroGood = [spikeSimRewAll.SimTNonZeroGood ...
%             meanSimTRew.meanGoodNZ(indSelSimT)];
%         spikeSimCueAll.SimTNonZeroGood = [spikeSimCueAll.SimTNonZeroGood ...
%             meanSimTCue.meanGoodNZ(indSelSimT)];
%         
%         spikeSimRunAll.SimTGood = [spikeSimRunAll.SimTGood ...
%             meanSimTRun.meanGood(indSelSimT)];
%         spikeSimRewAll.SimTGood = [spikeSimRewAll.SimTGood ...
%             meanSimTRew.meanGood(indSelSimT)];
%         spikeSimCueAll.SimTGood = [spikeSimCueAll.SimTGood ...
%             meanSimTCue.meanGood(indSelSimT)];
    end
    
    spikeSimRunAll.pRS_SimTNonZeroGood_RR = ranksum(spikeSimRunAll.SimTNonZeroGood,...
        spikeSimRewAll.SimTNonZeroGood);
    spikeSimRunAll.pRS_SimTNonZeroGood_RC = ranksum(spikeSimRunAll.SimTNonZeroGood,...
        spikeSimCueAll.SimTNonZeroGood);
    
    spikeSimRunAll.pRS_SimTGood_RR = ranksum(spikeSimRunAll.SimTGood,...
        spikeSimRewAll.SimTGood);
    spikeSimRunAll.pRS_SimTGood_RC = ranksum(spikeSimRunAll.SimTGood,...
        spikeSimCueAll.SimTGood);
    
    plotSimComp(spikeSimRunAll.SimTNonZeroGood,...
        spikeSimRewAll.SimTNonZeroGood,...
        spikeSimRunAll.pRS_SimTNonZeroGood_RR,...
        'Spike simT NZ aligned to run onset',...
        'Spike simT NZ aligned to reward onset');
    fileName1 = ['SpikeSimTAligned_NZGood_RunVsRew_AllRec'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');

    plotSimComp(spikeSimRunAll.SimTNonZeroGood,...
        spikeSimCueAll.SimTNonZeroGood,...
        spikeSimRunAll.pRS_SimTNonZeroGood_RC,...
        'Spike simT NZ aligned to run onset',...
        'Spike simT NZ aligned to cue onset');
    fileName1 = ['SpikeSimTAligned_NZGood_RunVsCue_AllRec'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');
    
    plotSimComp(spikeSimRunAll.SimTGood,...
        spikeSimRewAll.SimTGood,...
        spikeSimRunAll.pRS_SimTGood_RR,...
        'Spike simT aligned to run onset',...
        'Spike simT aligned to reward onset');
    fileName1 = ['SpikeSimTAligned_Good_RunVsRew_AllRec'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');

    plotSimComp(spikeSimRunAll.SimTGood,...
        spikeSimCueAll.SimTGood,...
        spikeSimRunAll.pRS_SimTGood_RC,...
        'Spike simT aligned to run onset',...
        'Spike simT aligned to cue onset');
    fileName1 = ['SpikeSimTAligned_Good_RunVsCue_AllRec'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');
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