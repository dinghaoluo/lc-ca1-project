function spikeTrainSimilarityAllRec(onlyRun)
% compare the spike train similarity for different alignment conditions over all the recordings 

    GlobalConst;
    intervalT = 20;
    intervalD = 1800;
    tc = 500; % ms
    cost = 2/(tc/1000*sampleFq);
    
    thrCorrT = 0; %0.03;
    
    RecordingList;
    
    spikeSimRunAll.CorrTNonZeroGood = [];
    spikeSimRunAll.CorrTGood = [];
    spikeSimRunAll.CorrTNeuSel = [];
    spikeSimRunAll.CorrTRecNo = [];
        
    spikeSimRewAll.CorrTNonZeroGood = [];
    spikeSimRewAll.CorrTGood = [];
    
    spikeSimCueAll.CorrTNonZeroGood = [];
    spikeSimCueAll.CorrTGood = [];
    for i = 1:size(listRecordingsActiveLickPath,1)
        path = listRecordingsActiveLickPath(i,:);
        fileName = listRecordingsActiveLickFileName(i,:);
        mazeSess = mazeSessionActiveLick(i);
        fullPath = [path,fileName, '_meanSpikesCorrTAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
        load(fullPath,'meanCorrTRun','meanCorrTRew','meanCorrTCue');
        
        indNeuSel = selPyrNeurons(path,fileName,onlyRun,minFR);
        indSelCorrT = meanCorrTRun.meanGoodNZ > thrCorrT;
        indSelCorrT = indNeuSel & indSelCorrT;
        spikeSimRunAll.CorrTNeuSel = [spikeSimRunAll.CorrTNeuSel ...
            find(indSelCorrT == 1)];
        spikeSimRunAll.CorrTRecNo = [spikeSimRunAll.CorrTRecNo ...
            i*ones(1,sum(indSelCorrT))];
        spikeSimRunAll.CorrTNonZeroGood = [spikeSimRunAll.CorrTNonZeroGood ...
            meanCorrTRun.meanGoodNZ(indSelCorrT)];
        spikeSimRewAll.CorrTNonZeroGood = [spikeSimRewAll.CorrTNonZeroGood ...
            meanCorrTRew.meanGoodNZ(indSelCorrT)];
        spikeSimCueAll.CorrTNonZeroGood = [spikeSimCueAll.CorrTNonZeroGood ...
            meanCorrTCue.meanGoodNZ(indSelCorrT)];
        
        spikeSimRunAll.CorrTGood = [spikeSimRunAll.CorrTGood ...
            meanCorrTRun.meanGood(indSelCorrT)];
        spikeSimRewAll.CorrTGood = [spikeSimRewAll.CorrTGood ...
            meanCorrTRew.meanGood(indSelCorrT)];
        spikeSimCueAll.CorrTGood = [spikeSimCueAll.CorrTGood ...
            meanCorrTCue.meanGood(indSelCorrT)];
        
%         path = listRecordingsActiveLickPath(i,:);
%         fileName = listRecordingsActiveLickFileName(i,:);
%         fullPath = [path,fileName, '_meanSpikesCorrDistAligned_Run' num2str(onlyRun) '_intD' ...
%             num2str(intervalT) '.mat'];
%         load(fullPath,'meanCorrDistRun','meanCorrDistRew','meanCorrDistCue');
%         
%         indSelCorrT = selNeuronsCorrD(path,fileName,onlyRun,intervalT,thrCorrD,minFR);
%         spikeSimRunAll.CorrTNeuSel = [spikeSimRunAll.CorrTNeuSel ...
%             find(indSelCorrT == 1)];
%         spikeSimRunAll.CorrTRecNo = [spikeSimRunAll.CorrTRecNo ...
%             i*ones(1,sum(indSelCorrT))];
%         spikeSimRunAll.CorrTNonZeroGood = [spikeSimRunAll.CorrTNonZeroGood ...
%             meanCorrTRun.meanGoodNZ(indSelCorrT)];
%         spikeSimRewAll.CorrTNonZeroGood = [spikeSimRewAll.CorrTNonZeroGood ...
%             meanCorrTRew.meanGoodNZ(indSelCorrT)];
%         spikeSimCueAll.CorrTNonZeroGood = [spikeSimCueAll.CorrTNonZeroGood ...
%             meanCorrTCue.meanGoodNZ(indSelCorrT)];
%         
%         spikeSimRunAll.CorrTGood = [spikeSimRunAll.CorrTGood ...
%             meanCorrTRun.meanGood(indSelCorrT)];
%         spikeSimRewAll.CorrTGood = [spikeSimRewAll.CorrTGood ...
%             meanCorrTRew.meanGood(indSelCorrT)];
%         spikeSimCueAll.CorrTGood = [spikeSimCueAll.CorrTGood ...
%             meanCorrTCue.meanGood(indSelCorrT)];
    end
    
    spikeSimRunAll.pRS_CorrTNonZeroGood_RR = ranksum(spikeSimRunAll.CorrTNonZeroGood,...
        spikeSimRewAll.CorrTNonZeroGood);
    spikeSimRunAll.pRS_CorrTNonZeroGood_RC = ranksum(spikeSimRunAll.CorrTNonZeroGood,...
        spikeSimCueAll.CorrTNonZeroGood);
    
    spikeSimRunAll.pRS_CorrTGood_RR = ranksum(spikeSimRunAll.CorrTGood,...
        spikeSimRewAll.CorrTGood);
    spikeSimRunAll.pRS_CorrTGood_RC = ranksum(spikeSimRunAll.CorrTGood,...
        spikeSimCueAll.CorrTGood);
    
    plotSimComp(spikeSimRunAll.CorrTNonZeroGood,...
        spikeSimRewAll.CorrTNonZeroGood,...
        spikeSimRunAll.pRS_CorrTNonZeroGood_RR,...
        'Spike corr. NZ aligned to run onset',...
        'Spike corr. NZ aligned to reward onset');
    fileName1 = ['SpikeCorrTAligned_NZGood_RunVsRew_AllRec'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');

    plotSimComp(spikeSimRunAll.CorrTNonZeroGood,...
        spikeSimCueAll.CorrTNonZeroGood,...
        spikeSimRunAll.pRS_CorrTNonZeroGood_RC,...
        'Spike corr. NZ aligned to run onset',...
        'Spike corr. NZ aligned to cue onset');
    fileName1 = ['SpikeCorrTAligned_NZGood_RunVsCue_AllRec'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');
    
    plotSimComp(spikeSimRunAll.CorrTGood,...
        spikeSimRewAll.CorrTGood,...
        spikeSimRunAll.pRS_CorrTGood_RR,...
        'Spike corr. aligned to run onset',...
        'Spike corr. aligned to reward onset');
    fileName1 = ['SpikeCorrTAligned_Good_RunVsRew_AllRec'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');

    plotSimComp(spikeSimRunAll.CorrTGood,...
        spikeSimCueAll.CorrTGood,...
        spikeSimRunAll.pRS_CorrTGood_RC,...
        'Spike corr. aligned to run onset',...
        'Spike corr. aligned to cue onset');
    fileName1 = ['SpikeCorrTAligned_Good_RunVsCue_AllRec'];
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