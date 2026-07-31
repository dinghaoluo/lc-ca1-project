function spikeTrainSimilarityVPAllRec(onlyRun)
% compare the spike train similarity for different alignment conditions over all the recordings 

    GlobalConst;
    intervalT = 20;
    intervalD = 1800;
    tc = 500; % ms
    cost = 2/(tc/1000*sampleFq);
    
    s = num2str(cost);
    ind = strfind(s,'.');
    s(ind) = 'p';
    
    thrCorrT = 0.03;
    
    RecordingList;
    
    spikeSimRunAll.SimVPNonZeroGood = [];
    spikeSimRunAll.SimVPGood = [];
    spikeSimRunAll.SimVPNeuSel = [];
    spikeSimRunAll.SimVPRecNo = [];
        
    spikeSimRewAll.SimVPNonZeroGood = [];
    spikeSimRewAll.SimVPGood = [];
    
    spikeSimCueAll.SimVPNonZeroGood = [];
    spikeSimCueAll.SimVPGood = [];
    for i = 1:size(listRecordingsActiveLickPath,1)
        path = listRecordingsActiveLickPath(i,:);
        fileName = listRecordingsActiveLickFileName(i,:);
        mazeSess = mazeSessionActiveLick(i);
        fullPath = [path,fileName, '_meanSpikeTrainSimVP_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_q' s '_intT' ...
            num2str(intervalT) '.mat'];
        load(fullPath,'meanSimVPRun','meanSimVPRew','meanSimVPCue');
        
        indNeuSel = selPyrNeurons(path,fileName,onlyRun,minFR);
        indSelSimVP = meanSimVPRun.meanGoodNZ > thrCorrT;
        indSelSimVP = indNeuSel & indSelSimVP;
        spikeSimRunAll.SimVPNeuSel = [spikeSimRunAll.SimVPNeuSel ...
            find(indSelSimVP == 1)];
        spikeSimRunAll.SimVPRecNo = [spikeSimRunAll.SimVPRecNo ...
            i*ones(1,sum(indSelSimVP))];
        spikeSimRunAll.SimVPNonZeroGood = [spikeSimRunAll.SimVPNonZeroGood ...
            meanSimVPRun.meanGoodNZ(indSelSimVP)];
        spikeSimRewAll.SimVPNonZeroGood = [spikeSimRewAll.SimVPNonZeroGood ...
            meanSimVPRew.meanGoodNZ(indSelSimVP)];
        spikeSimCueAll.SimVPNonZeroGood = [spikeSimCueAll.SimVPNonZeroGood ...
            meanSimVPCue.meanGoodNZ(indSelSimVP)];
        
        spikeSimRunAll.SimVPGood = [spikeSimRunAll.SimVPGood ...
            meanSimVPRun.meanGood(indSelSimVP)];
        spikeSimRewAll.SimVPGood = [spikeSimRewAll.SimVPGood ...
            meanSimVPRew.meanGood(indSelSimVP)];
        spikeSimCueAll.SimVPGood = [spikeSimCueAll.SimVPGood ...
            meanSimVPCue.meanGood(indSelSimVP)];
        
%         path = listRecordingsActiveLickPath(i,:);
%         fileName = listRecordingsActiveLickFileName(i,:);
%         fullPath = [path,fileName, '_meanSpikesCorrDistAligned_Run' num2str(onlyRun) '_intD' ...
%             num2str(intervalT) '.mat'];
%         load(fullPath,'meanCorrDistRun','meanCorrDistRew','meanCorrDistCue');
%         
%         indSelSimVP = selNeuronsCorrD(path,fileName,onlyRun,intervalT,thrCorrD,minFR);
%         spikeSimRunAll.SimVPNeuSel = [spikeSimRunAll.SimVPNeuSel ...
%             find(indSelSimVP == 1)];
%         spikeSimRunAll.SimVPRecNo = [spikeSimRunAll.SimVPRecNo ...
%             i*ones(1,sum(indSelSimVP))];
%         spikeSimRunAll.SimVPNonZeroGood = [spikeSimRunAll.SimVPNonZeroGood ...
%             meanSimVPRun.meanGoodNZ(indSelSimVP)];
%         spikeSimRewAll.SimVPNonZeroGood = [spikeSimRewAll.SimVPNonZeroGood ...
%             meanSimVPRew.meanGoodNZ(indSelSimVP)];
%         spikeSimCueAll.SimVPNonZeroGood = [spikeSimCueAll.SimVPNonZeroGood ...
%             meanSimVPCue.meanGoodNZ(indSelSimVP)];
%         
%         spikeSimRunAll.SimVPGood = [spikeSimRunAll.SimVPGood ...
%             meanSimVPRun.meanGood(indSelSimVP)];
%         spikeSimRewAll.SimVPGood = [spikeSimRewAll.SimVPGood ...
%             meanSimVPRew.meanGood(indSelSimVP)];
%         spikeSimCueAll.SimVPGood = [spikeSimCueAll.SimVPGood ...
%             meanSimVPCue.meanGood(indSelSimVP)];
    end
    
    spikeSimRunAll.pRS_SimVPNonZeroGood_RR = ranksum(spikeSimRunAll.SimVPNonZeroGood,...
        spikeSimRewAll.SimVPNonZeroGood);
    spikeSimRunAll.pRS_SimVPNonZeroGood_RC = ranksum(spikeSimRunAll.SimVPNonZeroGood,...
        spikeSimCueAll.SimVPNonZeroGood);
    
    spikeSimRunAll.pRS_SimVPGood_RR = ranksum(spikeSimRunAll.SimVPGood,...
        spikeSimRewAll.SimVPGood);
    spikeSimRunAll.pRS_SimVPGood_RC = ranksum(spikeSimRunAll.SimVPGood,...
        spikeSimCueAll.SimVPGood);
    
    plotSimComp(spikeSimRunAll.SimVPNonZeroGood,...
        spikeSimRewAll.SimVPNonZeroGood,...
        spikeSimRunAll.pRS_SimVPNonZeroGood_RR,...
        'Spike Dist VP NZ aligned to run onset',...
        'Spike Dist VP NZ aligned to reward onset');
    fileName1 = ['SpikeDistVPAligned_NZGood_RunVsRew_AllRec'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');

    plotSimComp(spikeSimRunAll.SimVPNonZeroGood,...
        spikeSimCueAll.SimVPNonZeroGood,...
        spikeSimRunAll.pRS_SimVPNonZeroGood_RC,...
        'Spike Dist VP NZ aligned to run onset',...
        'Spike Dist VP NZ aligned to cue onset');
    fileName1 = ['SpikeDistVPAligned_NZGood_RunVsCue_AllRec'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');
    
    plotSimComp(spikeSimRunAll.SimVPGood,...
        spikeSimRewAll.SimVPGood,...
        spikeSimRunAll.pRS_SimVPGood_RR,...
        'Spike Dist VP aligned to run onset',...
        'Spike Dist VP aligned to reward onset');
    fileName1 = ['SpikeDistVPAligned_Good_RunVsRew_AllRec'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');

    plotSimComp(spikeSimRunAll.SimVPGood,...
        spikeSimCueAll.SimVPGood,...
        spikeSimRunAll.pRS_SimVPGood_RC,...
        'Spike Dist VP aligned to run onset',...
        'Spike Dist VP aligned to cue onset');
    fileName1 = ['SpikeDistVPAligned_Good_RunVsCue_AllRec'];
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