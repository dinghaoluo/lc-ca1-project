function spikeTrainSimilarityAllRecCompExp(onlyRun)
% compare the spike train similarity for different experiments 
% (passive no cue, passive licking, and active licking) over all the recordings 

    GlobalConst;
    intervalT = 20;
    intervalD = 1800;
    tc = 500; % ms
    cost = 2/(tc/1000*sampleFq);
    
    thrCorrT = 0;
    
    RecordingList;
    
    spikeSimPassiveNoCueAll.CorrTNonZeroGoodRun = [];
    spikeSimPassiveNoCueAll.CorrTGoodRun = [];
    spikeSimPassiveNoCueAll.CorrTNonZeroGoodCue = [];
    spikeSimPassiveNoCueAll.CorrTGoodCue = [];
    spikeSimPassiveNoCueAll.CorrTNeuSel = [];
    spikeSimPassiveNoCueAll.CorrTRecNo = [];
        
    spikeSimPassiveAll.CorrTNonZeroGoodRun = [];
    spikeSimPassiveAll.CorrTGoodRun = [];
    spikeSimPassiveAll.CorrTNonZeroGoodCue = [];
    spikeSimPassiveAll.CorrTGoodCue = [];
    spikeSimPassiveAll.CorrTNeuSel = [];
    spikeSimPassiveAll.CorrTRecNo = [];
    
    spikeSimActiveAll.CorrTNonZeroGoodRun = [];
    spikeSimActiveAll.CorrTGoodRun = [];
    spikeSimActiveAll.CorrTNonZeroGoodCue = [];
    spikeSimActiveAll.CorrTGoodCue = [];
    spikeSimActiveAll.CorrTNeuSel = [];
    spikeSimActiveAll.CorrTRecNo = [];
    
    for i = 1:size(listRecordingsNoCuePath,1)
        path = listRecordingsNoCuePath(i,:);
        fileName = listRecordingsNoCueFileName(i,:);
        mazeSess = mazeSessionNoCue(i);
        fullPath = [path,fileName, '_meanSpikesCorrTAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
        load(fullPath,'meanCorrTRun','meanCorrTCue');
        
        indNeuSel = selPyrNeurons(path,fileName,onlyRun,minFR);
        indSelCorrT = meanCorrTRun.meanGoodNZ > thrCorrT;
        indSelCorrT = indNeuSel & indSelCorrT;
        spikeSimPassiveNoCueAll.CorrTNeuSel = [spikeSimPassiveNoCueAll.CorrTNeuSel ...
            find(indSelCorrT == 1)];
        spikeSimPassiveNoCueAll.CorrTRecNo = [spikeSimPassiveNoCueAll.CorrTRecNo ...
            i*ones(1,sum(indSelCorrT))];
        spikeSimPassiveNoCueAll.CorrTNonZeroGoodRun = [spikeSimPassiveNoCueAll.CorrTNonZeroGoodRun ...
            meanCorrTRun.meanGoodNZ(indSelCorrT)];
        spikeSimPassiveNoCueAll.CorrTGoodRun = [spikeSimPassiveNoCueAll.CorrTGoodRun ...
            meanCorrTRun.meanGood(indSelCorrT)];
        
       spikeSimPassiveNoCueAll.CorrTNonZeroGoodCue = [spikeSimPassiveNoCueAll.CorrTNonZeroGoodCue ...
            meanCorrTCue.meanGoodNZ(indSelCorrT)];
        spikeSimPassiveNoCueAll.CorrTGoodCue = [spikeSimPassiveNoCueAll.CorrTGoodCue ...
            meanCorrTCue.meanGood(indSelCorrT)];
    end
    
    for i = 1:size(listRecordingsPassiveLickPath,1)
        path = listRecordingsPassiveLickPath(i,:);
        fileName = listRecordingsPassiveLickFileName(i,:);
        mazeSess = mazeSessionPassiveLick(i);
        fullPath = [path,fileName, '_meanSpikesCorrTAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
        load(fullPath,'meanCorrTRun','meanCorrTCue');
        
        indNeuSel = selPyrNeurons(path,fileName,onlyRun,minFR);
        indSelCorrT = meanCorrTRun.meanGoodNZ > thrCorrT;
        indSelCorrT = indNeuSel & indSelCorrT;
        spikeSimPassiveAll.CorrTNeuSel = [spikeSimPassiveAll.CorrTNeuSel ...
            find(indSelCorrT == 1)];
        spikeSimPassiveAll.CorrTRecNo = [spikeSimPassiveAll.CorrTRecNo ...
            i*ones(1,sum(indSelCorrT))];
        spikeSimPassiveAll.CorrTNonZeroGoodRun = [spikeSimPassiveAll.CorrTNonZeroGoodRun ...
            meanCorrTRun.meanGoodNZ(indSelCorrT)];
        spikeSimPassiveAll.CorrTGoodRun = [spikeSimPassiveAll.CorrTGoodRun ...
            meanCorrTRun.meanGood(indSelCorrT)];
        
        spikeSimPassiveAll.CorrTNonZeroGoodCue = [spikeSimPassiveAll.CorrTNonZeroGoodCue ...
            meanCorrTCue.meanGoodNZ(indSelCorrT)];
        spikeSimPassiveAll.CorrTGoodCue = [spikeSimPassiveAll.CorrTGoodCue ...
            meanCorrTCue.meanGood(indSelCorrT)];
    end
    
    for i = 1:size(listRecordingsActiveLickPath,1)
        path = listRecordingsActiveLickPath(i,:);
        fileName = listRecordingsActiveLickFileName(i,:);
        mazeSess = mazeSessionActiveLick(i);
        fullPath = [path,fileName, '_meanSpikesCorrTAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
        load(fullPath,'meanCorrTRun','meanCorrTCue');
        
        indNeuSel = selPyrNeurons(path,fileName,onlyRun,minFR);
        indSelCorrT = meanCorrTRun.meanGoodNZ > thrCorrT;
        indSelCorrT = indNeuSel & indSelCorrT;
        spikeSimActiveAll.CorrTNeuSel = [spikeSimActiveAll.CorrTNeuSel ...
            find(indSelCorrT == 1)];
        spikeSimActiveAll.CorrTRecNo = [spikeSimActiveAll.CorrTRecNo ...
            i*ones(1,sum(indSelCorrT))];
        spikeSimActiveAll.CorrTNonZeroGoodRun = [spikeSimActiveAll.CorrTNonZeroGoodRun ...
            meanCorrTRun.meanGoodNZ(indSelCorrT)];
        spikeSimActiveAll.CorrTGoodRun = [spikeSimActiveAll.CorrTGoodRun ...
            meanCorrTRun.meanGood(indSelCorrT)];
        
        spikeSimActiveAll.CorrTNonZeroGoodCue = [spikeSimActiveAll.CorrTNonZeroGoodCue ...
            meanCorrTCue.meanGoodNZ(indSelCorrT)];
        spikeSimActiveAll.CorrTGoodCue = [spikeSimActiveAll.CorrTGoodCue ...
            meanCorrTCue.meanGood(indSelCorrT)];
    end
    
    %% statistics for aligning to run onset
    spikeSimRunAll.pRS_CorrTNonZeroGood_PNCvsP = ranksum(spikeSimPassiveNoCueAll.CorrTNonZeroGoodRun,...
        spikeSimPassiveAll.CorrTNonZeroGoodRun);
    spikeSimRunAll.pRS_CorrTNonZeroGood_PNCvsA = ranksum(spikeSimPassiveNoCueAll.CorrTNonZeroGoodRun,...
        spikeSimActiveAll.CorrTNonZeroGoodRun);
    spikeSimRunAll.pRS_CorrTNonZeroGood_PvsA = ranksum(spikeSimPassiveAll.CorrTNonZeroGoodRun,...
        spikeSimActiveAll.CorrTNonZeroGoodRun);
    
    spikeSimRunAll.pRS_CorrTGood_PNCvsP = ranksum(spikeSimPassiveNoCueAll.CorrTGoodRun,...
        spikeSimPassiveAll.CorrTGoodRun);
    spikeSimRunAll.pRS_CorrTGood_PNCvsA = ranksum(spikeSimPassiveNoCueAll.CorrTGoodRun,...
        spikeSimActiveAll.CorrTGoodRun);
    spikeSimRunAll.pRS_CorrTGood_PvsA = ranksum(spikeSimPassiveAll.CorrTGoodRun,...
        spikeSimActiveAll.CorrTGoodRun);
    
    [~,spikeSimRunAll.pKS_CorrTNonZeroGood_PNCvsP] = kstest2(spikeSimPassiveNoCueAll.CorrTNonZeroGoodRun,...
        spikeSimPassiveAll.CorrTNonZeroGoodRun);
    [~,spikeSimRunAll.pKS_CorrTNonZeroGood_PNCvsA] = kstest2(spikeSimPassiveNoCueAll.CorrTNonZeroGoodRun,...
        spikeSimActiveAll.CorrTNonZeroGoodRun);
    [~,spikeSimRunAll.pKS_CorrTNonZeroGood_PvsA] = kstest2(spikeSimPassiveAll.CorrTNonZeroGoodRun,...
        spikeSimActiveAll.CorrTNonZeroGoodRun);
    
    [~,spikeSimRunAll.pKS_CorrTGood_PNCvsP] = kstest2(spikeSimPassiveNoCueAll.CorrTGoodRun,...
        spikeSimPassiveAll.CorrTGoodRun);
    [~,spikeSimRunAll.pKS_CorrTGood_PNCvsA] = kstest2(spikeSimPassiveNoCueAll.CorrTGoodRun,...
        spikeSimActiveAll.CorrTGoodRun);
    [~,spikeSimRunAll.pKS_CorrTGood_PvsA] = kstest2(spikeSimPassiveAll.CorrTGoodRun,...
        spikeSimActiveAll.CorrTGoodRun);
    
    %% statistics for aligning to cue
    spikeSimCueAll.pRS_CorrTNonZeroGood_PNCvsP = ranksum(spikeSimPassiveNoCueAll.CorrTNonZeroGoodCue,...
        spikeSimPassiveAll.CorrTNonZeroGoodCue);
    spikeSimCueAll.pRS_CorrTNonZeroGood_PNCvsA = ranksum(spikeSimPassiveNoCueAll.CorrTNonZeroGoodCue,...
        spikeSimActiveAll.CorrTNonZeroGoodCue);
    spikeSimCueAll.pRS_CorrTNonZeroGood_PvsA = ranksum(spikeSimPassiveAll.CorrTNonZeroGoodCue,...
        spikeSimActiveAll.CorrTNonZeroGoodCue);
    
    spikeSimCueAll.pRS_CorrTGood_PNCvsP = ranksum(spikeSimPassiveNoCueAll.CorrTGoodCue,...
        spikeSimPassiveAll.CorrTGoodCue);
    spikeSimCueAll.pRS_CorrTGood_PNCvsA = ranksum(spikeSimPassiveNoCueAll.CorrTGoodCue,...
        spikeSimActiveAll.CorrTGoodCue);
    spikeSimCueAll.pRS_CorrTGood_PvsA = ranksum(spikeSimPassiveAll.CorrTGoodCue,...
        spikeSimActiveAll.CorrTGoodCue);
    
    [~,spikeSimCueAll.pKS_CorrTNonZeroGood_PNCvsP] = kstest2(spikeSimPassiveNoCueAll.CorrTNonZeroGoodCue,...
        spikeSimPassiveAll.CorrTNonZeroGoodCue);
    [~,spikeSimCueAll.pKS_CorrTNonZeroGood_PNCvsA] = kstest2(spikeSimPassiveNoCueAll.CorrTNonZeroGoodCue,...
        spikeSimActiveAll.CorrTNonZeroGoodCue);
    [~,spikeSimCueAll.pKS_CorrTNonZeroGood_PvsA] = kstest2(spikeSimPassiveAll.CorrTNonZeroGoodCue,...
        spikeSimActiveAll.CorrTNonZeroGoodCue);
    
    [~,spikeSimCueAll.pKS_CorrTGood_PNCvsP] = kstest2(spikeSimPassiveNoCueAll.CorrTGoodCue,...
        spikeSimPassiveAll.CorrTGoodCue);
    [~,spikeSimCueAll.pKS_CorrTGood_PNCvsA] = kstest2(spikeSimPassiveNoCueAll.CorrTGoodCue,...
        spikeSimActiveAll.CorrTGoodCue);
    [~,spikeSimCueAll.pKS_CorrTGood_PvsA] = kstest2(spikeSimPassiveAll.CorrTGoodCue,...
        spikeSimActiveAll.CorrTGoodCue);
    
    
    %% plot the distribution aligned to run onset
    maxVal = max([spikeSimPassiveNoCueAll.CorrTGoodRun,spikeSimPassiveAll.CorrTGoodRun,spikeSimActiveAll.CorrTGoodRun]);
    minVal = min([spikeSimPassiveNoCueAll.CorrTGoodRun,spikeSimPassiveAll.CorrTGoodRun,spikeSimActiveAll.CorrTGoodRun]);
    xbins = minVal-0.01:0.01:maxVal+0.01;
    countPNC = hist(spikeSimPassiveNoCueAll.CorrTGoodRun,xbins)/length(spikeSimPassiveNoCueAll.CorrTGoodRun);
    countP = hist(spikeSimPassiveAll.CorrTGoodRun,xbins)/length(spikeSimPassiveAll.CorrTGoodRun);
    countA = hist(spikeSimActiveAll.CorrTGoodRun,xbins)/length(spikeSimActiveAll.CorrTGoodRun);
    
    plotSimComp(xbins,countPNC,countP,countA,...
        spikeSimRunAll.pKS_CorrTGood_PNCvsP,spikeSimRunAll.pKS_CorrTGood_PNCvsA,...
        spikeSimRunAll.pRS_CorrTGood_PvsA,...
        'Spike CorrT. aligned to run onset',...
        'Count');
    fileName1 = ['SpikeCorrTAlignedRun_Good_AllRec_AllCond_Distr'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');

    maxVal = max([spikeSimPassiveNoCueAll.CorrTNonZeroGoodRun,spikeSimPassiveAll.CorrTNonZeroGoodRun,spikeSimActiveAll.CorrTNonZeroGoodRun]);
    minVal = min([spikeSimPassiveNoCueAll.CorrTNonZeroGoodRun,spikeSimPassiveAll.CorrTNonZeroGoodRun,spikeSimActiveAll.CorrTNonZeroGoodRun]);
    xbins = minVal-0.01:0.01:maxVal+0.01;
    countPNC = hist(spikeSimPassiveNoCueAll.CorrTNonZeroGoodRun,xbins)/length(spikeSimPassiveNoCueAll.CorrTNonZeroGoodRun);
    countP = hist(spikeSimPassiveAll.CorrTNonZeroGoodRun,xbins)/length(spikeSimPassiveAll.CorrTNonZeroGoodRun);
    countA = hist(spikeSimActiveAll.CorrTNonZeroGoodRun,xbins)/length(spikeSimActiveAll.CorrTNonZeroGoodRun);
    
    plotSimComp(xbins,countPNC,countP,countA,...
        spikeSimRunAll.pKS_CorrTNonZeroGood_PNCvsP,spikeSimRunAll.pKS_CorrTNonZeroGood_PNCvsA,...
        spikeSimRunAll.pRS_CorrTNonZeroGood_PvsA,...
        'Spike CorrT. NZ aligned to run onset',...
        'Count');
    fileName1 = ['SpikeCorrTAlignedRun_NZGood_AllRec_AllCond_Distr'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');
    
    %% plot the distribution aligned to cue
    maxVal = max([spikeSimPassiveNoCueAll.CorrTGoodCue,spikeSimPassiveAll.CorrTGoodCue,spikeSimActiveAll.CorrTGoodCue]);
    minVal = min([spikeSimPassiveNoCueAll.CorrTGoodCue,spikeSimPassiveAll.CorrTGoodCue,spikeSimActiveAll.CorrTGoodCue]);
    xbins = minVal-0.01:0.01:maxVal+0.01;
    countPNC = hist(spikeSimPassiveNoCueAll.CorrTGoodCue,xbins)/length(spikeSimPassiveNoCueAll.CorrTGoodCue);
    countP = hist(spikeSimPassiveAll.CorrTGoodCue,xbins)/length(spikeSimPassiveAll.CorrTGoodCue);
    countA = hist(spikeSimActiveAll.CorrTGoodCue,xbins)/length(spikeSimActiveAll.CorrTGoodCue);
    
    plotSimComp(xbins,countPNC,countP,countA,...
        spikeSimCueAll.pKS_CorrTGood_PNCvsP,spikeSimCueAll.pKS_CorrTGood_PNCvsA,...
        spikeSimCueAll.pRS_CorrTGood_PvsA,...
        'Spike CorrT. aligned to cue onset',...
        'Count');
    fileName1 = ['SpikeCorrTAlignedCue_Good_AllRec_AllCond_Distr'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');

    maxVal = max([spikeSimPassiveNoCueAll.CorrTNonZeroGoodCue,spikeSimPassiveAll.CorrTNonZeroGoodCue,spikeSimActiveAll.CorrTNonZeroGoodCue]);
    minVal = min([spikeSimPassiveNoCueAll.CorrTNonZeroGoodCue,spikeSimPassiveAll.CorrTNonZeroGoodCue,spikeSimActiveAll.CorrTNonZeroGoodCue]);
    xbins = minVal-0.01:0.01:maxVal+0.01;
    countPNC = hist(spikeSimPassiveNoCueAll.CorrTNonZeroGoodCue,xbins)/length(spikeSimPassiveNoCueAll.CorrTNonZeroGoodCue);
    countP = hist(spikeSimPassiveAll.CorrTNonZeroGoodCue,xbins)/length(spikeSimPassiveAll.CorrTNonZeroGoodCue);
    countA = hist(spikeSimActiveAll.CorrTNonZeroGoodCue,xbins)/length(spikeSimActiveAll.CorrTNonZeroGoodCue);
    
    plotSimComp(xbins,countPNC,countP,countA,...
        spikeSimCueAll.pKS_CorrTNonZeroGood_PNCvsP,spikeSimCueAll.pKS_CorrTNonZeroGood_PNCvsA,...
        spikeSimCueAll.pRS_CorrTNonZeroGood_PvsA,...
        'Spike CorrT. NZ aligned to cue onset',...
        'Count');
    fileName1 = ['SpikeCorrTAlignedCue_NZGood_AllRec_AllCond_Distr'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');
end

function plotSimComp(a,x,y,z,pxy,pxz,pyz,xl,yl)
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 280 280])
    h = plot(a,x,'-');
    set(h,'LineWidth',2,'Color',[0.5 0.5 0.9]);
    hold on;
    h = plot(a,y,'-');
    set(h,'LineWidth',2,'Color',[0.9 0.5 0.5]);
    h = plot(a,z,'-');
    set(h,'LineWidth',2,'Color',[0.5 0.9 0.5]);
    set(gca,'XLim',[min(a) max(a)],'YLim',[0 max([x,y,z]+0.1)]);
    title(['p12 = ' num2str(pxy) ',p13 = ' num2str(pxz) ',p23 = ' num2str(pyz)]);
    xlabel(xl)
    ylabel(yl)
end