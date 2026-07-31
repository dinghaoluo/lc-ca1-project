function spikeTrainSimilarityTAllRecCompExp(onlyRun)
% compare the spike train similarity for different experiments 
% (passive no cue, passive licking, and active licking) over all the recordings 

    GlobalConst;
    intervalT = 20;
    intervalD = 1800;
    tc = 500; % ms
    cost = 2/(tc/1000*sampleFq);
    
    thrCorrT = -1;
    
    RecordingList;
    
    spikeSimPassiveNoCueAll.SimTNonZeroGoodRun = [];
    spikeSimPassiveNoCueAll.SimTGoodRun = [];
    spikeSimPassiveNoCueAll.SimTNonZeroGoodCue = [];
    spikeSimPassiveNoCueAll.SimTGoodCue = [];
    spikeSimPassiveNoCueAll.SimTNeuSel = [];
    spikeSimPassiveNoCueAll.SimTRecNo = [];
        
    spikeSimPassiveAll.SimTNonZeroGoodRun = [];
    spikeSimPassiveAll.SimTGoodRun = [];
    spikeSimPassiveAll.SimTNonZeroGoodCue = [];
    spikeSimPassiveAll.SimTGoodCue = [];
    spikeSimPassiveAll.SimTNeuSel = [];
    spikeSimPassiveAll.SimTRecNo = [];
    
    spikeSimActiveAll.SimTNonZeroGoodRun = [];
    spikeSimActiveAll.SimTGoodRun = [];
    spikeSimActiveAll.SimTNonZeroGoodCue = [];
    spikeSimActiveAll.SimTGoodCue = [];
    spikeSimActiveAll.SimTNeuSel = [];
    spikeSimActiveAll.SimTRecNo = [];
    
    for i = 1:size(listRecordingsNoCuePath,1)
        path = listRecordingsNoCuePath(i,:);
        fileName = listRecordingsNoCueFileName(i,:);
        mazeSess = mazeSessionNoCue(i);
        fullPath = [path,fileName, '_meanSpikeTrainSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
        load(fullPath,'meanSimTRun','meanSimTCue');
        
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
        spikeSimPassiveNoCueAll.SimTNeuSel = [spikeSimPassiveNoCueAll.SimTNeuSel ...
            find(indSelSimT == 1)];
        spikeSimPassiveNoCueAll.SimTRecNo = [spikeSimPassiveNoCueAll.SimTRecNo ...
            i*ones(1,sum(indSelSimT))];
        spikeSimPassiveNoCueAll.SimTNonZeroGoodRun = [spikeSimPassiveNoCueAll.SimTNonZeroGoodRun ...
            meanSimTRun.meanGoodNZ(indSelSimT)];
        spikeSimPassiveNoCueAll.SimTGoodRun = [spikeSimPassiveNoCueAll.SimTGoodRun ...
            meanSimTRun.meanGood(indSelSimT)];
        
       spikeSimPassiveNoCueAll.SimTNonZeroGoodCue = [spikeSimPassiveNoCueAll.SimTNonZeroGoodCue ...
            meanSimTCue.meanGoodNZ(indSelSimT)];
        spikeSimPassiveNoCueAll.SimTGoodCue = [spikeSimPassiveNoCueAll.SimTGoodCue ...
            meanSimTCue.meanGood(indSelSimT)];
    end
    
    for i = 1:size(listRecordingsPassiveLickPath,1)
        path = listRecordingsPassiveLickPath(i,:);
        fileName = listRecordingsPassiveLickFileName(i,:);
        mazeSess = mazeSessionPassiveLick(i);
        fullPath = [path,fileName, '_meanSpikeTrainSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
        load(fullPath,'meanSimTRun','meanSimTCue');
        
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
        spikeSimPassiveAll.SimTNeuSel = [spikeSimPassiveAll.SimTNeuSel ...
            find(indSelSimT == 1)];
        spikeSimPassiveAll.SimTRecNo = [spikeSimPassiveAll.SimTRecNo ...
            i*ones(1,sum(indSelSimT))];
        spikeSimPassiveAll.SimTNonZeroGoodRun = [spikeSimPassiveAll.SimTNonZeroGoodRun ...
            meanSimTRun.meanGoodNZ(indSelSimT)];
        spikeSimPassiveAll.SimTGoodRun = [spikeSimPassiveAll.SimTGoodRun ...
            meanSimTRun.meanGood(indSelSimT)];
        
        spikeSimPassiveAll.SimTNonZeroGoodCue = [spikeSimPassiveAll.SimTNonZeroGoodCue ...
            meanSimTCue.meanGoodNZ(indSelSimT)];
        spikeSimPassiveAll.SimTGoodCue = [spikeSimPassiveAll.SimTGoodCue ...
            meanSimTCue.meanGood(indSelSimT)];
    end
    
    for i = 1:size(listRecordingsActiveLickPath,1)
        path = listRecordingsActiveLickPath(i,:);
        fileName = listRecordingsActiveLickFileName(i,:);
        mazeSess = mazeSessionActiveLick(i);
        fullPath = [path,fileName, '_meanSpikeTrainSimT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
        load(fullPath,'meanSimTRun','meanSimTCue');
        
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
        spikeSimActiveAll.SimTNeuSel = [spikeSimActiveAll.SimTNeuSel ...
            find(indSelSimT == 1)];
        spikeSimActiveAll.SimTRecNo = [spikeSimActiveAll.SimTRecNo ...
            i*ones(1,sum(indSelSimT))];
        spikeSimActiveAll.SimTNonZeroGoodRun = [spikeSimActiveAll.SimTNonZeroGoodRun ...
            meanSimTRun.meanGoodNZ(indSelSimT)];
        spikeSimActiveAll.SimTGoodRun = [spikeSimActiveAll.SimTGoodRun ...
            meanSimTRun.meanGood(indSelSimT)];
        
        spikeSimActiveAll.SimTNonZeroGoodCue = [spikeSimActiveAll.SimTNonZeroGoodCue ...
            meanSimTCue.meanGoodNZ(indSelSimT)];
        spikeSimActiveAll.SimTGoodCue = [spikeSimActiveAll.SimTGoodCue ...
            meanSimTCue.meanGood(indSelSimT)];
    end
    
    %% statistics for aligning to run onset
    spikeSimRunAll.pRS_SimTNonZeroGood_PNCvsP = ranksum(spikeSimPassiveNoCueAll.SimTNonZeroGoodRun,...
        spikeSimPassiveAll.SimTNonZeroGoodRun);
    spikeSimRunAll.pRS_SimTNonZeroGood_PNCvsA = ranksum(spikeSimPassiveNoCueAll.SimTNonZeroGoodRun,...
        spikeSimActiveAll.SimTNonZeroGoodRun);
    spikeSimRunAll.pRS_SimTNonZeroGood_PvsA = ranksum(spikeSimPassiveAll.SimTNonZeroGoodRun,...
        spikeSimActiveAll.SimTNonZeroGoodRun);
    
    spikeSimRunAll.pRS_SimTGood_PNCvsP = ranksum(spikeSimPassiveNoCueAll.SimTGoodRun,...
        spikeSimPassiveAll.SimTGoodRun);
    spikeSimRunAll.pRS_SimTGood_PNCvsA = ranksum(spikeSimPassiveNoCueAll.SimTGoodRun,...
        spikeSimActiveAll.SimTGoodRun);
    spikeSimRunAll.pRS_SimTGood_PvsA = ranksum(spikeSimPassiveAll.SimTGoodRun,...
        spikeSimActiveAll.SimTGoodRun);
    
    [~,spikeSimRunAll.pKS_SimTNonZeroGood_PNCvsP] = kstest2(spikeSimPassiveNoCueAll.SimTNonZeroGoodRun,...
        spikeSimPassiveAll.SimTNonZeroGoodRun);
    [~,spikeSimRunAll.pKS_SimTNonZeroGood_PNCvsA] = kstest2(spikeSimPassiveNoCueAll.SimTNonZeroGoodRun,...
        spikeSimActiveAll.SimTNonZeroGoodRun);
    [~,spikeSimRunAll.pKS_SimTNonZeroGood_PvsA] = kstest2(spikeSimPassiveAll.SimTNonZeroGoodRun,...
        spikeSimActiveAll.SimTNonZeroGoodRun);
    
    [~,spikeSimRunAll.pKS_SimTGood_PNCvsP] = kstest2(spikeSimPassiveNoCueAll.SimTGoodRun,...
        spikeSimPassiveAll.SimTGoodRun);
    [~,spikeSimRunAll.pKS_SimTGood_PNCvsA] = kstest2(spikeSimPassiveNoCueAll.SimTGoodRun,...
        spikeSimActiveAll.SimTGoodRun);
    [~,spikeSimRunAll.pKS_SimTGood_PvsA] = kstest2(spikeSimPassiveAll.SimTGoodRun,...
        spikeSimActiveAll.SimTGoodRun);
    
    %% statistics for aligning to cue
    spikeSimCueAll.pRS_SimTNonZeroGood_PNCvsP = ranksum(spikeSimPassiveNoCueAll.SimTNonZeroGoodCue,...
        spikeSimPassiveAll.SimTNonZeroGoodCue);
    spikeSimCueAll.pRS_SimTNonZeroGood_PNCvsA = ranksum(spikeSimPassiveNoCueAll.SimTNonZeroGoodCue,...
        spikeSimActiveAll.SimTNonZeroGoodCue);
    spikeSimCueAll.pRS_SimTNonZeroGood_PvsA = ranksum(spikeSimPassiveAll.SimTNonZeroGoodCue,...
        spikeSimActiveAll.SimTNonZeroGoodCue);
    
    spikeSimCueAll.pRS_SimTGood_PNCvsP = ranksum(spikeSimPassiveNoCueAll.SimTGoodCue,...
        spikeSimPassiveAll.SimTGoodCue);
    spikeSimCueAll.pRS_SimTGood_PNCvsA = ranksum(spikeSimPassiveNoCueAll.SimTGoodCue,...
        spikeSimActiveAll.SimTGoodCue);
    spikeSimCueAll.pRS_SimTGood_PvsA = ranksum(spikeSimPassiveAll.SimTGoodCue,...
        spikeSimActiveAll.SimTGoodCue);
    
    [~,spikeSimCueAll.pKS_SimTNonZeroGood_PNCvsP] = kstest2(spikeSimPassiveNoCueAll.SimTNonZeroGoodCue,...
        spikeSimPassiveAll.SimTNonZeroGoodCue);
    [~,spikeSimCueAll.pKS_SimTNonZeroGood_PNCvsA] = kstest2(spikeSimPassiveNoCueAll.SimTNonZeroGoodCue,...
        spikeSimActiveAll.SimTNonZeroGoodCue);
    [~,spikeSimCueAll.pKS_SimTNonZeroGood_PvsA] = kstest2(spikeSimPassiveAll.SimTNonZeroGoodCue,...
        spikeSimActiveAll.SimTNonZeroGoodCue);
    
    [~,spikeSimCueAll.pKS_SimTGood_PNCvsP] = kstest2(spikeSimPassiveNoCueAll.SimTGoodCue,...
        spikeSimPassiveAll.SimTGoodCue);
    [~,spikeSimCueAll.pKS_SimTGood_PNCvsA] = kstest2(spikeSimPassiveNoCueAll.SimTGoodCue,...
        spikeSimActiveAll.SimTGoodCue);
    [~,spikeSimCueAll.pKS_SimTGood_PvsA] = kstest2(spikeSimPassiveAll.SimTGoodCue,...
        spikeSimActiveAll.SimTGoodCue);
    
    save('Z:\Yingxue\DataAnalysisRaphi\SpikeCorrTAligned_AllRec_AllCond.mat',...
        'spikeSimRunAll','spikeSimCueAll');
    
    %% plot the distribution aligned to run onset
    maxVal = max([spikeSimPassiveNoCueAll.SimTGoodRun,spikeSimPassiveAll.SimTGoodRun,spikeSimActiveAll.SimTGoodRun]);
    minVal = min([spikeSimPassiveNoCueAll.SimTGoodRun,spikeSimPassiveAll.SimTGoodRun,spikeSimActiveAll.SimTGoodRun]);
    xbins = minVal-0.01:0.01:maxVal+0.01;
    countPNC = hist(spikeSimPassiveNoCueAll.SimTGoodRun,xbins)/length(spikeSimPassiveNoCueAll.SimTGoodRun);
    countP = hist(spikeSimPassiveAll.SimTGoodRun,xbins)/length(spikeSimPassiveAll.SimTGoodRun);
    countA = hist(spikeSimActiveAll.SimTGoodRun,xbins)/length(spikeSimActiveAll.SimTGoodRun);
    
    plotSimComp(xbins,countPNC,countP,countA,...
        spikeSimRunAll.pKS_SimTGood_PNCvsP,spikeSimRunAll.pKS_SimTGood_PNCvsA,...
        spikeSimRunAll.pRS_SimTGood_PvsA,...
        'Spike simT. aligned to run onset',...
        'Count');
    fileName1 = ['SpikeSimTAlignedRun_Good_AllRec_AllCond_Distr'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');

    maxVal = max([spikeSimPassiveNoCueAll.SimTNonZeroGoodRun,spikeSimPassiveAll.SimTNonZeroGoodRun,spikeSimActiveAll.SimTNonZeroGoodRun]);
    minVal = min([spikeSimPassiveNoCueAll.SimTNonZeroGoodRun,spikeSimPassiveAll.SimTNonZeroGoodRun,spikeSimActiveAll.SimTNonZeroGoodRun]);
    xbins = minVal-0.01:0.01:maxVal+0.01;
    countPNC = hist(spikeSimPassiveNoCueAll.SimTNonZeroGoodRun,xbins)/length(spikeSimPassiveNoCueAll.SimTNonZeroGoodRun);
    countP = hist(spikeSimPassiveAll.SimTNonZeroGoodRun,xbins)/length(spikeSimPassiveAll.SimTNonZeroGoodRun);
    countA = hist(spikeSimActiveAll.SimTNonZeroGoodRun,xbins)/length(spikeSimActiveAll.SimTNonZeroGoodRun);
    
    plotSimComp(xbins,countPNC,countP,countA,...
        spikeSimRunAll.pKS_SimTNonZeroGood_PNCvsP,spikeSimRunAll.pKS_SimTNonZeroGood_PNCvsA,...
        spikeSimRunAll.pRS_SimTNonZeroGood_PvsA,...
        'Spike simT. NZ aligned to run onset',...
        'Count');
    fileName1 = ['SpikeSimTAlignedRun_NZGood_AllRec_AllCond_Distr'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');
    
    %% plot the distribution aligned to cue
    maxVal = max([spikeSimPassiveNoCueAll.SimTGoodCue,spikeSimPassiveAll.SimTGoodCue,spikeSimActiveAll.SimTGoodCue]);
    minVal = min([spikeSimPassiveNoCueAll.SimTGoodCue,spikeSimPassiveAll.SimTGoodCue,spikeSimActiveAll.SimTGoodCue]);
    xbins = minVal-0.01:0.01:maxVal+0.01;
    countPNC = hist(spikeSimPassiveNoCueAll.SimTGoodCue,xbins)/length(spikeSimPassiveNoCueAll.SimTGoodCue);
    countP = hist(spikeSimPassiveAll.SimTGoodCue,xbins)/length(spikeSimPassiveAll.SimTGoodCue);
    countA = hist(spikeSimActiveAll.SimTGoodCue,xbins)/length(spikeSimActiveAll.SimTGoodCue);
    
    plotSimComp(xbins,countPNC,countP,countA,...
        spikeSimCueAll.pKS_SimTGood_PNCvsP,spikeSimCueAll.pKS_SimTGood_PNCvsA,...
        spikeSimCueAll.pRS_SimTGood_PvsA,...
        'Spike simT. aligned to cue onset',...
        'Count');
    fileName1 = ['SpikeSimTAlignedCue_Good_AllRec_AllCond_Distr'];
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');

    maxVal = max([spikeSimPassiveNoCueAll.SimTNonZeroGoodCue,spikeSimPassiveAll.SimTNonZeroGoodCue,spikeSimActiveAll.SimTNonZeroGoodCue]);
    minVal = min([spikeSimPassiveNoCueAll.SimTNonZeroGoodCue,spikeSimPassiveAll.SimTNonZeroGoodCue,spikeSimActiveAll.SimTNonZeroGoodCue]);
    xbins = minVal-0.01:0.01:maxVal+0.01;
    countPNC = hist(spikeSimPassiveNoCueAll.SimTNonZeroGoodCue,xbins)/length(spikeSimPassiveNoCueAll.SimTNonZeroGoodCue);
    countP = hist(spikeSimPassiveAll.SimTNonZeroGoodCue,xbins)/length(spikeSimPassiveAll.SimTNonZeroGoodCue);
    countA = hist(spikeSimActiveAll.SimTNonZeroGoodCue,xbins)/length(spikeSimActiveAll.SimTNonZeroGoodCue);
    
    plotSimComp(xbins,countPNC,countP,countA,...
        spikeSimCueAll.pKS_SimTNonZeroGood_PNCvsP,spikeSimCueAll.pKS_SimTNonZeroGood_PNCvsA,...
        spikeSimCueAll.pRS_SimTNonZeroGood_PvsA,...
        'Spike simT. NZ aligned to cue onset',...
        'Count');
    fileName1 = ['SpikeSimTAlignedCue_NZGood_AllRec_AllCond_Distr'];
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