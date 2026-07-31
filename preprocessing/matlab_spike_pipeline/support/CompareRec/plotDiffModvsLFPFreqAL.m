function plotDiffModvsLFPFreqAL(pathAnal,type1,type2)
% compare the passive no cue recordings with the active licking no cue
% behavior, and look at the population activity correlation over trials
% default analysis path ['Z:\Raphael_tests\mice_expdata\Analysis\'];
% e.g.  plotDiffModvsLFPFreqAL('Z:\Raphael_tests\mice_expdata\Analysis\',1,3,2)
% type1/type2 = 1, passive licking; 2, active licking; 3, passive no cue

    load([pathAnal 'popActivityCorr_RecList' num2str(type1) '_Run1.mat']);
    recList1 = recListPopActivityCorr;
    
    load([pathAnal 'popActivityCorr_RecList' num2str(type2) '_Run1.mat']);
    recList2 = recListPopActivityCorr;
    
    RecordingList;
    pathAnal1 = 'Z:\Yingxue\DataAnalysisRaphi\Pyramidal\';
    
    if(exist([pathAnal1 'autoCorrPyrAllRec.mat']))
        load([pathAnal1 'autoCorrPyrAllRec.mat'],...
            'autoCorrPyrNoCue','autoCorrPyrAL','autoCorrPyrPL',...
            'modPyrNoCue','modPyrAL','modPyrPL');
    end
    
    diffNeuronLFPFreqNoCueField = [];
    diffNeuronLFPFreqNoCueNoField = [];
    for i = 1:length(recList1.indPyr)
        indRec = autoCorrPyrNoCue.indRec == i;
        [~,indNeu] = intersect(autoCorrPyrNoCue.indNeu(indRec),recList1.indPyr{i});
        diffNeuronLFPFreqNoCueField = [diffNeuronLFPFreqNoCueField ...
            modPyrNoCue.thetaModFreq3(indNeu)-modPyrNoCue.thetaFreqHMean(indNeu)];
        [~,indNeu] = setdiff(autoCorrPyrNoCue.indNeu(indRec),recList1.indPyr{i});
        diffNeuronLFPFreqNoCueNoField = [diffNeuronLFPFreqNoCueNoField ...
            modPyrNoCue.thetaModFreq3(indNeu)-modPyrNoCue.thetaFreqHMean(indNeu)];
    end
    
    diffNeuronLFPFreqALField = [];
    diffNeuronLFPFreqALNoField = [];
    for i = 1:length(recList2.indPyr)
        indRec = autoCorrPyrAL.indRec == i;
        [~,indNeu] = intersect(autoCorrPyrAL.indNeu(indRec),recList2.indPyr{i});
        diffNeuronLFPFreqALField = [diffNeuronLFPFreqALField ...
            modPyrAL.thetaModFreq3(indNeu)-modPyrAL.thetaFreqHMean(indNeu)];
        [~,indNeu] = setdiff(autoCorrPyrAL.indNeu(indRec),recList2.indPyr{i});
        diffNeuronLFPFreqALNoField = [diffNeuronLFPFreqALNoField ...
            modPyrAL.thetaModFreq3(indNeu)-modPyrAL.thetaFreqHMean(indNeu)];
    end
    
    meanDiffNeuronLFPFreqNoCueField = mean(diffNeuronLFPFreqNoCueField);
    meanDiffNeuronLFPFreqNoCueNoField = mean(diffNeuronLFPFreqNoCueNoField);
    semDiffNeuronLFPFreqNoCueField = std(diffNeuronLFPFreqNoCueField)/sqrt(length(diffNeuronLFPFreqNoCueField));
    semDiffNeuronLFPFreqNoCueNoField = std(diffNeuronLFPFreqNoCueNoField)/sqrt(length(diffNeuronLFPFreqNoCueNoField));
    pRSDiffNeuronLFPNoCueFreq = ranksum(diffNeuronLFPFreqNoCueField,diffNeuronLFPFreqNoCueNoField);
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels') 
    figTitle = '';
    set(figure(figNew),'OuterPosition',...
        [pos(1) pos(2) 280 280],'Name',figTitle)
   
    plotStat(meanDiffNeuronLFPFreqNoCueNoField,meanDiffNeuronLFPFreqNoCueField,...
        semDiffNeuronLFPFreqNoCueNoField,semDiffNeuronLFPFreqNoCueField,...
                pRSDiffNeuronLFPNoCueFreq,'','Neuron mod. freq. - LFP freq.','NeuronModFreq_LFPFreqNoCue',pathAnal,[]);
            
    meanDiffNeuronLFPFreqALField = mean(diffNeuronLFPFreqALField);
    meanDiffNeuronLFPFreqALNoField = mean(diffNeuronLFPFreqALNoField);
    semDiffNeuronLFPFreqALField = std(diffNeuronLFPFreqALField)/sqrt(length(diffNeuronLFPFreqALField));
    semDiffNeuronLFPFreqALNoField = std(diffNeuronLFPFreqALNoField)/sqrt(length(diffNeuronLFPFreqALNoField));
    pRSDiffNeuronLFPALFreq = ranksum(diffNeuronLFPFreqALField,diffNeuronLFPFreqALNoField);
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels') 
    figTitle = '';
    set(figure(figNew),'OuterPosition',...
        [pos(1) pos(2) 280 280],'Name',figTitle)
   
    plotStat(meanDiffNeuronLFPFreqALNoField,meanDiffNeuronLFPFreqALField,...
        semDiffNeuronLFPFreqALNoField,semDiffNeuronLFPFreqALField,...
                pRSDiffNeuronLFPALFreq,'','Neuron mod. freq. - LFP freq.','NeuronModFreq_LFPFreqAL',pathAnal,[]);
end

function plotStat(meanX,meanY,semX,semY,...
                rankXY,xlab,ylab,tit,path,ylimit)    
    meanArr = [meanX meanY];
    errBar = [semX semY];

    h = bar(1:2,meanArr);
    set(h,'faceColor',[0.7,0.7,0.7],'lineWidth',0.5);
    hold on;
    h = errorbar(1:2,meanArr,errBar,'.');
    set(h,'Color',[0 0 0],'lineWidth',1.5);
    
    text(1,1.1*max(meanArr+errBar),num2str(rankXY,'p = %f'));
    xlabel(xlab);
    ylabel(ylab);
    if(~isempty(tit))
        title(tit);
    end
    if(~isempty(ylimit))
        set(gca,'YLim',ylimit);
    end
    print('-painters', '-dpdf', path, '-r600')
    savefig([path 'tit.fig']);
end
