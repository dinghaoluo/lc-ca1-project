function plotDeltaOscLFPFreqFieldNoField()

    onlyRun = 1;
    RecordingList;
    pathAnal = 'Z:\Yingxue\DataAnalysisRaphi\Pyramidal\';
    
    if(exist([pathAnal 'autoCorrPyrAllRec.mat']))
        load([pathAnal 'autoCorrPyrAllRec.mat'],...
            'autoCorrPyrNoCue','autoCorrPyrAL','autoCorrPyrPL',...
            'modPyrNoCue','modPyrAL','modPyrPL');
    end
    mod.isNeuWithField = [modPyrNoCue.isNeuWithField modPyrAL.isNeuWithField modPyrPL.isNeuWithField];
    mod.diffNeuronLFPFreq = [modPyrNoCue.thetaModFreq3-modPyrNoCue.thetaFreqHMean...
        modPyrAL.thetaModFreq3-modPyrAL.thetaFreqHMean ... 
        modPyrPL.thetaModFreq3-modPyrPL.thetaFreqHMean];
    ind = mod.isNeuWithField == 1;
    ind2 = mod.isNeuWithField == 0;
    diffNeuronLFPFreqField = mod.diffNeuronLFPFreq(ind);
    diffNeuronLFPFreqNoField = mod.diffNeuronLFPFreq(ind2);
    meanDiffNeuronLFPFreqField = mean(diffNeuronLFPFreqField);
    meanDiffNeuronLFPFreqNoField = mean(diffNeuronLFPFreqNoField);
    semDiffNeuronLFPFreqField = std(diffNeuronLFPFreqField)/sqrt(length(diffNeuronLFPFreqField));
    semDiffNeuronLFPFreqNoField = std(diffNeuronLFPFreqNoField)/sqrt(length(diffNeuronLFPFreqNoField));
    pRSDiffNeuronLFPFreq = ranksum(diffNeuronLFPFreqField,diffNeuronLFPFreqNoField);
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels') 
    figTitle = '';
    set(figure(figNew),'OuterPosition',...
        [pos(1) pos(2) 280 280],'Name',figTitle)
   
    plotStat(meanDiffNeuronLFPFreqField,meanDiffNeuronLFPFreqNoField,...
        semDiffNeuronLFPFreqField,semDiffNeuronLFPFreqNoField,...
                pRSDiffNeuronLFPFreq,'','Neuron mod. freq. - LFP freq.','NeuronModFreq_LFPFreq',pathAnal,[]);  
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