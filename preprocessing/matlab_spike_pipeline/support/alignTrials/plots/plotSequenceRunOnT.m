function plotSequenceRunOnT(path, fileName, onlyRun, mazeSess, intervalT,neuSel)
% plot the sequence from one recording, align the sequence to different
% behavioral features

    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    indBadBeh = behPar.indTrBadBeh;
    
    goodTrials = find(indBadBeh == 0);
    badTrials = find(indBadBeh == 1);
    disp(['Number of good trials = ' num2str(length(goodTrials))]);
    
    fileNameCorr = [fileName '_meanSpikesCorrTAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    fullPath = [path fileNameCorr];
    if(exist(fullPath) == 0)
        disp('The _meanSpikesCorrTAligned_Run file does not exist');
        return;
    end    
    load(fullPath,'meanCorrTRun');
    
    GlobalConst;
    if(neuSel == 1)
        thrCorrT = 0.25;
        indSelCorrT = selPyrNeurons(path,fileName,onlyRun,minFR);
        indSelCorrT = meanCorrTRun.meanGood > thrCorrT & indSelCorrT;
        neuronNo = find(indSelCorrT == 1);
        filenameEnd = '_CorrT';
    elseif(neuSel == 2)
        neuronNo = [2 6 8 18 25 30 40 41 47 51 55 57 67 68 69 71 72 73 77 78 84];
        filenameEnd = '_Sel';
    else
        fileNameFW = [fileName '_FieldSpCorrAligned_Run' num2str(mazeSess) ...
                            '_Run' num2str(onlyRun) '.mat'];
        fullPath = [path fileNameFW];
        if(exist(fullPath) == 0)
            disp('The _FieldSpCorrAligned_Run file does not exist');
            return;
        end
        load(fullPath,'fieldSpCorrSessNonStimGood');
        neuronNo = fieldSpCorrSessNonStimGood.indNeuron;
        filenameEnd = '_Field';
    end
    
    fullPath = [path fileName '_convSpikesAligned_msess' num2str(mazeSess) ...
        '_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _convSpikesAligned_Run file does not exist');
        return;
    end
    load(fullPath,'filteredSpikeArrayRun','filteredSpikeArrayRew','filteredSpikeArrayCue',...
            'timeStep');
    timeStep1 = timeStep;
    GlobalConst;
    
    avgFilteredSpikeArrayRunGood...
                = avgFilteredSpikeArray(filteredSpikeArrayRun,goodTrials,neuronNo);
    avgFilteredSpikeArrayRewGood...
                = avgFilteredSpikeArray(filteredSpikeArrayRew,goodTrials,neuronNo);
    avgFilteredSpikeArrayCueGood...
                = avgFilteredSpikeArray(filteredSpikeArrayCue,goodTrials,neuronNo);
    
    avgFilteredSpikeArrayRunBad...
                = avgFilteredSpikeArray(filteredSpikeArrayRun,badTrials,neuronNo);
    avgFilteredSpikeArrayRewBad...
                = avgFilteredSpikeArray(filteredSpikeArrayRew,badTrials,neuronNo);
    avgFilteredSpikeArrayCueBad...
                = avgFilteredSpikeArray(filteredSpikeArrayCue,badTrials,neuronNo);
            
    indPeak = zeros(1,length(neuronNo));
    for i = 1:length(neuronNo)
        [~,indPeakTmp] = max(avgFilteredSpikeArrayRunGood(i,:));
        indPeak(i) = indPeakTmp(1);
    end
    [~,indNeuronOrder] = sort(indPeak);
    avgFilteredSpikeArrayRunGood = avgFilteredSpikeArrayRunGood(indNeuronOrder,:);
    avgFilteredSpikeArrayRewGood = avgFilteredSpikeArrayRewGood(indNeuronOrder,:);
    avgFilteredSpikeArrayCueGood = avgFilteredSpikeArrayCueGood(indNeuronOrder,:);
    
    avgFilteredSpikeArrayRunBad = avgFilteredSpikeArrayRunBad(indNeuronOrder,:);
    avgFilteredSpikeArrayRewBad = avgFilteredSpikeArrayRewBad(indNeuronOrder,:);
    avgFilteredSpikeArrayCueBad = avgFilteredSpikeArrayCueBad(indNeuronOrder,:);
            
    plotSequenceRunOnT_aligned(avgFilteredSpikeArrayRunGood,...
        avgFilteredSpikeArrayRewGood,avgFilteredSpikeArrayCueGood,...
        timeStep1/sampleFq,'Time (s)','Neuron No.','Firing rate good trials');
    indstr = strfind(fileName,'_');
    fileName1 = [fileName(1:indstr(1)-1) '-Seq-GoodTr' filenameEnd];
    savefig([fileName1 '.fig'])
    print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');
    
    if(~isempty(badTrials))
        plotSequenceRunOnT_aligned(avgFilteredSpikeArrayRunBad,...
            avgFilteredSpikeArrayRewBad,avgFilteredSpikeArrayCueBad,...
            timeStep1/sampleFq,'Time (s)','Neuron No.','Firing rate bad trials');
        indstr = strfind(fileName,'_');
        fileName1 = [fileName(1:indstr(1)-1) '-Seq-BadTr' filenameEnd];
        savefig([fileName1 '.fig'])
        print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');
    end
end

function [avgFilteredSpikeArr]...
                = avgFilteredSpikeArray(filteredSpikeArr,indTr,neuronNo)
    lenTr = size(filteredSpikeArr{1},2);
    nTrials = length(indTr);
    nNeurons = length(neuronNo);
    avgFilteredSpikeArr = zeros(nNeurons,lenTr);
    for i = 1:nNeurons
        avgFilteredSpikeTmp = zeros(1,lenTr);
        for j = 1:nTrials
            avgFilteredSpikeTmp = avgFilteredSpikeTmp + filteredSpikeArr{neuronNo(i)}(indTr(j),:);
        end
        avgFilteredSpikeArr(i,:) = avgFilteredSpikeTmp/nTrials/max(avgFilteredSpikeTmp);
    end      
end

function plotSequenceRunOnT_aligned(arr1,arr2,arr3,timeStep,xl,yl,title)
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 600   280],'Name',title)
    nNeurons = size(arr1,1);
    subplot(1,3,1)
    imagesc(timeStep,1:nNeurons,arr1);
    colormap jet
    set(gca,'XLim',[0 10])
    xlabel(xl);
    ylabel(yl);
    
    subplot(1,3,2)
    imagesc(timeStep,1:nNeurons,arr2);
    colormap jet
    set(gca,'XLim',[0 10])
    xlabel(xl);
    
    subplot(1,3,3)
    imagesc(timeStep,1:nNeurons,arr3);
    colormap jet
    set(gca,'XLim',[0 10])
    xlabel(xl);
end