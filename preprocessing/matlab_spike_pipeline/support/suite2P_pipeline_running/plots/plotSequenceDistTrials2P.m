function plotSequenceDistTrials2P(path, fileName, onlyRun, mazeSess, intervalT, spaceBin, neuSel, neuArr, trialNo, seqLabel)
% plot the sequence from one recording, align the sequence to different
% behavioral features
% examples: 
% 1. plot all the neurons with fields: plotSequenceDistTrials2P('./','A576-20210914-02_DataStructure_mazeSection1_TrialType1',1,1,10,20,0,[],1:5,[])
% 2. plot a selected array of neurons: plotSequenceDistTrials2P('./','A576-20210914-02_DataStructure_mazeSection1_TrialType1',1,1,10,20,2,[1 3 5],1:5,'PC1')

    GlobalConst2P;
     
    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
   
    fileNameCorr = [fileName '_meanSpikesCorrTAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalT) '.mat'];
    fullPath = [path fileNameCorr];
    if(exist(fullPath) == 0)
        disp('The _meanSpikesCorrTAligned_Run file does not exist');
        return;
    end    
    load(fullPath,'meanCorrTRun');
    
    if(neuSel == 1)
        thrCorrT = 0.2;
        indSelCorrT = selPyrNeurons2P(path,fileName,onlyRun,minFR);
        indSelCorrT = meanCorrTRun.meanGood > thrCorrT & indSelCorrT;
        neuronNo = find(indSelCorrT == 1);
        filenameEnd = '_CorrT';
    elseif(neuSel == 2)
        if(~isempty(neuArr))
            neuronNo = neuArr;
        else
            neuronNo = [1:5];
        end
        %neuronNo = [2 6 8 18 25 30 40 41 47 51 55 57 67 68 69 71 72 73 77 78 84];
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
    
    fullPath = [path fileName '_convSpikesDist' num2str(spaceBin) ...
                    'mm_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _convSpikesDist file does not exist');
        return;
    end
    load(fullPath,'filteredSpikeArrayNormT','paramC');
    spaceSteps =  paramC.spaceSteps;
    
    if(~isempty(trialNo))
        
        avgFilteredSpikeArrayRunGood...
                    = avgFilteredSpikeArray(filteredSpikeArrayNormT,trialNo,neuronNo);
       
        indPeak = zeros(1,length(neuronNo));
        for i = 1:length(neuronNo)
            [~,indPeakTmp] = max(avgFilteredSpikeArrayRunGood(i,:));
            indPeak(i) = indPeakTmp(1);
        end
        [~,indNeuronOrder] = sort(indPeak);
        avgFilteredSpikeArrayRunGood = avgFilteredSpikeArrayRunGood(indNeuronOrder,:);

        plotSequenceDist1(avgFilteredSpikeArrayRunGood,...
            spaceSteps{1}/10,'Dist (cm)','Neuron No.','Firing rate good trials');
        indstr = strfind(fileName,'_');
        fileName1 = [fileName(1:indstr(1)-1) '-SeqDist-SelTr' filenameEnd];
        if(~isempty(seqLabel))
            fileName1 = [fileName1 '-' seqLabel];
        end
        savefig([path fileName1 '.fig'])
        print('-painters','-dpdf',[path fileName1],'-r600');
    %     savefig(['Z:\Yingxue\DataAnalysisXiaoliang\2PData\' fileName1 '.fig'])
    %     print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisXiaoliang\2PData\' fileName1],'-r600');
    else
        disp('Number of trials = 0')
    end
end

function [avgFilteredSpikeArr]...
                = avgFilteredSpikeArray(filteredSpikeArr,indTr,neuronNo)
    if(isempty(indTr))
        avgFilteredSpikeArr = [];
        return;
    end
    lenTr = size(filteredSpikeArr{indTr(1)},2);
    nTrials = length(indTr);
    nNeurons = length(neuronNo);
    avgFilteredSpikeArr = zeros(nNeurons,lenTr);
    for i = 1:nNeurons
        avgFilteredSpikeTmp = zeros(1,lenTr);
        for j = 1:nTrials
            avgFilteredSpikeTmp = avgFilteredSpikeTmp + filteredSpikeArr{indTr(j)}(neuronNo(i),:);
        end
        avgFilteredSpikeArr(i,:) = avgFilteredSpikeTmp/nTrials/max(avgFilteredSpikeTmp);
    end      
end

function plotSequenceDist1(arr1,spaceSteps,xl,yl,title)
    [figNew,pos] = CreateFig2P();
    set(0,'Units','pixels')
    set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 200 280],'Name',title)
    nNeurons = size(arr1,1);
    imagesc(spaceSteps,1:nNeurons,arr1);
    colormap jet
%     set(gca,'XLim',[0 20])
    xlabel(xl);
    ylabel(yl);
   
end