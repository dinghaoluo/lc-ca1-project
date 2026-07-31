function plotSpikeRasterEg(path, fileName, onlyRun, trialNo, neuronNo)
% plot spikes rasters across run segments
% E.G.: plotSpikeRasterEg('./','A007-20190116-01_DataStructure_mazeSection1_TrialType1',1,1,[1:49],36)

    if(nargin == 3)
        trialNo = [];
        neuronNo = [];
    elseif(nargin == 4)
        neuronNo = [];
    end

    %%%%%%%%% load recording file
    fullPath = [path fileName '.mat']; 
    if(exist(fullPath) == 0)
        disp('The recording file does not exist');
        return;
    end
    load(fullPath,'cluList','lap','trials');
    if(isempty(neuronNo))
        neuronNo = 1:length(cluList.all);
    end
    if(isempty(trialNo))
        trialNo = 1:length(trials);
    end
    
    GlobalConst;
    
    segLen = 300;
    trialLen = 180; %cm
    count = 0;
    indTrialCut = find(diff(trialNo)<0);
    if(~isempty(indTrialCut))
        indTrialCut = indTrialCut + 1;
    end
    for i = neuronNo
        disp(['Neuron ' num2str(i)]);
        [figNew,pos] = CreateFig();
        set(0,'Units','pixels') 
        figTitle = 'Spikes vs Time';
%             figTitle = 'Spikes vs Dist';
        set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 280 280],'Name',figTitle)
        
        hold on;
        for j = 1:length(trialNo)   
            indTr = trialNo(j);
            if(onlyRun == 1)
                ind = trials{indTr}.speed(trials{indTr}.spikes{i}) > minSpeed;
            else
                ind = trials{indTr}.speed(trials{indTr}.spikes{i}) >= 0;
            end
            plot(trials{indTr}.xMM(trials{indTr}.spikes{i}(ind))'/10,...
                 j*ones(1,sum(ind)),'k.',...
                 'MarkerSize',4);    
             
%             plot([trialsRunSpikes.TimeBef{i,j}'/sampleFq ...
%                  trialsRunSpikes.Time{i,j}'/sampleFq],...
%                  [trialsRunSpikes.thPhaseInterpSpikeBef{i,j}' ...
%                  trialsRunSpikes.thPhaseInterpSpike{i,j}'],'k.',...
%                  'MarkerSize',3);     
        end
        if(~isempty(indTrialCut))
            plot([0 trialLen],[indTrialCut indTrialCut],'r');
        end
        set(gca, 'XLim', [0 trialLen]);
        ylabel('Trial no.');
        xlabel('Distance (cm)');
        figTitle = ['Neu ' num2str(i) '(' num2str(cluList.shank(i))...
                    ' ' num2str(cluList.localClu(i)) ')'];
        title(figTitle);
        
        ind = strfind(fileName,'_');
        fileName1 = ['Z:\Yingxue\DataAnalysisRaphi\' ...
            fileName(1:ind(1)) 'Neu' num2str(i) '_SpikeRaster'];
        saveas(gcf,fileName1);
        print('-painters', '-dpdf', fileName1, '-r600')
    end