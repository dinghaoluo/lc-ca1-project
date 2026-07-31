function plotSpikeRaster_alignedRunVsCueOff(path, fileName, onlyRun, mazeSess, trialNo, neuronNo)
% plot spikes rasters across run segments
% E.G.: plotSpikeRaster_aligned('./','A002-20181005-01_DataStructure_mazeSection1_TrialType',1)

    if(nargin == 4)
        trialNo = [];
        neuronNo = [];
    elseif(nargin == 5)
        neuronNo = [];
    end

    %%%%%%%%% load recording file
    fullPath = [path fileName '.mat']; 
    if(exist(fullPath) == 0)
        disp('The recording file does not exist');
        return;
    end
    load(fullPath,'cluList','lap');
    if(isempty(neuronNo))
        neuronNo = 1:length(cluList.all);
    end
       
    fullPath = [path fileName '_alignedSpikesPerNPerT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _alignedSpikesPerNPerT file does not exist');
        return;
    end
    load(fullPath);
    if(isempty(trialNo))
        trialNo = 1:size(trialsRunSpikes.Time,2);
    end
    
    GlobalConst;
    
    segLen = 300;
    trialLenT = 20; %sec
    count = 0;
    indTrialCut = find(diff(trialNo)<0);
    if(~isempty(indTrialCut))
        indTrialCut = indTrialCut + 1;
    end
    for i = neuronNo
        disp(['Neuron ' num2str(i)]);
        count = count + 1;
        if(mod(count-1,18) == 0)
            [figNew,pos] = CreateFig();
            set(0,'Units','pixels') 
            figTitle = 'Spikes vs Time';
%             figTitle = 'Spikes vs Dist';
            set(figure(figNew),'OuterPosition',...
                [pos(1) pos(2)-300 pos(3)*1.5 pos(4)*1.7],'Name',figTitle)
        end
        
        subplot(6,3,mod(count-1,18)+1)
        hold on;
        for j = trialNo   
%             disp(['Trial ' num2str(j)])
%             if j == 103
%                 a = 1;
%             end
            % spikes over distance
            plot([trialsRunSpikes.TimeBef{i,j}'/sampleFq ...
                 trialsRunSpikes.Time{i,j}'/sampleFq],...
                 j*ones(1,length([trialsRunSpikes.TimeBef{i,j}' ...
                 trialsRunSpikes.Time{i,j}'])),'k.',...
                 'MarkerSize',3);    
             
%             plot([trialsRunSpikes.TimeBef{i,j}'/sampleFq ...
%                  trialsRunSpikes.Time{i,j}'/sampleFq],...
%                  [trialsRunSpikes.thPhaseInterpSpikeBef{i,j}' ...
%                  trialsRunSpikes.thPhaseInterpSpike{i,j}'],'k.',...
%                  'MarkerSize',3);     
        end
        if(~isempty(indTrialCut))
            plot([0 trialLenT],[indTrialCut indTrialCut],'r');
        end
        set(gca, 'XLim', [-3 trialLenT]);
        if(mod(count-1,3) == 0)
            ylabel('Trial no.');
        end
        if(mod(count-1,18) > 14)
            xlabel('Time (s)');
        end
        figTitle = ['Neu ' num2str(i) '(' num2str(cluList.shank(i))...
                    ' ' num2str(cluList.localClu(i)) ')'];
        title(figTitle);
        
        count = count + 1;
        subplot(6,3,mod(count-1,18)+1)
        hold on;
        for j = trialNo      
            % spikes over distance
            plot(trialsCueSpikes.Time{i,j}'/sampleFq,...
                 j*ones(1,length(trialsCueSpikes.Time{i,j})),'k.',...
                 'MarkerSize',3);    
        end
        if(~isempty(indTrialCut))
            plot([0 trialLenT],[indTrialCut indTrialCut],'r');
        end
        set(gca, 'XLim', [-3 trialLenT]);
        if(mod(count-1,3) == 0)
            ylabel('Trial no.');
        end
        if(mod(count-1,18) > 14)
            xlabel('Time (s)');
        end
        figTitle = ['Neu ' num2str(i) '(' num2str(cluList.shank(i))...
                    ' ' num2str(cluList.localClu(i)) ')'];
        title(figTitle);
        
        count = count + 1;
        subplot(6,3,mod(count-1,18)+1)
        hold on;
        for j = trialNo      
            % spikes over distance
            plot([trialsCueOffSpikes.TimeBef{i,j}'/sampleFq ...
                 trialsCueOffSpikes.Time{i,j}'/sampleFq],...
                 j*ones(1,length([trialsCueOffSpikes.TimeBef{i,j}' ...
                 trialsCueOffSpikes.Time{i,j}'])),'k.',...
                 'MarkerSize',3);    
        end
        if(~isempty(indTrialCut))
            plot([0 trialLenT],[indTrialCut indTrialCut],'r');
        end
        set(gca, 'XLim', [-3 trialLenT]);
        if(mod(count-1,3) == 0)
            ylabel('Trial no.');
        end
        if(mod(count-1,18) > 14)
            xlabel('Time (s)');
        end
        figTitle = ['Neu ' num2str(i) '(' num2str(cluList.shank(i))...
                    ' ' num2str(cluList.localClu(i)) ')'];
        title(figTitle);
    end