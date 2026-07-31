function plotSpikeRasterAndRunOnT_aligned(path, fileName, onlyRun, mazeSess, trialNum, neuronNo, nGoodTr, cue2RunT, rew2RunT, tit)
% plot spikes rasters across run segments
% E.G.: plotSpikeRaster_aligned('./','A002-20181005-01_DataStructure_mazeSection1_TrialType',1)

    if(nargin == 3)
        trialNum = [];
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
    load(fullPath,'cluList','lap');
    if(isempty(neuronNo))
        neuronNo = 1:length(cluList.all);
    end
    if(isempty(trialNum))
        trialNum = 1:length(lap.trackLen);
    end
    
    fullPath = [path fileName '_alignedSpikesPerNPerT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _alignedSpikesPerNPerT file does not exist');
        return;
    end
    load(fullPath);
    
    GlobalConst;
    
    segLen = 300;
    trialLenT = 10; %sec
    count = 0;
    indTrialCut = find(diff(trialNum)<0);
    if(~isempty(indTrialCut))
        indTrialCut = indTrialCut + 1;
    end
    for i = neuronNo
        [figNew,pos] = CreateFig();
        set(0,'Units','pixels') 
        figTitle = 'Spikes vs Time';
%             figTitle = 'Spikes vs Dist';
        set(figure(figNew),'OuterPosition',...
            [pos(1) pos(2) 600   280],'Name',figTitle)
        
        subplot(1,3,1)
        hold on;
        for j = 1:length(trialNum)      
            % spikes over distance
            plot(trialsRunSpikes.Time{i,trialNum(j)}/sampleFq,...
                 j*ones(1,length(trialsRunSpikes.Time{i,trialNum(j)})),'k.',...
                 'MarkerSize',4.5);    
        end
        if(nGoodTr < length(trialNum))
            plot([0 trialLenT],[nGoodTr+0.5 nGoodTr+0.5],'r');
        end
        set(gca, 'XLim', [0 trialLenT],'YLim',[0 length(trialNum)]);
        ylabel('Trial no.');
        xlabel('Time (s)');
        figTitle = ['Neu ' num2str(i) '(' num2str(cluList.shank(i))...
                    ' ' num2str(cluList.localClu(i)) ')'];
        title(figTitle);
        
        subplot(1,3,2)
        hold on;
        for j = 1:length(trialNum)      
            % spikes over distance
            plot(trialsRewSpikes.Time{i,trialNum(j)}/sampleFq,...
                 j*ones(1,length(trialsRewSpikes.Time{i,trialNum(j)})),'k.',...
                 'MarkerSize',4.5);    
            h = plot([rew2RunT(j),rew2RunT(j)],[j-0.5 j+0.5],'-');
%             set(h,'Color',[0.3 0.9 0.3],'LineWidth',1);
            set(h,'Color',[224 161 224]/255,'LineWidth',1);
        end
        if(nGoodTr < length(trialNum))
            plot([min(rew2RunT) trialLenT],[nGoodTr+0.5 nGoodTr+0.5],'r');
        end
        set(gca, 'XLim', [0 trialLenT],'YLim',[0 length(trialNum)]);
        xlabel('Time (s)');
        
        subplot(1,3,3)
        hold on;
        for j = 1:length(trialNum)      
            % spikes over distance
            plot(trialsCueSpikes.Time{i,trialNum(j)}/sampleFq,...
                 j*ones(1,length(trialsCueSpikes.Time{i,trialNum(j)})),'k.',...
                 'MarkerSize',4.5);    
            h = plot([cue2RunT(j),cue2RunT(j)],[j-0.5 j+0.5],'-');
%             set(h,'Color',[0.3 0.9 0.3],'LineWidth',1);
            set(h,'Color',[224 161 224]/255,'LineWidth',1);
        end
        if(nGoodTr < length(trialNum))
            plot([min(cue2RunT) trialLenT],[nGoodTr+0.5 nGoodTr+0.5],'r');
        end
        set(gca, 'XLim', [0 trialLenT],'YLim',[0 length(trialNum)]);
        xlabel('Time (s)');
        
        pause;
        
        indstr = strfind(fileName,'_');
        fileName1 = [fileName(1:indstr(1)-1) '-Neu' num2str(i) '-Aligned' tit];
        print('-painters','-dpdf',['Z:\Yingxue\DataAnalysisRaphi\' fileName1],'-r600');
    end
