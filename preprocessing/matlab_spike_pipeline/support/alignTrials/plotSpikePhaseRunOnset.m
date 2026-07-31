function plotSpikePhaseRunOnset(path, fileName, onlyRun, mazeSess)
% plot spike rasters and separate trials based on the animal behavior

    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    fullPath = [path fileName '_Info.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to run file does not exist');
        return;
    end
    load(fullPath,'beh');
        
    GlobalConst;
    ind = beh.mazeSess == mazeSess;
    
    if(isfield(beh,'indStimLap'))
        trialNo = find(beh.indStimLap(ind) == 0);
    else
        trialNo = 1:length(behPar.indTrBadBeh);
    end
    
    indBadBeh = behPar.indTrBadBeh;
    
    goodTrials = intersect(find(indBadBeh == 0),trialNo);
    badTrials = intersect(find(indBadBeh == 1),trialNo);
    disp(['Number of good trials = ' num2str(length(goodTrials))]);
    
    plotSpikePhase_aligned(path, fileName, onlyRun, mazeSess, goodTrials, badTrials);
    
%     pause;
%     close all;
end

function plotSpikePhase_aligned(path, fileName, onlyRun, mazeSess, trialNoG, trialNoB, neuronNo)
% plot spikes rasters across run segments
% E.G.: plotSpikeRaster_aligned('./','A002-20181005-01_DataStructure_mazeSection1_TrialType',1)

    if(nargin == 4)
        trialNoG = [];
        trialNoB = [];
        neuronNo = [];
    elseif(nargin == 5)
        trialNoB = [];
        neuronNo = [];
    elseif(nargin == 6)
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
    
    fullPath = [path fileName '_PeakFR_msess' num2str(mazeSess) '_RunOnset' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The peak firing rate aligned to run file does not exist');
        return;
    end
    load(fullPath,'pFRNonStimGoodStruct','pFRNonStimBadStruct');
    
    fullPath = [path fileName '_convSpikesAligned_msess' num2str(mazeSess) '_BefRun' num2str(onlyRun) '.mat'];
    if(exist(fullPath,'file') == 0)
        disp(['The firing profile file does not exist. Try to run the',...
                    'function again with fileState = 0.']);
    end
    load(fullPath,'timeStepRun','paramC');
            
    if(isempty(trialNoG))
        trialNoG = 1:size(trialsRunSpikes.Time,2);
        trialNoB = [];
    end
    
    GlobalConst;
    
    segLen = 300;
    trialLenT = 20; %sec
    count = 0;

    for i = neuronNo
        disp(['Neuron ' num2str(i)]);
        count = count + 1;
        if(mod(count-1,12) == 0)
            [figNew,pos] = CreateFig();
            set(0,'Units','pixels') 
            figTitle = 'Spikes vs Time';
%             figTitle = 'Spikes vs Dist';
            set(figure(figNew),'OuterPosition',...
                [pos(1) pos(2)-300 pos(3)*1.5 pos(4)*1.7],'Name',figTitle)
        end
        
        subplot(6,2,mod(count-1,12)+1)
        yyaxis left
        hold on;
        for j = trialNoG   
            disp(['Trial ' num2str(j)])
            plot([trialsRunSpikes.TimeBef{i,j}'/sampleFq ...
                 trialsRunSpikes.Time{i,j}'/sampleFq],...
                 [trialsRunSpikes.thPhaseInterpSpikeBef{i,j}' ...
                 trialsRunSpikes.thPhaseInterpSpike{i,j}'],'k.',...
                 'MarkerSize',3);     
        end
        if(mod(count-1,2) == 0)
            ylabel('Theta phase');
        end
        yyaxis right
        h = plot(timeStepRun/sampleFq,pFRNonStimGoodStruct.avgFRProfile(i,:),'r-');
        set(h,'LineWidth',0.5);
        if(mod(count-1,2) == 0)
            ylabel('FR (Hz)');
        end
        set(gca, 'XLim', [-3 trialLenT]);
        
        if(mod(count-1,12) > 9)
            xlabel('Time (s)');
        end
        figTitle = ['Neu ' num2str(i) '(' num2str(cluList.shank(i))...
                    ' ' num2str(cluList.localClu(i)) ')'];
        title(figTitle);
        
        count = count + 1;
        subplot(6,2,mod(count-1,12)+1)
        yyaxis left
        hold on;
        for j = trialNoB   
            disp(['Trial ' num2str(j)])
            plot([trialsRunSpikes.TimeBef{i,j}'/sampleFq ...
                 trialsRunSpikes.Time{i,j}'/sampleFq],...
                 [trialsRunSpikes.thPhaseInterpSpikeBef{i,j}' ...
                 trialsRunSpikes.thPhaseInterpSpike{i,j}'],'k.',...
                 'MarkerSize',3);     
        end
        if(mod(count-1,2) == 0)
            ylabel('Theta phase');
        end
        yyaxis right
        if(isfield(pFRNonStimBadStruct,'avgFRProfile'))
            h = plot(timeStepRun/sampleFq,pFRNonStimBadStruct.avgFRProfile(i,:),'r-');
            set(h,'LineWidth',0.5);
            if(mod(count-1,2) == 0)
                ylabel('FR (Hz)');
            end
        end        
        set(gca, 'XLim', [-3 trialLenT]);
        if(mod(count-1,2) == 0)
            ylabel('Theta phase');
        end
        if(mod(count-1,12) > 9)
            xlabel('Time (s)');
        end
        figTitle = ['Neu ' num2str(i) '(' num2str(cluList.shank(i))...
                    ' ' num2str(cluList.localClu(i)) ')'];
        title(figTitle);
        
    end
end
