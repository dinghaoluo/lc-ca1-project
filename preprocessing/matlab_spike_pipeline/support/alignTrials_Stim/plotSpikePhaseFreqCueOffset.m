function plotSpikePhaseFreqCueOffset(path, fileName, onlyRun, mazeSess)
% plot spike phase and theta frequency over time

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
    load(fullPath,'cluList');
    if(isempty(neuronNo))
        neuronNo = 1:length(cluList.all);
    end
       
    fullPath = [path fileName '_alignedSpikesPerNPerT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _alignedSpikesPerNPerT file does not exist');
        return;
    end
    load(fullPath);
    fullPath = [path fileName '_alignCueOff_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to run file does not exist');
        return;
    end
    load(fullPath,'trialsCueOff');
            
    if(isempty(trialNoG))
        trialNoG = 1:size(trialsCueOffSpikes.Time,2);
        trialNoB = [];
    end
    
    GlobalConst;

    trialLenT = 10; %sec
    [figNew,pos] = CreateFig();
    set(0,'Units','pixels') 
    figTitle = 'Spikes vs Time';
%             figTitle = 'Spikes vs Dist';
    set(figure(figNew),'OuterPosition',...
        [pos(1) pos(2)-300 pos(3)*1.5 pos(4)*1.7],'Name',figTitle)
        
    subplot(2,2,1)
    hold on;
    for j = trialNoG   
%         disp(['Trial ' num2str(j)])
        if(trialsCueOff.startLfpInd{j} == -1)
            continue;
        end
        if(isempty(trialsCueOff.ThetaFreq{j}))
            continue;
        end
        plot([-nSampBef+1:trialsCueOff.numSamples(j)]/sampleFq,...
             [trialsCueOff.ThetaFreqBef{j}' trialsCueOff.ThetaFreq{j}'],'k.','MarkerSize',3);     
    end
    set(gca, 'XLim', [-nSampBef/sampleFq trialLenT],'YLim',[-10 40]);
    ylabel('Theta freq (Hz)');
       
    subplot(2,2,2)
    hold on;
    for j = trialNoB   
%         disp(['Trial ' num2str(j)])
        if(trialsCueOff.startLfpInd{j} == -1)
            continue;
        end
        if(isempty(trialsCueOff.ThetaFreq{j}))
            continue;
        end
        plot([-nSampBef+1:trialsCueOff.numSamples(j)]/sampleFq,...
             [trialsCueOff.ThetaFreqBef{j}' trialsCueOff.ThetaFreq{j}'],'k.','MarkerSize',3);    
    end
    set(gca, 'XLim', [-nSampBef/sampleFq trialLenT],'YLim',[-10 40]);
    
    subplot(2,2,3)
    hold on;
    for j = trialNoG   
%         disp(['Trial ' num2str(j)])
        if(trialsCueOff.startLfpInd{j} == -1)
            continue;
        end
        if(isempty(trialsCueOff.ThetaFreq{j}))
            continue;
        end
        plot([-nSampBef+1:trialsCueOff.numSamples(j)]/sampleFq,...
             [trialsCueOff.ThetaAmpBef{j}' trialsCueOff.ThetaAmp{j}'],'k.','MarkerSize',3);     
    end
    set(gca, 'XLim', [-nSampBef/sampleFq trialLenT]);
    ylabel('Theta amp');
    xlabel('Time (s)');
    
    subplot(2,2,4)
    hold on;
    for j = trialNoB   
%         disp(['Trial ' num2str(j)])
        if(trialsCueOff.startLfpInd{j} == -1)
            continue;
        end
        if(isempty(trialsCueOff.ThetaFreq{j}))
            continue;
        end
        plot([-nSampBef+1:trialsCueOff.numSamples(j)]/sampleFq,...
             [trialsCueOff.ThetaAmpBef{j}' trialsCueOff.ThetaAmp{j}'],'k.','MarkerSize',3);     
    end
    set(gca, 'XLim', [-nSampBef/sampleFq trialLenT]);
    ylabel('Theta amp');
    xlabel('Time (s)');
end
