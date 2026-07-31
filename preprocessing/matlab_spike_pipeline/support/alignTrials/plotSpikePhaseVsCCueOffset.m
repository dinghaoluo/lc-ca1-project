function plotSpikePhaseVsCCueOffset(path, fileName, onlyRun, mazeSess, methodTheta)
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
    
    plotSpikePhase_aligned(path, fileName, onlyRun, mazeSess, methodTheta, goodTrials, badTrials);
end

function plotSpikePhase_aligned(path, fileName, onlyRun, mazeSess, methodTheta, trialNoG, trialNoB, neuronNo)
% plot spikes rasters across run segments
% E.G.: plotSpikeRaster_aligned('./','A002-20181005-01_DataStructure_mazeSection1_TrialType',1)

    if(nargin == 5)
        trialNoG = [];
        trialNoB = [];
        neuronNo = [];
    elseif(nargin == 6)
        trialNoB = [];
        neuronNo = [];
    elseif(nargin == 7)
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
       
    if(methodTheta == 0)
        fileNameThetaPhase = [fileName '_ThetaPhaseH_msess' num2str(mazeSess) '_CueOffRun' num2str(onlyRun) '.mat'];
    else
        fileNameThetaPhase = [fileName '_ThetaPhaseL_msess' num2str(mazeSess) '_CueOffRun' num2str(onlyRun) '.mat'];
    end
    fullPath = [path fileNameThetaPhase];
    if(exist(fullPath,'file') == 0)
        disp(['The _ThetaPhase_CueOff does not exist.']);
    end
    load(fullPath,'spikeThetaPhaseCueOffNoStimGood','spikeThetaPhaseCueOffNoStimBad');
        
    if(isempty(trialNoG))
        trialNoG = 1:size(trialsCueOffSpikes.Time,2);
        trialNoB = [];
    end
    
    GlobalConst;
  
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
                [pos(1) pos(2) pos(3)*2 pos(4)*2.2],'Name',figTitle)
        end
        
        subplot(6,2,mod(count-1,12)+1)
        hold on;
        plot(spikeThetaPhaseCueOffNoStimGood.spCycleRunOnPerNeuron{i}, ...
                 spikeThetaPhaseCueOffNoStimGood.spPhaseVsCPerNeuron{i},'k.',...
                 'MarkerSize',3);     
        if(isempty(spikeThetaPhaseCueOffNoStimGood.spCycleRunOnPerNeuron{i})) 
            continue;
        end
        minCycle = min(spikeThetaPhaseCueOffNoStimGood.spCycleRunOnPerNeuron{i});
        maxCycle = max(spikeThetaPhaseCueOffNoStimGood.spCycleRunOnPerNeuron{i});
        h = plot(minCycle:maxCycle,spikeThetaPhaseCueOffNoStimGood.meanPhaseVsCRunOn{i},'ro');
        set(h,'LineWidth',0.5);
        h = errorbar(minCycle:maxCycle,spikeThetaPhaseCueOffNoStimGood.meanPhaseVsCRunOn{i},...
            spikeThetaPhaseCueOffNoStimGood.stdPhaseVsCRunOn{i},'r');
        set(h,'LineWidth',0.3);
        set(gca, 'XLim', [minCycle maxCycle]);
        if(mod(count-1,2) == 0)
            ylabel('Theta phase');
        end
        if(mod(count-1,12) > 9)
            xlabel('Cycle no.');
        end
        figTitle = ['Neu ' num2str(i) '(' num2str(cluList.shank(i))...
                    ' ' num2str(cluList.localClu(i)) ')'];
        title(figTitle);
        
        count = count + 1;
        subplot(6,2,mod(count-1,12)+1)
        hold on;
        plot(spikeThetaPhaseCueOffNoStimBad.spCycleRunOnPerNeuron{i}, ...
                 spikeThetaPhaseCueOffNoStimBad.spPhaseVsCPerNeuron{i},'k.',...
                 'MarkerSize',3);     
             
        minCycle1 = min(spikeThetaPhaseCueOffNoStimBad.spCycleRunOnPerNeuron{i});
        maxCycle1 = max(spikeThetaPhaseCueOffNoStimBad.spCycleRunOnPerNeuron{i});
        h = plot(minCycle1:maxCycle1,spikeThetaPhaseCueOffNoStimBad.meanPhaseVsCRunOn{i},'ro');
        set(h,'LineWidth',0.5);
        h = errorbar(minCycle1:maxCycle1,spikeThetaPhaseCueOffNoStimBad.meanPhaseVsCRunOn{i},...
            spikeThetaPhaseCueOffNoStimBad.stdPhaseVsCRunOn{i},'r');
        set(h,'LineWidth',0.3);
        set(gca, 'XLim', [minCycle maxCycle]);
        if(mod(count-1,2) == 0)
            ylabel('Theta phase');
        end
        if(mod(count-1,12) > 9)
            xlabel('Cycle no.');
        end
        figTitle = ['Neu ' num2str(i) '(' num2str(cluList.shank(i))...
                    ' ' num2str(cluList.localClu(i)) ')'];
        title(figTitle);
        
    end
end
