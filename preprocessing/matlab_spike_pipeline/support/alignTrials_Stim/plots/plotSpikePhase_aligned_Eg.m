function plotSpikePhase_aligned_Eg(path, fileName, onlyRun, mazeSess, neuronNo)
% plot spikes phase vs time, and the instantaneous firing rate change
% E.G.: plotSpikePhase_aligned_Eg('./','A016-20190531-01_DataStructure_mazeSection1_TrialType1',0,5,62)
% plotSpikePhase_aligned_Eg('./','A022-20191107-01_DataStructure_mazeSection1_TrialType1',0,1,91)
% plotSpikePhase_aligned_Eg('./','A016-20190603-01_DataStructure_mazeSection1_TrialType1',0,3,49)
% plotSpikePhase_aligned_Eg('./','A016-20190531-01_DataStructure_mazeSection1_TrialType1',0,5,19)

    if(nargin == 4)
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
    trialNoG = intersect(find(indBadBeh == 0),trialNo);
    trialNoB = [];   
    timeBin = 0.0025; %sec
        
    GlobalConst;
    
    segLen = 300;
    trialLenT = 3; %sec
    count = 0;

    for i = neuronNo
        disp(['Neuron ' num2str(i)]);
        [figNew,pos] = CreateFig();
        set(0,'Units','pixels') 
        set(figure(figNew),'OuterPosition',...
            [500 500 350 350])
                
        yyaxis left
        hold on;
        time = [];
        thPhase = [];
        for j = trialNoG   
            disp(['Trial ' num2str(j)])
            time = [time trialsRunSpikes.TimeBef{i,j}'/sampleFq ...
                 trialsRunSpikes.Time{i,j}'/sampleFq];
            tmp = [trialsRunSpikes.thPhaseInterpSpikeBef{i,j}' ...
                 trialsRunSpikes.thPhaseInterpSpike{i,j}'];
            thPhase = [thPhase tmp];
%             plot([trialsRunSpikes.TimeBef{i,j}'/sampleFq ...
%                  trialsRunSpikes.Time{i,j}'/sampleFq],...
%                  [trialsRunSpikes.thPhaseInterpSpikeBef{i,j}' ...
%                  trialsRunSpikes.thPhaseInterpSpike{i,j}'],'k.',...
%                  'MarkerSize',3);     
        end
        [timeDouble,thPhaseDouble] = ...
                getMultiCycles(time',thPhase',3);
        h = plot(timeDouble,thPhaseDouble,'.');
        set(h,'MarkerSize',3,'Color',[0 0 0]);
        set(gca,'FontSize',14.0,'YLim', [0 720]); 
            
        ylabel('Theta phase (deg.)');
        
        yyaxis right
        h = plot(timeStepRun/sampleFq,pFRNonStimGoodStruct.avgFRProfile(i,:),'r-');
        set(h,'LineWidth',0.5);
        ylabel('FR (Hz)');        
        set(gca, 'XLim', [-1 trialLenT]);       
        xlabel('Time (s)');
        figTitle = ['Neu ' num2str(i) '(' num2str(cluList.shank(i))...
                    ' ' num2str(cluList.localClu(i)) ')'];
        title(figTitle);
        
        ind = findstr(fileName,'_');
        print ('-painters', '-dpdf', ['spikePhaseAligned_' fileName(1:ind(1)-1) 'Neu' num2str(i)], '-r600')
        
    end
end