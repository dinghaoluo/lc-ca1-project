function PeakFiringRate_Aligned(path,fileName,fileState,onlyRun,mazeSess,figureState)
% Calculate the smoothed mean firing rate curve over trials and the peak
% firing rate for each recorded neuron
%             if fileState == 1, function "ConvSpikeTrain" should be
%             executed first
% path:         the path of the recording file
% fileName:     name of the recording file
% spaceBin:      2SD of the Gaussian filter used to obtain the firing rate
%               profile (in mm), default value is 10 mm
% fileState:    0: calculate using the recorded data (default)
%               1: load the firing rate profile if the file exists, and
%               calculate the peak firing rate from there
% onlyRun:      1: only consider the time period when the animal is running 
% figureState:  0: figure off
%               else: figure on
%
% Example: 
% PeakFiringRateVR('./','A111-20150301-01_DataStructure_mazeSection1_TrialType1',1,1,1,1)

    %%%%%%%%% check arguments
    if nargin<2
        disp('At least two arguments are needed for this function.');
        return;
    elseif nargin == 2
        fileState = 0;
        onlyRun = 1;
        mazeSess = 1;
        figureState = 0;
    elseif nargin == 3
        onlyRun = 1;
        mazeSess = 1;
        figureState = 0;
    elseif nargin == 4
        mazeSess = 1;
        figureState = 0;
    elseif nargin == 5
        figureState = 0;
    elseif nargin > 6
        disp('Too many input arguments.');
        return;
    end
    
    GlobalConst;
    
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    fileNamePeakFR = [fileName '_PeakFRAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];        
    if(fileState == 0)
        fileNameFull = [fileName '.mat'];
    else
        fileNameFull = [fileName '_convSpikesAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    end
    
    fullPath = [path fileNameFull];
    if(exist(fullPath,'file') == 0)
        if(fileState == 0)
            disp('File does not exist.');
        else
            disp(['The firing profile file does not exist. Try to run the',...
                    'function again with fileState = 0.']);
        end
        return;
    end
    load(fullPath,'filteredSpikeArrayRun','filteredSpikeArrayCue','timeStep');
    timeStepRun = timeStep;
    
    fullPath = [path fileName '_Info.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to run file does not exist');
        return;
    end
    load(fullPath,'beh');
    
    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
            
    GlobalConst;
    %% changed by Yingxue on 2/14/2021
    if(isfield(beh,'indStimLap'))
        ind = find(beh.mazeSess == mazeSess);
        firstStimTr = find(beh.indStimLap(ind) == 1 & beh.pulseMethod(ind) > 0,1,'first');
        if(~isempty(firstStimTr))
            indNonStim = ind(1:firstStimTr-1); 
        else
            indNonStim = ind;
        end
        trialNoNonStimGood = intersect(indNonStim,ind(behPar.indTrBadBeh == 0));
        [~,trialNoNonStimGood] = intersect(ind,trialNoNonStimGood);
        trialNoNonStimBad = intersect(indNonStim,ind(behPar.indTrBadBeh == 1));
        [~,trialNoNonStimBad] = intersect(ind,trialNoNonStimBad);
        pulseMeth = unique(beh.pulseMethod(beh.pulseMethod > 0));
        trialNoStim = [];
        trialNoStimCtrl = [];
        for i = 1:length(pulseMeth)
            firstStimTr = find(beh.indStimLap(ind) == 1 & beh.pulseMethod(ind) == pulseMeth(i),1,'first');
            lastStimTr = find(beh.indStimLap(ind) == 1 & beh.pulseMethod(ind) == pulseMeth(i),1,'last');
            indStim = ind(firstStimTr:lastStimTr);
            trialNoStim{i} = indStim(beh.indStimLap(indStim) == 1);
            trialNoStimCtrl{i} = indStim(beh.indStimLap(indStim) == 0);
            [~,trialNoStim{i}] = intersect(ind,trialNoStim{i});
            [~,trialNoStimCtrl{i}] = intersect(ind,trialNoStimCtrl{i});
        end
    else
        trialNoNonStimGood = find(behPar.indTrBadBeh == 0);
        trialNoNonStimBad = find(behPar.indTrBadBeh == 1);
        trialNoStim = [];
        trialNoStimCtrl = [];
    end
    %%
     
    %%%%%%%%%% calculate the peak firing rate
    neuronNo = length(filteredSpikeArrayRun);
    numSamples = length(timeStepRun);
    disp('calculate peak firing rate for non-stimulated good trials - run onset'); 
    pFRNonStimGoodStruct = PeakFRAligned(filteredSpikeArrayRun,trialNoNonStimGood,neuronNo,numSamples);
    disp('calculate peak firing rate for non-stimulated bad trials - run onset');
    pFRNonStimBadStruct = PeakFRAligned(filteredSpikeArrayRun,trialNoNonStimBad,neuronNo,numSamples);
    disp('calculate peak firing rate for stimulated trials - run onset');
    %% changed by Yingxue on 2/14/2021
    pFRStimStruct = [];
    pFRStimCtrlStruct = [];
    for i = 1:length(pulseMeth)
        pFRStimStruct{i} = PeakFRAligned(filteredSpikeArrayRun,trialNoStim{i},neuronNo,numSamples);
       
        %% added by Yingxue on 2/14/2021
        pFRStimCtrlStruct{i} = PeakFRAligned(filteredSpikeArrayRun,trialNoStimCtrl{i},neuronNo,numSamples);
    end
    
    %% added by Yingxue on 1/21/2022
    disp('calculate peak firing rate for non-stimulated good trials - cue'); 
    pFRNonStimCueGoodStruct = PeakFRAligned(filteredSpikeArrayCue,trialNoNonStimGood,neuronNo,numSamples);
    disp('calculate peak firing rate for non-stimulated bad trials - cue');
    pFRNonStimCueBadStruct = PeakFRAligned(filteredSpikeArrayCue,trialNoNonStimBad,neuronNo,numSamples);
    disp('calculate peak firing rate for stimulated trials - cue');
    %% changed by Yingxue on 2/14/2021
    pFRStimCueStruct = [];
    pFRStimCtrlCueStruct = [];
    for i = 1:length(pulseMeth)
        pFRStimCueStruct{i} = PeakFRAligned(filteredSpikeArrayCue,trialNoStim{i},neuronNo,numSamples);
       
        %% added by Yingxue on 2/14/2021
        pFRStimCtrlCueStruct{i} = PeakFRAligned(filteredSpikeArrayCue,trialNoStimCtrl{i},neuronNo,numSamples);
    end
    %%
        
    save([path fileNamePeakFR], 'pFRNonStimGoodStruct','pFRNonStimBadStruct','pFRStimStruct','pFRStimCtrlStruct',...
        'pFRNonStimCueGoodStruct','pFRNonStimCueBadStruct','pFRStimCueStruct','pFRStimCtrlCueStruct',...
        'trialNoNonStimGood','trialNoNonStimBad','trialNoStim','trialNoStimCtrl','pulseMeth','-v7.3');
                       
    %%%%%%%%% draw figure is the state is on
    if(figureState ~= 0)
        % Ensure root units are pixels and get the size of the screen and create a
        % figure window
        set(0,'Units','pixels') 
        
        %%%% plot peak and mean instantaneous firing rate
        plotPFR(neuronNo,pFRNonStimGoodStruct.peakFR,pFRNonStimGoodStruct.meanInstFR);
        title('All neurons, non-stimulated good trials')
        
        %%%% plot neurons peak instantaneous firing rate vs mean firing rate
        plotPFRVsMInstFR(pFRNonStimGoodStruct.peakFR,pFRNonStimGoodStruct.meanInstFR);
        title('All neurons, non-stimulated good trials')
        
        %%%% plot neurons peak/mean instantaneous firing rate ratio vs mean
        %%%% firing rate
%         plotMInstFRVsP2M(pFRNonStimStruct.meanInstFR,pFRNonStimStruct.p2MInstRatio);
%         title('All neurons, non-stimulated trials')

        %%%% plot peak and mean instantaneous firing rate
        plotPFR(neuronNo,pFRStimStruct.peakFR,pFRStimStruct.meanInstFR);
        title('All neurons, stimulated trials')
        
        %%%% plot neurons peak instantaneous firing rate vs mean firing rate
        plotPFRVsMInstFR(pFRStimStruct.peakFR,pFRStimStruct.meanInstFR);
        title('All neurons, stimulated trials')
       
    end
    
    clear mydata;
    
end
