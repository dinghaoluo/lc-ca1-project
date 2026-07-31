function PeakFiringRate_Aligned2P(path,fileName,fileState,onlyRun,mazeSess,figureState)
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
    
    GlobalConst2P;
    
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
               
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
            
    GlobalConst2P;
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

    
    %% peak firing rate aligned to run
    if(fileState == 0)
        fileNameFull = [fileName '.mat'];
    else
        fileNameFull = [fileName '_convSpikesAlignedRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
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
    load(fullPath,'filteredSpikeArrayRun','dFFArrayRun');
    
    %% peak firing rate aligned to run
    if(fileState == 0)
        fileNameFull = [fileName '.mat'];
    else
        fileNameFull = [fileName '_convSpikesAlignedCue_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
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
    load(fullPath,'filteredSpikeArrayCue');
    
    %%%%%%%%%% calculate the peak firing rate (align to run)
    neuronNo = length(filteredSpikeArrayRun);
    numSamples = size(filteredSpikeArrayRun{end},2);
    disp('Align to run onset')
    disp('calculate peak firing rate for non-stimulated good trials'); 
    pFRNonStimGoodStruct = PeakFRAligned2P(filteredSpikeArrayRun,trialNoNonStimGood,neuronNo,numSamples,startTrNo);
    pdFFNonStimGoodStruct = PeakFRAligned2P(dFFArrayRun,trialNoNonStimGood,neuronNo,numSamples,startTrNo);
    disp('calculate peak firing rate for non-stimulated bad trials');
    pFRNonStimBadStruct = PeakFRAligned2P(filteredSpikeArrayRun,trialNoNonStimBad,neuronNo,numSamples,startTrNo);
    pdFFNonStimBadStruct = PeakFRAligned2P(dFFArrayRun,trialNoNonStimBad,neuronNo,numSamples,startTrNo);
    disp('calculate peak firing rate for stimulated trials');
    %% changed by Yingxue on 2/14/2021
    pFRStimStruct = [];
    pFRStimCtrlStruct = [];
    pdFFStimStruct = [];
    pdFFStimCtrlStruct = [];
    for i = 1:length(pulseMeth)
        pFRStimStruct{i} = PeakFRAligned2P(filteredSpikeArrayRun,trialNoStim{i},neuronNo,numSamples,startTrNo);
        pdFFStimStruct{i} = PeakFRAligned2P(dFFArrayRun,trialNoStim{i},neuronNo,numSamples,startTrNo);
        
        %% added by Yingxue on 2/14/2021
        pFRStimCtrlStruct{i} = PeakFRAligned2P(filteredSpikeArrayRun,trialNoStimCtrl{i},neuronNo,numSamples,startTrNo);
        pdFFStimCtrlStruct{i} = PeakFRAligned2P(dFFArrayRun,trialNoStimCtrl{i},neuronNo,numSamples,startTrNo);
    end
        
    fileNamePeakFR = [fileName '_PeakFRAlignedRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];    
    save([path fileNamePeakFR], 'pFRNonStimGoodStruct','pFRNonStimBadStruct','pFRStimStruct','pFRStimCtrlStruct',...
        'pdFFNonStimGoodStruct','pdFFNonStimBadStruct','pdFFStimStruct','pdFFStimCtrlStruct',...
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
    
    clear filteredSpikeArrayRun dFFArrayRun ...
        pFRNonStimGoodStruct pFRNonStimBadStruct pFRStimStruct pFRStimCtrlStruct ...
        pdFFNonStimGoodStruct pdFFNonStimBadStruct pdFFStimStruct pdFFStimCtrlStruct
    
    
    %% peak firing rate aligned to cue   
    if(fileState == 0)
        fileNameFull = [fileName '.mat'];
    else
        fileNameFull = [fileName '_convSpikesAlignedCue_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
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
    load(fullPath,'filteredSpikeArrayCue','dFFArrayCue');
    
    %%%%%%%%%% calculate the peak firing rate (aligned to cue)
    neuronNo = length(filteredSpikeArrayCue);
    numSamples = size(filteredSpikeArrayCue{end},2);
    disp('Align to cue')
    disp('calculate peak firing rate for non-stimulated good trials'); 
    pFRNonStimGoodStructCue = PeakFRAligned2P(filteredSpikeArrayCue,trialNoNonStimGood,neuronNo,numSamples,startTrNo);
    pdFFNonStimGoodStructCue = PeakFRAligned2P(dFFArrayCue,trialNoNonStimGood,neuronNo,numSamples,startTrNo);
    disp('calculate peak firing rate for non-stimulated bad trials');
    pFRNonStimBadStructCue = PeakFRAligned2P(filteredSpikeArrayCue,trialNoNonStimBad,neuronNo,numSamples,startTrNo);
    pdFFNonStimBadStructCue = PeakFRAligned2P(dFFArrayCue,trialNoNonStimBad,neuronNo,numSamples,startTrNo);
    disp('calculate peak firing rate for stimulated trials');
    %% changed by Yingxue on 2/14/2021
    pFRStimStructCue = [];
    pFRStimCtrlStructCue = [];
    pdFFStimStructCue = [];
    pdFFStimCtrlStructCue = [];
    for i = 1:length(pulseMeth)
        pFRStimStructCue{i} = PeakFRAligned2P(filteredSpikeArrayCue,trialNoStim{i},neuronNo,numSamples,startTrNo);
        pdFFStimStructCue{i} = PeakFRAligned2P(dFFArrayCue,trialNoStim{i},neuronNo,numSamples,startTrNo);
        
        %% added by Yingxue on 2/14/2021
        pFRStimCtrlStructCue{i} = PeakFRAligned2P(filteredSpikeArrayCue,trialNoStimCtrl{i},neuronNo,numSamples,startTrNo);
        pdFFStimCtrlStructCue{i} = PeakFRAligned2P(dFFArrayCue,trialNoStimCtrl{i},neuronNo,numSamples,startTrNo);
    end
    fileNamePeakFR = [fileName '_PeakFRAlignedCue_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];    
    save([path fileNamePeakFR], 'pFRNonStimGoodStructCue','pFRNonStimBadStructCue','pFRStimStructCue','pFRStimCtrlStructCue',...
        'pdFFNonStimGoodStructCue','pdFFNonStimBadStructCue','pdFFStimStructCue','pdFFStimCtrlStructCue','-v7.3');
    clear filteredSpikeArrayCue dFFArrayCue ...
        pFRNonStimGoodStructCue pFRNonStimBadStructCue pFRStimStructCue pFRStimCtrlStructCue ...
        pdFFNonStimGoodStructCue pdFFNonStimBadStructCue pdFFStimStructCue pdFFStimCtrlStructCue 
    
    
    %% peak firing rate aligned to reward   
    if(fileState == 0)
        fileNameFull = [fileName '.mat'];
    else
        fileNameFull = [fileName '_convSpikesAlignedRew_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
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
    load(fullPath,'filteredSpikeArrayRew','dFFArrayRew');
    
    %%%%%%%%%% calculate the peak firing rate (aligned to reward)
    neuronNo = length(filteredSpikeArrayRew);
    numSamples = size(filteredSpikeArrayRew{end},2);
    disp('Align to reward')
    disp('calculate peak firing rate for non-stimulated good trials'); 
    pFRNonStimGoodStructRew = PeakFRAligned2P(filteredSpikeArrayRew,trialNoNonStimGood,neuronNo,numSamples,startTrNo);
    pdFFNonStimGoodStructRew = PeakFRAligned2P(dFFArrayRew,trialNoNonStimGood,neuronNo,numSamples,startTrNo);
    disp('calculate peak firing rate for non-stimulated bad trials');
    pFRNonStimBadStructRew = PeakFRAligned2P(filteredSpikeArrayRew,trialNoNonStimBad,neuronNo,numSamples,startTrNo);
    pdFFNonStimBadStructRew = PeakFRAligned2P(dFFArrayRew,trialNoNonStimBad,neuronNo,numSamples,startTrNo);
    disp('calculate peak firing rate for stimulated trials');
    %% changed by Yingxue on 2/14/2021
    pFRStimStructRew = [];
    pFRStimCtrlStructRew = [];
    pdFFStimStructRew = [];
    pdFFStimCtrlStructRew = [];
    for i = 1:length(pulseMeth)
        pFRStimStructRew{i} = PeakFRAligned2P(filteredSpikeArrayRew,trialNoStim{i},neuronNo,numSamples,startTrNo);
        pdFFStimStructRew{i} = PeakFRAligned2P(dFFArrayRew,trialNoStim{i},neuronNo,numSamples,startTrNo);
        
        %% added by Yingxue on 2/14/2021
        pFRStimCtrlStructRew{i} = PeakFRAligned2P(filteredSpikeArrayRew,trialNoStimCtrl{i},neuronNo,numSamples,startTrNo);
        pdFFStimCtrlStructRew{i} = PeakFRAligned2P(dFFArrayRew,trialNoStimCtrl{i},neuronNo,numSamples,startTrNo);
    end
    fileNamePeakFR = [fileName '_PeakFRAlignedRew_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];        
    save([path fileNamePeakFR], 'pFRNonStimGoodStructRew','pFRNonStimBadStructRew','pFRStimStructRew','pFRStimCtrlStructRew',...
        'pdFFNonStimGoodStructRew','pdFFNonStimBadStructRew','pdFFStimStructRew','pdFFStimCtrlStructRew','-v7.3');
    clear filteredSpikeArrayRew dFFArrayRew ...
        pFRNonStimGoodStructRew pFRNonStimBadStructRew pFRStimStructRew pFRStimCtrlStructRew ...
        pdFFNonStimGoodStructRew pdFFNonStimBadStructRew pdFFStimStructRew pdFFStimCtrlStructRew 
    
    clear mydata;
    
end
