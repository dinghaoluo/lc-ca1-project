function PeakFiringRate_AlignedCtrlOnly(path,fileName,fileState,onlyRun,mazeSess,figureState)
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
    fileNamePeakFR = [fileName '_PeakFRAlignedCtrl_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];        
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
            disp(['The firing profile file does not exist. Try to run the ',...
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
        [~,trialNoNonStim] = intersect(ind,indNonStim);
    else
        trialNoNonStim = 1:length(behPar.indTrBadBeh);
    end
    %%
     
    %%%%%%%%%% calculate the peak firing rate
    neuronNo = length(filteredSpikeArrayRun);
    numSamples = length(timeStepRun);
    disp('calculate peak firing rate for non-stimulated trials - run onset'); 
    pFRNonStimStruct = PeakFRAligned(filteredSpikeArrayRun,trialNoNonStim,neuronNo,numSamples);
    
    disp('calculate peak firing rate for non-stimulated trials - cue'); 
    pFRNonStimCueStruct = PeakFRAligned(filteredSpikeArrayCue,trialNoNonStim,neuronNo,numSamples);
            
    save([path fileNamePeakFR], 'pFRNonStimStruct','pFRNonStimCueStruct','trialNoNonStim','-v7.3');
        
    clear mydata;
    
end
