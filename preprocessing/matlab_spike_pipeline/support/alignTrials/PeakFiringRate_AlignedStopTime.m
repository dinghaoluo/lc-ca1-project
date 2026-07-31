function PeakFiringRate_AlignedStopTime(path,fileName,fileState,onlyRun,mazeSess,figureState)
% Calculate the smoothed mean firing rate curve over trials and the peak
% firing rate for each recorded neuron, cluster the trials based on stop
% time before run
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
% PeakFiringRate_AlignedStopTime('./','A111-20150301-01_DataStructure_mazeSection1_TrialType1',1,1,1,1)

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
    fileNamePeakFR = [fileName '_PeakFRAlignedStopTime_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];        
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
    load(fullPath,'filteredSpikeArrayRun','timeStep');
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
        indST1 = behPar.rewardToRun < stopTimeRew2RunThrL & behPar.indTrBadBeh == 0;
        trialNoNonStimST1 = intersect(indNonStim,ind(indST1));
        [~,trialNoNonStimST1] = intersect(ind,trialNoNonStimST1);
        indST2 = behPar.rewardToRun >= stopTimeRew2RunThrL &...
            behPar.rewardToRun < stopTimeRew2RunThrH & behPar.indTrBadBeh == 0;
        trialNoNonStimST2 = intersect(indNonStim,ind(indST2));
        [~,trialNoNonStimST2] = intersect(ind,trialNoNonStimST2);
        indST3 = behPar.rewardToRun >= stopTimeRew2RunThrH & behPar.indTrBadBeh == 0;
        trialNoNonStimST3 = intersect(indNonStim,ind(indST3));
        [~,trialNoNonStimST3] = intersect(ind,trialNoNonStimST3);
    else
        trialNoNonStimST1 = find(behPar.rewardToRun < stopTimeRew2RunThrL ...
            && behPar.indTrBadBeh == 0);
        trialNoNonStimST2 = find(behPar.rewardToRun >= stopTimeRew2RunThrL &&...
            behPar.rewardToRun < stopTimeRew2RunThrH && behPar.indTrBadBeh == 0);
        trialNoNonStimST2 = find(behPar.rewardToRun >= stopTimeRew2RunThrH ...
            && behPar.indTrBadBeh == 0);
    end
    %%
     
    %%%%%%%%%% calculate the peak firing rate
    neuronNo = length(filteredSpikeArrayRun);
    numSamples = length(timeStepRun);
    disp('calculate peak firing rate for trials with short stop time'); 
    pFRNonStimST1Struct = PeakFRAligned(filteredSpikeArrayRun,trialNoNonStimST1,neuronNo,numSamples);
    disp('calculate peak firing rate for trials with intermediate stop time');
    pFRNonStimST2Struct = PeakFRAligned(filteredSpikeArrayRun,trialNoNonStimST2,neuronNo,numSamples);
    disp('calculate peak firing rate for trials with long stop time');
    pFRNonStimST3Struct = PeakFRAligned(filteredSpikeArrayRun,trialNoNonStimST3,neuronNo,numSamples);
            
    save([path fileNamePeakFR], 'pFRNonStimST1Struct','pFRNonStimST2Struct','pFRNonStimST3Struct',...
        'trialNoNonStimST1','trialNoNonStimST2','trialNoNonStimST3','-v7.3');

    clear mydata;
    
end
