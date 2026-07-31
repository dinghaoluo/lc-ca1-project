function PeakFiringRate_AlignedRunOnset2P(path,fileName,fileState,onlyRun,mazeSess,figureState)
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
    fileNamePeakFR = [fileName '_PeakFR_msess' num2str(mazeSess) '_RunOnset' num2str(onlyRun) '.mat'];        
    if(fileState == 0)
        fileNameFull = [fileName '.mat'];
    else
        fileNameFull = [fileName '_convSpikesAligned_msess' num2str(mazeSess) '_BefRun' num2str(onlyRun) '.mat'];
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
    load(fullPath);
    
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
    
    fullPath = [path fileName '_PeakFRAlignedRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];    
    if(exist(fullPath) == 0)
        disp('The _PeakFRAligned file does not exist');
        return;
    end
    load(fullPath,'trialNoNonStimGood','trialNoNonStimBad','trialNoStim','trialNoStimCtrl','pulseMeth');
            
    GlobalConst2P;
         
    %%%%%%%%%% calculate the peak firing rate
    neuronNo = length(filteredSpikeArrayRunOnset);
    numSamples = size(filteredSpikeArrayRunOnset{end},2);
    disp('calculate peak firing rate for non-stimulated good trials');
    pFRNonStimGoodStruct = PeakFRAligned2P(filteredSpikeArrayRunOnset,trialNoNonStimGood,neuronNo,numSamples,startTrNo);
    pdFFNonStimGoodStruct = PeakFRAligned2P(dFFArrayRunOnset,trialNoNonStimGood,neuronNo,numSamples,startTrNo);
    disp('calculate peak firing rate for non-stimulated bad trials');
    pFRNonStimBadStruct = PeakFRAligned2P(filteredSpikeArrayRunOnset,trialNoNonStimBad,neuronNo,numSamples,startTrNo);
    pdFFNonStimBadStruct = PeakFRAligned2P(dFFArrayRunOnset,trialNoNonStimBad,neuronNo,numSamples,startTrNo);
       
    %% changed by Yingxue on 7/15/2021
    pFRStimStruct = [];
    pFRStimCtrlStruct = [];
    pdFFStimStruct = [];
    pdFFStimCtrlStruct = [];
    for i = 1:length(pulseMeth)
        disp('calculate peak firing rate for stimulated trials');
        pFRStimStruct{i} = PeakFRAligned2P(filteredSpikeArrayRunOnset,trialNoStim{i},neuronNo,numSamples,startTrNo);
        pdFFStimStruct{i} = PeakFRAligned2P(dFFArrayRunOnset,trialNoStim{i},neuronNo,numSamples,startTrNo);
        
        disp('calculate peak firing rate for stimulated ctrl trials');
        pFRStimCtrlStruct{i} = PeakFRAligned2P(filteredSpikeArrayRunOnset,trialNoStimCtrl{i},neuronNo,numSamples,startTrNo);
        pdFFStimCtrlStruct{i} = PeakFRAligned2P(dFFArrayRunOnset,trialNoStimCtrl{i},neuronNo,numSamples,startTrNo);
    end
        
    save([path fileNamePeakFR], 'pFRNonStimGoodStruct','pFRNonStimBadStruct','pFRStimStruct',...
        'pFRStimCtrlStruct','pdFFNonStimGoodStruct','pdFFNonStimBadStruct','pdFFStimStruct','pdFFStimCtrlStruct',...
        '-v7.3');
                       
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
