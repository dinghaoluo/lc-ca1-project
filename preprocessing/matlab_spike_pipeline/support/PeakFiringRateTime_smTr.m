function PeakFiringRateTime_smTr(path,fileName,timeBin,fileState,figureState)
% Calculate the smoothed mean firing rate curve over trials and the peak
% firing rate for each recorded neuron
%             if fileState == 1, function "ConvSpikeTrain" should be
%             executed first
% path:         the path of the recording file
% fileName:     name of the recording file
% timeBin:      2SD of the Gaussian filter used to obtain the firing rate
%               profile (in mm), default value is 10 mm
% fileState:    0: calculate using the recorded data (default)
%               1: load the firing rate profile if the file exists, and
%               calculate the peak firing rate from there
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
        timeBin = 0.2; %(in s)
        fileState = 0;
        figureState = 0;
    elseif nargin == 3
        fileState = 0;
        figureState = 0;
    elseif nargin == 4
        figureState = 0;
    elseif nargin > 5
        disp('Too many input arguments.');
        return;
    end
    
    GlobalConst;
    paramP.maxTrialLenLimit = 4.5; %s
    
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    fileNamePeakFR = [fileName '_PeakFR' num2str(timeBin*1000) ...
                      'ms.mat'];        
    if(fileState == 0)
        fileNameFull = [fileName '.mat'];
    else
        fileNameFull = [fileName '_convSpikesTime' num2str(timeBin*1000) ...
                        'ms.mat'];
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
    
    load([path fileName '.mat']);
    numTrials = length(trials);
    trialLenArr = [];
    for i = 1:numTrials
        trialLenArr = [trialLenArr trials{i}.Nsamples];
    end
    indGoodTrial = find(trialLenArr <= paramP.maxTrialLenLimit*sampleFq);
       
    %%%%%%%% calculate the normalizing factor
    %%%%%%%% gaussian filter is normalized to a total energy = 1. Thus in
    %%%%%%%% time domain, each spike translates into a value of sum(gaussFilt).
    %%%%%%%% To guarantee that the mean inst firing rate is closest to the 
    %%%%%%%% mean firing rate of the neuron,
    %%%%%%%% we divide the profile array by sum(gaussFilt) to estimate the
    %%%%%%%% mean inst firing at each sampling point, and then * sampleFq
    %%%%%%%% to estimate the mean inst firing rate per second
    
    %%%%%%%%% calculate the smoothed mean firing rate curve for each neuron
    if(fileState == 0)
        [filteredSpikeArray, filteredSpikeArrayNormAmp] = ...
            ConvSpikeTrainTimePar_smTr(path,fileName,timeBin,1); 
    end
    
    numTrials = length(filteredSpikeArray);
    numNeurons = size(filteredSpikeArray{1},1);
    numSamples = size(filteredSpikeArray{1},2);
    %%%%%%%%%% calculate the peak firing rate
    disp('calculate peak firing rate for all the trials');
    pFRStruct = PeakFR(filteredSpikeArray,1:numTrials,numNeurons,numSamples);
    disp('calculate peak firing rate for all the trials with time normalization');
    pFRStructNormAmp = PeakFR(filteredSpikeArrayNormAmp,1:numTrials,...
                            numNeurons,numSamples); 
                       
    pFRStructGoodTr = PeakFR(filteredSpikeArray,indGoodTrial,numNeurons,numSamples);
    disp('calculate peak firing rate for all the trials with time normalization');
    pFRStructNormAmpGoodTr = PeakFR(filteredSpikeArrayNormAmp,indGoodTrial,...
                            numNeurons,numSamples); 
                        
    save([path fileNamePeakFR], 'pFRStruct','pFRStructNormAmp',...
            'pFRStructGoodTr','pFRStructNormAmpGoodTr','paramP',...
            '-v7.3');
                       
    %%%%%%%%% draw figure is the state is on
    if(figureState ~= 0)
        % Ensure root units are pixels and get the size of the screen and create a
        % figure window
        set(0,'Units','pixels') 
        
        %%%% plot peak and mean instantaneous firing rate
        plotPFR(numNeurons,pFRStruct.peakFR,pFRStruct.meanInstFR);
        title('All neuron')
        
        %%%% plot neurons peak instantaneous firing rate vs mean firing rate
        plotPFRVsMInstFR(pFRStruct.peakFR,pFRStruct.meanInstFR);
        title('All neuron')
        
        %%%% plot neurons peak/mean instantaneous firing rate ratio vs mean
        %%%% firing rate
%         plotMInstFRVsP2M(pFRStruct.meanInstFR,pFRStruct.p2MInstRatio);
%         title('All neuron')
        
        %%%% plot peak and mean instantaneous firing rate
        plotPFR(numNeurons,pFRStructNormAmp.peakFR,pFRStructNormAmp.meanInstFR);
        title('All neuron normT')
        
        %%%% plot neurons peak instantaneous firing rate vs mean firing rate
        plotPFRVsMInstFR(pFRStructNormAmp.peakFR,pFRStructNormAmp.meanInstFR);
        title('All neuron normT')
        
        %%%% plot neurons peak/mean instantaneous firing rate ratio vs mean
        %%%% firing rate
%         plotMInstFRVsP2M(pFRStructNormT.meanInstFR,pFRStructNormT.p2MInstRatio);
%         title('All neuron normT')
               
    end
    
    clear mydata;
    
end
