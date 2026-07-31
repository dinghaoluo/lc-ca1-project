function PeakFiringRate_smTr(path,fileName,spaceBin,fileState,onlyRun,figureState)
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
        spaceBin = 20; %(in mm)
        fileState = 0;
        onlyRun = 1;
        figureState = 0;
    elseif nargin == 3
        fileState = 0;
        onlyRun = 1;
        figureState = 0;
    elseif nargin == 4
        onlyRun = 1;
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
    fileNamePeakFR = [fileName '_PeakFR' num2str(spaceBin) ...
                      'mm_Run' num2str(onlyRun) '.mat'];        
    if(fileState == 0)
        fileNameFull = [fileName '.mat'];
    else
        fileNameFull = [fileName '_convSpikesDist' num2str(spaceBin) ...
                        'mm_Run' num2str(onlyRun) '.mat'];
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
        
    %%%%%%%%% initialize constants
    if(isempty(indexFileName))
        fileNameInfo = [fileName '_Info.mat'];
    else
        fileNameInfo = [fileName(1:indexFileName(end)-1) '_Info.mat'];
    end 

    fullPath = [path fileNameInfo];
    if(exist(fullPath,'file') == 0)
        BasicInfo_smTr(path,fileName);
    end
    load(fullPath);
    mazeSess = beh.mazeSessAll;
    trackLen = unique(beh.trackLen);
    numSamples = zeros(1,length(trackLen));
    for i = 1:length(trackLen)
        numSamples(i) = length(paramC.spaceSteps{i});
    end
    
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
        [filteredSpikeArray, filteredSpikeArrayNormT] = ...
            ConvSpikeTrainDistPar_smTr(path,fileName,spaceBin,onlyRun,1); 
    end
    
    %%%%%%%%%% calculate the peak firing rate
    disp('calculate peak firing rate for all the trials');
    pFRStruct = PeakFR(filteredSpikeArray,beh.indGoodLap,rec.numNeurons,max(numSamples));
    disp('calculate peak firing rate for all the trials with time normalization');
    pFRStructNormT = PeakFR(filteredSpikeArrayNormT,beh.indGoodLap,...
                            rec.numNeurons,max(numSamples)); 
                        
    disp('calculate peak firing rate for each session')
    pFRStructSess = cell(length(mazeSess),1);
    pFRStructNormTSess = cell(length(mazeSess),1);  
    if(length(mazeSess)>1)
        for i = 1:length(mazeSess) 
            fprintf('\nSession %d\n',i);
            indLaps = find(beh.mazeSess == mazeSess(i)); 
            indLaps = intersect(indLaps,beh.indGoodLap);
            pFRStructSess{i} = PeakFR(filteredSpikeArray,indLaps,rec.numNeurons,...
                        max(numSamples));
            pFRStructNormTSess{i} = PeakFR(filteredSpikeArrayNormT,indLaps,...
                        rec.numNeurons,max(numSamples));   
        end
    end
     
    save([path fileNamePeakFR], 'pFRStruct','pFRStructNormT',...
            'pFRStructSess','pFRStructNormTSess',...
            '-v7.3');
                       
    %%%%%%%%% draw figure is the state is on
    if(figureState ~= 0)
        % Ensure root units are pixels and get the size of the screen and create a
        % figure window
        set(0,'Units','pixels') 
        
        %%%% plot peak and mean instantaneous firing rate
        plotPFR(rec.numNeurons,pFRStruct.peakFR,pFRStruct.meanInstFR);
        title('All neuron')
        
        %%%% plot neurons peak instantaneous firing rate vs mean firing rate
        plotPFRVsMInstFR(pFRStruct.peakFR,pFRStruct.meanInstFR);
        title('All neuron')
        
        %%%% plot neurons peak/mean instantaneous firing rate ratio vs mean
        %%%% firing rate
%         plotMInstFRVsP2M(pFRStruct.meanInstFR,pFRStruct.p2MInstRatio);
%         title('All neuron')
        
        %%%% plot peak and mean instantaneous firing rate
        plotPFR(rec.numNeurons,pFRStructNormT.peakFR,pFRStructNormT.meanInstFR);
        title('All neuron normT')
        
        %%%% plot neurons peak instantaneous firing rate vs mean firing rate
        plotPFRVsMInstFR(pFRStructNormT.peakFR,pFRStructNormT.meanInstFR);
        title('All neuron normT')
        
        %%%% plot neurons peak/mean instantaneous firing rate ratio vs mean
        %%%% firing rate
%         plotMInstFRVsP2M(pFRStructNormT.meanInstFR,pFRStructNormT.p2MInstRatio);
%         title('All neuron normT')
               
    end
    
    clear mydata;
    
end
