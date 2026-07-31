function PeakFiringRate_smTr_GoodTr(path,fileName,spaceBin,fileState,onlyRun)
% Calculate the smoothed mean firing rate curve over trials and the peak
% firing rate for each recorded neuron, separating good vs bad trials
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
%
% Example: 
% PeakFiringRate_smTr_GoodTr('./','A111-20150301-01_DataStructure_mazeSection1_TrialType1',1,1,1,1)

    %%%%%%%%% check arguments
    if nargin<2
        disp('At least two arguments are needed for this function.');
        return;
    elseif nargin == 2
        spaceBin = 20; %(in mm)
        fileState = 0;
        onlyRun = 1;
    elseif nargin == 3
        fileState = 0;
        onlyRun = 1;
    elseif nargin == 4
        onlyRun = 1;
    elseif nargin > 5
        disp('Too many input arguments.');
        return;
    end
    
    GlobalConst;
    
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    fileNamePeakFR = [fileName '_PeakFR_GoodTr' num2str(spaceBin) ...
                      'mm_Run' num2str(onlyRun) '.mat'];        
   
    fileNameFull = [fileName '_convSpikesDist' num2str(spaceBin) ...
                        'mm_Run' num2str(onlyRun) '.mat'];   
    fullPath = [path fileNameFull];
    if(exist(fullPath,'file') == 0)
        disp(['The firing profile file does not exist. Try to run the',...
                    'function again with fileState = 0.']);
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
    numSamples = length(paramC.spaceSteps{mazeSess});
    
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
    disp('calculate peak firing rate for each session')
    pFRStructSessGoodTr = cell(length(mazeSess),1);
    pFRStructNormTSessGoodTr = cell(length(mazeSess),1);
    pFRStructSessBadTr = cell(length(mazeSess),1);
    pFRStructNormTSessBadTr = cell(length(mazeSess),1);
    pFRStructSessOKTr = cell(length(mazeSess),1);
    pFRStructNormTSessOKTr = cell(length(mazeSess),1);
    if(length(mazeSess)>1)
        for i = 1:length(mazeSess) 
            fprintf('\nSession %d\n',i);
            %%% calculate peak firing rate for good, bad and ok trials
            %%% within each session
            %%% added by Yingxue on 12/16/2020
            indLaps = find(beh.mazeSess == mazeSess(i));
            indLaps = intersect(indLaps,beh.indTrCtrl); 
            indLapsGoodTr = intersect(indLaps,beh.indGoodTrCtrl); 
            pFRStructSessGoodTr{i} = PeakFR(filteredSpikeArray,indLapsGoodTr,...
                        rec.numNeurons,max(numSamples));
            pFRStructNormTSessGoodTr{i} = PeakFR(filteredSpikeArrayNormT,indLapsGoodTr,...
                        rec.numNeurons,max(numSamples)); 
            indLapsBadTr = intersect(indLaps,beh.indBadTrCtrl); 
            pFRStructSessBadTr{i} = PeakFR(filteredSpikeArray,indLapsBadTr,...
                        rec.numNeurons,max(numSamples));
            pFRStructNormTSessBadTr{i} = PeakFR(filteredSpikeArrayNormT,indLapsBadTr,...
                        rec.numNeurons,max(numSamples));
            indLapsOKTr = setdiff(indLaps,[beh.indGoodTrCtrl,beh.indBadTrCtrl]); 
            pFRStructSessOKTr{i} = PeakFR(filteredSpikeArray,indLapsOKTr,...
                        rec.numNeurons,max(numSamples));
            pFRStructNormTSessOKTr{i} = PeakFR(filteredSpikeArrayNormT,indLapsOKTr,...
                        rec.numNeurons,max(numSamples));
        end
    else
        indLapsGoodTr = beh.indGoodTrCtrl; 
        pFRStructSessGoodTr{1} = PeakFR(filteredSpikeArray,indLapsGoodTr,...
                    rec.numNeurons,max(numSamples));
        pFRStructNormTSessGoodTr{1} = PeakFR(filteredSpikeArrayNormT,indLapsGoodTr,...
                    rec.numNeurons,max(numSamples)); 
        indLapsBadTr = beh.indBadTrCtrl; 
        pFRStructSessBadTr{1} = PeakFR(filteredSpikeArray,indLapsBadTr,...
                    rec.numNeurons,max(numSamples));
        pFRStructNormTSessBadTr{1} = PeakFR(filteredSpikeArrayNormT,indLapsBadTr,...
                    rec.numNeurons,max(numSamples));
        indLapsOKTr = setdiff(beh.indTrCtrl,[beh.indGoodTrCtrl,beh.indBadTrCtrl]);
        pFRStructSessOKTr{1} = PeakFR(filteredSpikeArray,indLapsOKTr,...
                    rec.numNeurons,max(numSamples));
        pFRStructNormTSessOKTr{1} = PeakFR(filteredSpikeArrayNormT,indLapsOKTr,...
                    rec.numNeurons,max(numSamples));
    end
     
    save([path fileNamePeakFR],...
            'pFRStructSessGoodTr','pFRStructNormTSessGoodTr',...
            'pFRStructSessBadTr','pFRStructNormTSessBadTr',...
            'pFRStructSessOKTr','pFRStructNormTSessOKTr',...
            '-v7.3');
                       
    clear mydata;
    
end
