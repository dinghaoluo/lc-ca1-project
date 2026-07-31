function PeakFiringRate_smTrCtrlOnly(path,fileName,spaceBin,fileState,onlyRun,mazeSess)
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
        mazeSess = 1;
    elseif nargin == 3
        fileState = 0;
        onlyRun = 1;
        mazeSess = 1;
    elseif nargin == 4
        onlyRun = 1;
        mazeSess = 1;
    elseif nargin == 5
        mazeSess = 1;
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
    fileNameFR = [fileName '_FR_Ctrl_Run' num2str(onlyRun) '_mazeSess' num2str(mazeSess) '.mat'];
    fileNamePeakFR = [fileName '_PeakFRCtrl' num2str(spaceBin) ...
                      'mm_Run' num2str(onlyRun) '_mazeSess' num2str(mazeSess) '.mat'];        
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
    
    fullPath = [path fileNameFR];
    if(exist(fullPath,'file') == 0)
        disp('_FR_Ctrl file does not exist. Please run MeanFiringRateCtrlOnly function first');
    end
    load(fullPath,'mFRStructSessCtrl');
        
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
    indLaps = mFRStructSessCtrl.indLapList;
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
    
        
    %%%%%%%%%% calculate the peak firing rate                 
    disp('calculate peak firing rate for each session')
    pFRStructSessCtrl = PeakFR(filteredSpikeArray,indLaps,...
                    rec.numNeurons,max(numSamples));
    pFRStructNormTSessCtrl = PeakFR(filteredSpikeArrayNormT,indLaps,...
                        rec.numNeurons,max(numSamples)); 
                 
    save([path fileNamePeakFR],...
            'pFRStructSessCtrl','pFRStructNormTSessCtrl',...
            '-v7.3');
                       
    clear mydata;
    
end
