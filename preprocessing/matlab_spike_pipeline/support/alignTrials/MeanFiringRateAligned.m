function MeanFiringRateAligned(path,fileName,onlyRun,mazeSess,figureState)
% Calculate the mean firing rate for each recorded neuron
% path:         the path of the recording file
% fileName:     name of the recording file
% onlyRun:      1: only consider the time period when the animal is running 
% figureState:  0: figure off
%               1: plot the mean and std of the mean firing rate of each
%                  neuron
%               2: plot the histogram of mean firing rate
%
% Example:
% MeanFiringRate('./','A111-20150301-01_DataStructure_mazeSection1_TrialType1',1,1,0)

    %%%%%%%%% check arguments
    if nargin<2
        disp('At least two arguments are needed for this function.');
        return;
    elseif nargin == 2
        onlyRun = 1;
        mazeSess = 1;
        figureState = 0;
    elseif nargin == 3
        figureState = 0;
        mazeSess = 1;
    elseif nargin == 4
        figureState = 0;
    elseif nargin > 5
        disp('Too many input arguments.');
        return;
    end
        
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    fileNameFR = [fileName '_FRAlignedRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    
    fullPath = [path fileName '_alignedSpikesPerNPerT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _alignedSpikesPerNPerT file does not exist');
        return;
    end
    load(fullPath,'trialsRunSpikes');
    
    fullPath = [path fileName '_alignRun_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to run file does not exist');
        return;
    end
    load(fullPath,'trialsRun');
    
    fullPath = [path fileName '_PeakFRAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];    
    if(exist(fullPath) == 0)
        disp('The _PeakFRAligned file does not exist');
        return;
    end
    load(fullPath,'trialNoNonStimGood','trialNoNonStimBad','trialNoStim','trialNoStimCtrl','pulseMeth');
    
    fullPath = [path fileName '_Info.mat']; 
    if(exist(fullPath) == 0)
        disp('The behavioral information file does not exist');
        return;
    end
    load(fullPath,'rec');
    
    %%%%%%%%% initialize constants
    GlobalConst;    
    
    %%%%%%%%% calculate firing rate for each neuron over each trial
    %%% extract the spikes from trials data structure  
    spikes = trialsRunSpikes.Time;
    
    disp('Calculate mean firing rate for non-stimulated good laps');
    mFRStructNonStimGood = MFRAligned(spikes,trialNoNonStimGood,rec.numNeurons,...
                    trialsRun.numSamples/sampleFq,[]); 
    
    disp('Calculate mean firing rate for non-stimulated bad laps');
    mFRStructNonStimBad = MFRAligned(spikes,trialNoNonStimBad,rec.numNeurons,...
                    trialsRun.numSamples/sampleFq,[]); 
         
    %% added by Yingxue on 7/14/2021
    mFRStructStim = [];
    mFRStructStimCtrl = [];
    for i = 1:length(pulseMeth)
        disp('Calculate mean firing rate for stimulated laps');
        mFRStructStim{i} = MFRAligned(spikes,trialNoStim{i},rec.numNeurons,...
                        trialsRun.numSamples/sampleFq,[]); 

        disp('Calculate mean firing rate for stimulated control laps');
        mFRStructStimCtrl{i} = MFRAligned(spikes,trialNoStimCtrl{i},rec.numNeurons,...
                        trialsRun.numSamples/sampleFq,[]); 
    end
    %% 
    
    save([path fileNameFR], 'mFRStructNonStimGood','mFRStructNonStimBad',...
        'mFRStructStim','mFRStructStimCtrl');
