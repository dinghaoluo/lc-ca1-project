function MeanFiringRateCtrlOnly(path,fileName,onlyRun,mazeSess)
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
% MeanFiringRate('./','A111-20150301-01_DataStructure_mazeSection1_TrialType1',1,0)

    %%%%%%%%% check arguments
    if nargin<2
        disp('At least two arguments are needed for this function.');
        return;
    elseif nargin == 2
        onlyRun = 1;
        mazeSess = 1;
    elseif nargin == 3
        mazeSess = 1;
    elseif nargin > 4
        disp('Too many input arguments.');
        return;
    end
        
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    fileNameFR = [fileName '_FR_Ctrl_Run' num2str(onlyRun) '_mazeSess' num2str(mazeSess) '.mat'];
    fileNameExt = [fileName '_ext.mat'];
    fileNameInfo = [fileName '_Info.mat'];
    fileName = [fileName '.mat'];
    
    if(onlyRun == 0)
        fullPath = [path fileName];
        if(exist(fullPath) == 0)
            disp('File does not exist.');
            return;
        end
        load(fullPath);
    else
        fullPath = [path fileName];
        if(exist(fullPath) == 0)
            disp('File does not exist.');
            return;
        end
        load(fullPath,'lapList');
        
        fullPath = [path fileNameExt];
        if(exist(fullPath) == 0)
            disp(['Extended file does not exist.',...
                'Please run function SpikeDuringRun first']);
            return;
        end
        load(fullPath);
    end
    
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo_smTr(path,fileName);
    end
    load(fullPath);
   
    %%%%%%%%% initialize constants
    GlobalConst;    
    
    %%%%%%%%% calculate firing rate for each neuron over each trial
    %%% extract the spikes from trials data structure  
    if(onlyRun == 0)
        spikes = getRecField(trials,'spikes',1:length(lapList));
        indRunInLap = [];
    else
        spikes = getRecField(trialsExt,'spikes',1:length(lapList));
        indRunInLap = beh.indRunInLap;
    end
    
    indLaps = find(beh.mazeSess == mazeSess);
    indLapsCtrl = intersect(indLaps,beh.indTrCtrl);

    mFRStructSessCtrl = MFR(spikes,indLapsCtrl,rec.numNeurons,...
                beh.lenTrials/sampleFq,indRunInLap);
    
    save([path fileNameFR], 'mFRStructSessCtrl');
       
    clear trials spikes rec beh indRunInLap indLaps indLapsStim indLapsCtrl...
        mFRStructSessCtrl 