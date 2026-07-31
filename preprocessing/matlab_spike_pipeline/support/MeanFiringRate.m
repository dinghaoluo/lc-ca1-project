function MeanFiringRate(path,fileName,onlyRun,figureState)
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
        figureState = 0;
    elseif nargin == 3
        figureState = 0;
    elseif nargin > 4
        disp('Too many input arguments.');
        return;
    end
        
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    fileNameFR = [fileName '_FR_Run' num2str(onlyRun) '.mat'];
    fileNameExt = [fileName '_ext.mat'];
    fileNameInfo = [fileName '_Info.mat'];
    fileName = [fileName '.mat'];
    
    %% changed on 1/22/2022 by Yingxue, load files depending on onlyRun
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
    mazeSess = beh.mazeSessAll;
    
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
    disp('Calculate mean firing rate for all the good laps');
    mFRStruct = MFR(spikes,beh.indGoodLap,rec.numNeurons,...
                    beh.lenTrials/sampleFq,indRunInLap); 
    disp('Calculate mean firing rate for each session');
    mFRStructSess = cell(length(mazeSess),1);
    if(length(mazeSess)>1)
        for i = 1:length(mazeSess) 
            fprintf('\nSession %d\n',i);
            indLaps = find(beh.mazeSess == mazeSess(i));
            indLaps = intersect(indLaps,beh.indGoodLap); 
            %%% calculate mean firing rate for each neuron over specified trials
            mFRStructSess{i} = MFR(spikes,indLaps,rec.numNeurons,...
                        beh.lenTrials/sampleFq,indRunInLap);
        end
    end
    
    %%%%%%%%%% separate excitatory neurons and inhibitory neurons    
%     rec.indInhNeurons = find(mFRStruct.mFR > maxFR ...
%         & cluList.isIntern == 1);
%     rec.indExcNeurons = find(mFRStruct.mFR > minFR ...
%         & cluList.refracViolPercent < refracViolPercentThre...
%         & cluList.mahalDist > mahalDistThre ...
%         & cluList.centerMax < centerMaxThre ...
%         & cluList.isIntern == 0);
    
    save([path fileNameFR], 'mFRStruct','mFRStructSess');
    
    %%%%%%%%% draw figure is the state is on
    if(figureState ~= 0)
        % Ensure root units are pixels and get the size of the screen and create a
        % figure window
        set(0,'Units','pixels') 
        
        % plot mean firing rate
        plotMFR(rec.numNeurons,mFRStruct.mFR,mFRStruct.stdMFR)
        title('mean FR, all neurons')
        
        % plot the histogram of mean firing rate
        plotMFRHist(mFRStruct.mFR);
        title('mean FR, all neurons')
        
        % plot mean firing rate for exc neurons
        plotMFRHist(mFRStruct.mFR(autoCorr.isPyrNeuron));
        if(onlyRun == 0)
            title('mean FR, exc neurons')
        else
            % plot mean firing rate for running period and exc neurons
            title('mean FR during running, exc neurons')
        end
    end
    
    clear trials spikes mFRStruct mFRStructSess rec beh indRunInLap indLaps
