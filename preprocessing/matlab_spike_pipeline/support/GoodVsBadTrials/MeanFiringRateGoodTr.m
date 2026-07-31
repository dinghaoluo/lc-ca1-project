function MeanFiringRateGoodTr(path,fileName,onlyRun)
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
    elseif nargin > 3
        disp('Too many input arguments.');
        return;
    end
        
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    fileNameFR = [fileName '_FR_GoodTr_Run' num2str(onlyRun) '.mat'];
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
    mFRStructSessGoodTr = cell(length(mazeSess),1);
    mFRStructSessBadTr = cell(length(mazeSess),1);
    mFRStructSessOKTr = cell(length(mazeSess),1);
    if(length(mazeSess)>1)
        for i = 1:length(mazeSess) 
            fprintf('\nSession %d\n',i);
            %%% calculate mean firing rate for good, bad and ok trials
            %%% added by Yingxue on 12/16/2020
            indLaps = find(beh.mazeSess == mazeSess(i));
            indLaps = intersect(indLaps,beh.indTrCtrl); 
            indLapsGoodTr = intersect(indLaps,beh.indGoodTrCtrl); 
            mFRStructSessGoodTr{i} = MFR(spikes,indLapsGoodTr,rec.numNeurons,...
                        beh.lenTrials/sampleFq,indRunInLap);
            indLapsBadTr = intersect(indLaps,beh.indBadTrCtrl); 
            mFRStructSessBadTr{i} = MFR(spikes,indLapsBadTr,rec.numNeurons,...
                        beh.lenTrials/sampleFq,indRunInLap);
            indLapsOKTr = setdiff(indLaps,[beh.indGoodTrCtrl,beh.indBadTrCtrl]); 
            mFRStructSessOKTr{i} = MFR(spikes,indLapsOKTr,rec.numNeurons,...
                        beh.lenTrials/sampleFq,indRunInLap);
        end
    else
        indLapsGoodTr = beh.indGoodTrCtrl; 
        mFRStructSessGoodTr{1} = MFR(spikes,indLapsGoodTr,rec.numNeurons,...
                    beh.lenTrials/sampleFq,indRunInLap);
        indLapsBadTr = beh.indBadTrCtrl; 
        mFRStructSessBadTr{1} = MFR(spikes,indLapsBadTr,rec.numNeurons,...
                    beh.lenTrials/sampleFq,indRunInLap);
        indLapsOKTr = setdiff(beh.indTrCtrl,[beh.indGoodTrCtrl,beh.indBadTrCtrl]);
        mFRStructSessOKTr{1} = MFR(spikes,indLapsOKTr,rec.numNeurons,...
                    beh.lenTrials/sampleFq,indRunInLap);
    end
    
    %%%%%%%%%% separate excitatory neurons and inhibitory neurons    
%     rec.indInhNeurons = find(mFRStruct.mFR > maxFR ...
%         & cluList.isIntern == 1);
%     rec.indExcNeurons = find(mFRStruct.mFR > minFR ...
%         & cluList.refracViolPercent < refracViolPercentThre...
%         & cluList.mahalDist > mahalDistThre ...
%         & cluList.centerMax < centerMaxThre ...
%         & cluList.isIntern == 0);
    
    save([path fileNameFR], 'mFRStructSessGoodTr',...
        'mFRStructSessBadTr','mFRStructSessOKTr');
       
    clear trials spikes rec beh indRunInLap indLapsGoodTr indLapsBadTr indLapsOKTr...
        mFRStructSessGoodTr mFRStructSessBadTr mFRStructSessOKTr
