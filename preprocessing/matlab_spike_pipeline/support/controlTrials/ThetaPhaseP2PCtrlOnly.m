 function ThetaPhaseP2PCtrlOnly(path,fileName,method,onlyRun,mazeSess)
% Plot the theta phase and the smoothed mean firing rate curves for the
% selected neurons (a theta cycle is detected as peak to peak distance)
% path:         path of the recording file
% fileName:     name of the recording file
% method:       0: hilbert transform
%               1: linear interpolation
% onlyRun:      1: only consider the time period when the animal is running 
%
% Example:
% ThetaPhaseP2PCtrlOnly('./','A111-20150301-01_DataStructure_mazeSection1_TrialType1',1,1,1)

    %%%%%%%%% check arguments
    if nargin<2
        disp('At least three arguments are needed for this function.');
        return;
    elseif nargin == 2
        onlyRun = 1;
        method = 0;
        mazeSess = 1;
    elseif nargin == 3
        onlyRun = 1;
        mazeSess = 1;
    elseif nargin == 4
        mazeSess = 1;
    elseif nargin > 5
        disp('Too many input arguments');        
        return;
    end
    figureState = 0;
        
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    fileNameFR = [fileName '_FR_Ctrl_Run' num2str(onlyRun) '_mazeSess' num2str(mazeSess) '.mat'];
    if(method == 0)
        fileNameThetaPhase = [fileName '_ThetaPhaseH_Ctrl_Run' num2str(onlyRun) '_mazeSess' num2str(mazeSess) '.mat'];
    else
        fileNameThetaPhase = [fileName '_ThetaPhaseL_Ctrl_Run' num2str(onlyRun) '_mazeSess' num2str(mazeSess) '.mat'];
    end
    fileNameExt = [fileName '_ext.mat'];
    fileName = [fileName '.mat'];
    
    %% changed by Yingxue on 1/22/2022 to load file based on onlyRun
    fullPath = [path fileName]; 
    if(exist(fullPath) == 0)
        disp('The file does not exist');
        return;
    end
    load(fullPath);
    
    fullPath = [path fileNameExt];
    if(exist(fullPath) == 0)
        disp(['Extended file does not exist. Please run function',...
              ' SpikeDuringRun first']);
        return;
    end
    load(fullPath);
    
    fullPath = [path fileNameFR];
    if(exist(fullPath,'file') == 0)
        disp('_FR_Ctrl file does not exist. Please run "MeanFiringRateCtrlOnly" function first');
    end
    load(fullPath,'mFRStructSessCtrl');
                       
    %%%%%%%%% initialize constants
    GlobalConst;
    fileNameInfo = [fileName(1:end-4) '_Info.mat'];
        
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo_smTr(path,fileName);
    end
    load(fullPath);
     
    %%%%%%%%%% extract trials
    
    %%%%%%%%% extract theta phase from the data structure
    if(method == 0)
        theta = getRecField(trials,'thetaHil',1:length(lapList));
    else
        theta = getRecField(trials,'thetaLin',1:length(lapList));
    end
    if(onlyRun == 0)
        tmp = trials;
    else
        tmp = trialsExt;
    end
    spikes = getRecField(tmp,'spikes',1:length(lapList));    
    spikesDist = getRecField(tmp,'spikesMM',1:length(lapList));    
    if(method == 0)
        spikeTheta = getRecField(tmp,'spikesThetaHil',1:length(lapList));        
    else
        spikeTheta = getRecField(tmp,'spikesThetaLin',1:length(lapList));
    end
                
    disp('Calculate spike theta phase for control trials')
    indLaps = mFRStructSessCtrl.indLapList;
    thetaCycleStructTmp = ThetaCycle1(theta, indLaps, ...
                                minSamplePerCycle, thetaPhaseJump);
    spikeThetaPhaseStructSessCtrl = ...
        SpikeThetaPhase(spikes, spikesDist, spikeTheta,...
        thetaCycleStructTmp, indLaps, rec.numNeurons, ...
        timeStep, figureState);
           
    fullPath = [path fileNameThetaPhase];
    save(fullPath, 'spikeThetaPhaseStructSessCtrl');
        
    clear mydata;
    clear all;
