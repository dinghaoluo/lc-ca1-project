 function ThetaPhaseP2P_GoodTr(path,fileName,method,onlyRun)
% Plot the theta phase and the smoothed mean firing rate curves for the
% selected neurons (a theta cycle is detected as peak to peak distance)
% path:         path of the recording file
% fileName:     name of the recording file
% method:       0: hilbert transform
%               1: linear interpolation
% onlyRun:      1: only consider the time period when the animal is running 
%
% Example:
% ThetaPhaseP2P_GoodTr('./','A111-20150301-01_DataStructure_mazeSection1_TrialType1',1,1,1)

    %%%%%%%%% check arguments
    if nargin<2
        disp('At least three arguments are needed for this function.');
        return;
    elseif nargin == 2
        onlyRun = 1;
        method = 0;
    elseif nargin == 3
        onlyRun = 1;
    elseif nargin > 4
        disp('Too many input arguments');        
        return;
    end
    figureState = 0;
        
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    if(method == 0)
        fileNameThetaPhase = [fileName '_ThetaPhaseH_GoodTr_Run' num2str(onlyRun) '.mat'];
    else
        fileNameThetaPhase = [fileName '_ThetaPhaseL_GoodTr_Run' num2str(onlyRun) '.mat'];
    end
    fileNameExt = [fileName '_ext.mat'];
    fileName = [fileName '.mat'];
    
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
                       
    %%%%%%%%% initialize constants
    GlobalConst;
    fileNameInfo = [fileName(1:end-4) '_Info.mat'];
        
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo_smTr(path,fileName);
    end
    load(fullPath);
    mazeSess = beh.mazeSessAll;
       
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
                
    spikeThetaPhaseStructSessGoodTr = cell(length(mazeSess),1);
    spikeThetaPhaseStructSessBadTr = cell(length(mazeSess),1);
    spikeThetaPhaseStructSessOKTr = cell(length(mazeSess),1);
    if(length(mazeSess)>1)
        disp('Calculate spike theta phase for each session')
        for i = 1:length(mazeSess) 
            fprintf('\nSession %d\n',i);
            indLaps = find(beh.mazeSess == mazeSess(i));
            indLaps = intersect(indLaps,beh.indTrCtrl); 
            %%% calculate peak firing rate for good, bad and ok trials
            %%% within each session
            %%% added by Yingxue on 12/16/2020
            indLapsGoodTr = intersect(indLaps,beh.indGoodTrCtrl);
            thetaCycleStructTmp = ThetaCycle1(theta, indLapsGoodTr, ...
                                    minSamplePerCycle, thetaPhaseJump);
            spikeThetaPhaseStructSessGoodTr{i} = ...
                SpikeThetaPhase(spikes, spikesDist, spikeTheta,...
                thetaCycleStructTmp, indLapsGoodTr, rec.numNeurons, ...
                timeStep, figureState);
            indLapsBadTr = intersect(indLaps,beh.indBadTrCtrl); 
            thetaCycleStructTmp = ThetaCycle1(theta, indLapsBadTr, ...
                                    minSamplePerCycle, thetaPhaseJump);
            spikeThetaPhaseStructSessBadTr{i} = ...
                SpikeThetaPhase(spikes, spikesDist, spikeTheta,...
                thetaCycleStructTmp, indLapsBadTr, rec.numNeurons, ...
                timeStep, figureState);
            indLapsOKTr = setdiff(indLaps,[beh.indGoodTrCtrl,beh.indBadTrCtrl]); 
            thetaCycleStructTmp = ThetaCycle1(theta, indLapsOKTr, ...
                                    minSamplePerCycle, thetaPhaseJump);
            spikeThetaPhaseStructSessOKTr{i} = ...
                SpikeThetaPhase(spikes, spikesDist, spikeTheta,...
                thetaCycleStructTmp, indLapsOKTr, rec.numNeurons, ...
                timeStep, figureState);            
        end
    else
        indLapsGoodTr = beh.indGoodTrCtrl; 
        thetaCycleStructTmp = ThetaCycle1(theta, indLapsGoodTr, ...
                                    minSamplePerCycle, thetaPhaseJump);
        spikeThetaPhaseStructSessGoodTr{1} = ...
            SpikeThetaPhase(spikes, spikesDist, spikeTheta,...
            thetaCycleStructTmp, indLapsGoodTr, rec.numNeurons, ...
            timeStep, figureState);
        indLapsBadTr = beh.indBadTrCtrl; 
        thetaCycleStructTmp = ThetaCycle1(theta, indLapsBadTr, ...
                                    minSamplePerCycle, thetaPhaseJump);
        spikeThetaPhaseStructSessBadTr{1} = ...
            SpikeThetaPhase(spikes, spikesDist, spikeTheta,...
            thetaCycleStructTmp, indLapsBadTr, rec.numNeurons, ...
            timeStep, figureState);
        indLapsOKTr = setdiff(beh.indTrCtrl,[beh.indGoodTrCtrl,beh.indBadTrCtrl]);
        thetaCycleStructTmp = ThetaCycle1(theta, indLapsOKTr, ...
                                    minSamplePerCycle, thetaPhaseJump);
        spikeThetaPhaseStructSessOKTr{1} = ...
            SpikeThetaPhase(spikes, spikesDist, spikeTheta,...
            thetaCycleStructTmp, indLapsOKTr, rec.numNeurons, ...
            timeStep, figureState);  
    end
    
    fullPath = [path fileNameThetaPhase];
    save(fullPath, 'spikeThetaPhaseStructSessGoodTr',...
            'spikeThetaPhaseStructSessBadTr',...
            'spikeThetaPhaseStructSessOKTr');
                   
    clear mydata;
    clear all;
