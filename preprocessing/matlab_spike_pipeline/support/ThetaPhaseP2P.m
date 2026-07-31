 function ThetaPhaseP2P(path,fileName,method,onlyRun,figureState)
% Plot the theta phase and the smoothed mean firing rate curves for the
% selected neurons (a theta cycle is detected as peak to peak distance)
% path:         path of the recording file
% fileName:     name of the recording file
% method:       0: hilbert transform
%               1: linear interpolation
% onlyRun:      1: only consider the time period when the animal is running 
% figureState:  0: figure off
%               1: figure on
%               2: plot the phase histogram before the analysis of the starting phase of each neuron 
%               (Since the phase is periodic, the phases are actually distributed on the surface of the unit cylinder. 
%               By first estimating the phase where the cycle starts and ends for individual neuron, 
%               we can then use linear regression to obtain the phase
%               precession slope)
%
% Example:
% ThetaPhaseP2P('./','A111-20150301-01_DataStructure_mazeSection1_TrialType1',1,1,1)

    %%%%%%%%% check arguments
    if nargin<2
        disp('At least three arguments are needed for this function.');
        return;
    elseif nargin == 2
        figureState = 0;
        onlyRun = 1;
        method = 0;
    elseif nargin == 3
        onlyRun = 1;
        figureState = 0;
    elseif nargin == 4
        figureState = 0;
    elseif nargin > 5
        disp('Too many input arguments');        
        return;
    end
        
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    if(method == 0)
        fileNameThetaPhase = [fileName '_ThetaPhaseH_Run' num2str(onlyRun) '.mat'];
    else
        fileNameThetaPhase = [fileName '_ThetaPhaseL_Run' num2str(onlyRun) '.mat'];
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
                
    %%%%%%%%% collect information about each theta cycle during each trial
    disp('Calculate theta cycle for all the trials')
    thetaCycleStruct = ThetaCycle1(theta, beh.indGoodLap, ...
                                    minSamplePerCycle, thetaPhaseJump);
            
    %%%%%%%%% collect information about the spike theta phases
    %%% calculate theta phase with the numSamples = the shortest trial
    disp('Calculate spike theta phase for all the trials')
    spikeThetaPhaseStruct = SpikeThetaPhase(spikes, spikesDist, spikeTheta,...
        thetaCycleStruct, beh.indGoodLap, rec.numNeurons, timeStep, figureState);
        
    spikeThetaPhaseStructSess = cell(length(mazeSess),1);
    if(length(mazeSess)>1)
        disp('Calculate spike theta phase for each session')
        for i = 1:length(mazeSess) 
            fprintf('\nSession %d\n',i);
            indLaps = find(beh.mazeSess == mazeSess(i)); 
            indLaps = intersect(indLaps,beh.indGoodLap);
            %% added by Yingxue on 12/19/2020
            thetaCycleStructTmp = ThetaCycle1(theta, indLaps, ...
                                    minSamplePerCycle, thetaPhaseJump);
            spikeThetaPhaseStructSess{i} = ...
                SpikeThetaPhase(spikes, spikesDist, spikeTheta,...
                thetaCycleStructTmp, indLaps, rec.numNeurons, ...
                timeStep, figureState);         
        end
    end
    
    fullPath = [path fileNameThetaPhase];
    save(fullPath, 'thetaCycleStruct', 'spikeThetaPhaseStruct',...
            'spikeThetaPhaseStructSess');
                   
    clear mydata;
    clear all;
    
    %%%%%%%%% draw figure (theta phase together with the smoothed mean
    %%%%%%%%% firing rate curve)
%     plotThetaCycle(thetaCycleStruct,timeStep); 
%     plotSpikeThetaPhase(spikeThetaPhaseStruct,timeStep,numSamples);
