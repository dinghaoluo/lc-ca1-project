 function ThetaPhaseP2PAlignRunSegBefRun(path,fileName,method,onlyRun,mazeSess)
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
% ThetaPhaseP2PRunOnset('./','A111-20150301-01_DataStructure_mazeSection1_TrialType1',1,1,1)

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
        
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    if(method == 0)
        fileNameThetaPhase = [fileName '_ThetaPhaseH_msess' num2str(mazeSess) '_RunOnset' num2str(onlyRun) '.mat'];
        fileNameThetaPhaseSeg = [fileName '_ThetaPhaseHSeg_msess' num2str(mazeSess) '_RunOnset' num2str(onlyRun) '.mat'];
    else
        fileNameThetaPhase = [fileName '_ThetaPhaseL_msess' num2str(mazeSess) '_RunOnset' num2str(onlyRun) '.mat'];
        fileNameThetaPhaseSeg = [fileName '_ThetaPhaseLSeg_msess' num2str(mazeSess) '_RunOnset' num2str(onlyRun) '.mat'];
    end
    
    %%%%%%%%% initialize constants
    GlobalConst;
    fullPath = [path fileNameThetaPhase]; 
    if(exist(fullPath) == 0)
        disp('The _ThetaPhaseL_RunOnset file does not exist');
        return;
    end
    load(fullPath,'spikeThetaPhaseRunNoStimGood');
    
    disp('Calculate spike theta phase for non-stimulated good trials between 0.5s before run')
    spikeThetaPhaseRunNoStimGoodbefRun = SpikeThetaPhaseAlignedRunSeg(spikeThetaPhaseRunNoStimGood,...
        [-0.5 0],sampleFq);
    
    fullPath = [path fileNameThetaPhaseSeg];
    save(fullPath, 'spikeThetaPhaseRunNoStimGoodbefRun');
                   
    clear mydata;
    clear all;
    
    %%%%%%%%% draw figure (theta phase together with the smoothed mean
    %%%%%%%%% firing rate curve)
%     plotThetaCycle(thetaCycleStruct,timeStep); 
%     plotSpikeThetaPhase(spikeThetaPhaseStruct,timeStep,numSamples);
