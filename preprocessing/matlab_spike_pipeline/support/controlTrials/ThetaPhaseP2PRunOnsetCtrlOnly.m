 function ThetaPhaseP2PRunOnsetCtrlOnly(path,fileName,method,onlyRun,mazeSess)
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
        fileNameThetaPhaseCtrl = [fileName '_ThetaPhaseHCtrl_msess' num2str(mazeSess) '_RunOnset' num2str(onlyRun) '.mat'];
        fileNameThetaPhase = [fileName '_ThetaPhaseH_msess' num2str(mazeSess) '_RunOnset' num2str(onlyRun) '.mat'];
    else
        fileNameThetaPhaseCtrl = [fileName '_ThetaPhaseLCtrl_msess' num2str(mazeSess) '_RunOnset' num2str(onlyRun) '.mat'];
        fileNameThetaPhase = [fileName '_ThetaPhaseL_msess' num2str(mazeSess) '_RunOnset' num2str(onlyRun) '.mat'];
    end
    fullPath = [path fileNameThetaPhase];
    if(exist(fullPath) == 0)
        disp('The _ThetaPhase_RunOnset file does not exist');
        return;
    end
    load(fullPath,'thetaCycleRun');
    
    %%%%%%%%% initialize constants
    GlobalConst;
    
    fullPath = [path fileName '_alignedSpikesPerNPerT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _alignedSpikesPerNPerT file does not exist');
        return;
    end
    load(fullPath,'trialsRunSpikes');
       
    fullPath = [path fileName '_PeakFRAlignedCtrl_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];    
    if(exist(fullPath) == 0)
        disp('The _PeakFRAlignedCtrl file does not exist');
        return;
    end
    load(fullPath,'trialNoNonStim');
            
    
    %%%%%%%%% collect information about the spike theta phases
    %%% calculate theta phase with the numSamples = the shortest trial
    disp('Calculate spike theta phase for non-stimulated trials')
    noNeurons = size(trialsRunSpikes.Time,1);
    spikeThetaPhaseRunNoStim = SpikeThetaPhaseRunOnset(trialsRunSpikes.TimeBef,...
        trialsRunSpikes.Time, trialsRunSpikes.thPhaseInterpSpikeBef,...
        trialsRunSpikes.thPhaseInterpSpike, thetaCycleRun, trialNoNonStim, noNeurons);
    
    fullPath = [path fileNameThetaPhaseCtrl];
    save(fullPath, 'spikeThetaPhaseRunNoStim');
                   
    clear mydata;
    clear all;
