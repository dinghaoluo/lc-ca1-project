 function ThetaPhaseP2PAlignRun(path,fileName,method,onlyRun,mazeSess)
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
        fileNameThetaPhase = [fileName '_ThetaPhaseHAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    else
        fileNameThetaPhase = [fileName '_ThetaPhaseLAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    end
    
    %%%%%%%%% initialize constants
    GlobalConst;
    fileNameInfo = [fileName '_Info.mat'];
                
    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
                
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
    
    trialNo = 1:length(behPar.indTrBadBeh);         
    thetaPh = cell(1,length(trialNo));
    
    for i = trialNo
        if(method == 0) % hilbert
            thetaPh{i} = trialsRun.thetaPhHilb{i}';
            
        else % linear
            thetaPh{i} = trialsRun.thetaPhLinInterp{i}';
        end
    end
    
    if(method == 0)
        thPhase = trialsRunSpikes.thPhaseHilbSpike;
    else
        thPhase = trialsRunSpikes.thPhaseInterpSpike;
    end
            
    %%%%%%%%% collect information about each theta cycle during each trial
    disp('Calculate theta cycle for all the trials')
    thetaCycleRun = ThetaCycle1(thetaPh, trialNo, ...
                                    minSamplePerCycle, thetaPhaseJump);
            
    %%%%%%%%% collect information about the spike theta phases
    %%% calculate theta phase with the numSamples = the shortest trial
    disp('Calculate spike theta phase for non-stimulated good trials')
    noNeurons = size(trialsRunSpikes.Time,1);
    spikeThetaPhaseRunNoStimGood = SpikeThetaPhaseAlignedRun(trialsRunSpikes.Time,...
        thPhase, thetaCycleRun, trialNoNonStimGood, noNeurons);
    
    disp('Calculate spike theta phase for non-stimulated bad trials')
    spikeThetaPhaseRunNoStimBad = SpikeThetaPhaseAlignedRun(trialsRunSpikes.Time,...
        thPhase, thetaCycleRun, trialNoNonStimBad, noNeurons);
    
    %% added by Yingxue on 7/15/2021
    spikeThetaPhaseRunStim = [];
    spikeThetaPhaseRunStimCtrl = [];
    for i = 1:length(pulseMeth)
        disp('Calculate spike theta phase for stimulated trials')
        spikeThetaPhaseRunStim{i} = SpikeThetaPhaseAlignedRun(trialsRunSpikes.Time,...
            thPhase, thetaCycleRun, trialNoStim{i}, noNeurons);

        %% added by Yingxue on 7/14/2021
        disp('Calculate spike theta phase for stimulated control trials')
        spikeThetaPhaseRunStimCtrl{i} = SpikeThetaPhaseAlignedRun(trialsRunSpikes.Time,...
            thPhase, thetaCycleRun, trialNoStimCtrl{i}, noNeurons);
        %%
    end
    
    fullPath = [path fileNameThetaPhase];
    save(fullPath, 'thetaCycleRun', 'spikeThetaPhaseRunNoStimGood',...
            'spikeThetaPhaseRunNoStimBad','spikeThetaPhaseRunStim',...
            'spikeThetaPhaseRunStimCtrl');
                   
    clear mydata;
    clear all;
    
    %%%%%%%%% draw figure (theta phase together with the smoothed mean
    %%%%%%%%% firing rate curve)
%     plotThetaCycle(thetaCycleStruct,timeStep); 
%     plotSpikeThetaPhase(spikeThetaPhaseStruct,timeStep,numSamples);
