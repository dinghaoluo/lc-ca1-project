 function ThetaPhaseP2PCueOffset(path,fileName,method,onlyRun,mazeSess)
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
        fileNameThetaPhase = [fileName '_ThetaPhaseH_msess' num2str(mazeSess) '_CueOffRun' num2str(onlyRun) '.mat'];
    else
        fileNameThetaPhase = [fileName '_ThetaPhaseL_msess' num2str(mazeSess) '_CueOffRun' num2str(onlyRun) '.mat'];
    end
    
    %%%%%%%%% initialize constants
    GlobalConst;
    fileNameInfo = [fileName '_Info.mat'];
        
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo(path,fileName);
    end
    load(fullPath);
        
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
    load(fullPath,'trialsCueOffSpikes');
    
    fullPath = [path fileName '_alignCueOff_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to run file does not exist');
        return;
    end
    load(fullPath,'trialsCueOff');
       
    trialNo = find(trialsCueOff.startLfpInd ~= -1);
    if(isfield(beh,'indStimLap'))
        ind = beh.mazeSess == mazeSess;
        trialNoNonStimGood = find(beh.indStimLap(ind) == 0 & behPar.indTrBadBeh == 0);
        trialNoNonStimBad = find(beh.indStimLap(ind) == 0 & behPar.indTrBadBeh == 1);
        trialNoStim = find(beh.indStimLap(ind) == 1);
    else
        trialNoNonStimGood = find(behPar.indTrBadBeh == 0);
        trialNoNonStimBad = find(behPar.indTrBadBeh == 1);
        trialNoStim = [];
    end
     
    thetaPh = cell(1,length(trialNo));
    
    for i = trialNo
        if(method == 0) % hilbert
            thetaPh{i} = [trialsCueOff.thetaPhHilbBef{i}',trialsCueOff.thetaPhHilb{i}'];
            
        else % linear
            thetaPh{i} = [trialsCueOff.thetaPhLinInterpBef{i}',trialsCueOff.thetaPhLinInterp{i}'];
        end
    end
            
    %%%%%%%%% collect information about each theta cycle during each trial
    disp('Calculate theta cycle for all the trials')
    thetaCycleCueOff = ThetaCycle1(thetaPh, trialNo, ...
                                    minSamplePerCycle, thetaPhaseJump);
            
    %%%%%%%%% collect information about the spike theta phases
    %%% calculate theta phase with the numSamples = the shortest trial
    disp('Calculate spike theta phase for non-stimulated good trials')
    noNeurons = size(trialsCueOffSpikes.Time,1);
    spikeThetaPhaseCueOffNoStimGood = SpikeThetaPhaseRunOnset(trialsCueOffSpikes.TimeBef,...
        trialsCueOffSpikes.Time, trialsCueOffSpikes.thPhaseInterpSpikeBef,...
        trialsCueOffSpikes.thPhaseInterpSpike, thetaCycleCueOff, trialNoNonStimGood, noNeurons);
    
    disp('Calculate spike theta phase for non-stimulated bad trials')
    spikeThetaPhaseCueOffNoStimBad = SpikeThetaPhaseRunOnset(trialsCueOffSpikes.TimeBef,...
        trialsCueOffSpikes.Time, trialsCueOffSpikes.thPhaseInterpSpikeBef,...
        trialsCueOffSpikes.thPhaseInterpSpike, thetaCycleCueOff, trialNoNonStimBad, noNeurons);
    
    disp('Calculate spike theta phase for stimulated trials')
    spikeThetaPhaseCueOffStim = SpikeThetaPhaseRunOnset(trialsCueOffSpikes.TimeBef,...
        trialsCueOffSpikes.Time, trialsCueOffSpikes.thPhaseInterpSpikeBef,...
        trialsCueOffSpikes.thPhaseInterpSpike, thetaCycleCueOff, trialNoStim, noNeurons);
    
    fullPath = [path fileNameThetaPhase];
    save(fullPath, 'thetaCycleCueOff', 'spikeThetaPhaseCueOffNoStimGood',...
            'spikeThetaPhaseCueOffNoStimBad','spikeThetaPhaseCueOffStim');
                   
    clear mydata;
    clear all;
    
    %%%%%%%%% draw figure (theta phase together with the smoothed mean
    %%%%%%%%% firing rate curve)
%     plotThetaCycle(thetaCycleStruct,timeStep); 
%     plotSpikeThetaPhase(spikeThetaPhaseStruct,timeStep,numSamples);
