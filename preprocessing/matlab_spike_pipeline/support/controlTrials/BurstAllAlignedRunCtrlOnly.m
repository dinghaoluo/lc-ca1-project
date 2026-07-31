function BurstAllAlignedRunCtrlOnly(path,fileName,methodTheta,figureState,onlyRun,mazeSess)
% Find all the bursts which satisfy the criteria (the isi < 6 ms)
% Dependence:  function "ThetaPhaseLR" should be executed first
% path:         the path of the recording file
% fileName:     name of the recording file
% methodTheta:  method used for theta phase estimation
%               0: Hilbert transform
%               1: Linear interpolation
% figureState:  0: figure off
%               1: plot the histogram of the burst time
%               2: plot the burst time and the burst phase during the analysis
% onlyRun:      1: only consider the time period when the animal is running 
%
% Example:
% BurstAllVR('./','xzvr_PR1-20170727-01_DataStructure_mazeSection1_TrialType1',1,0,1)

    %%%%%%%%% check arguments
    if nargin<2
        disp('At least two arguments are needed for this function.');
        return;
    elseif nargin == 2
        methodTheta = 1;
        figureState = 1;
        onlyRun = 1;
        mazeSess = 1;
    elseif nargin == 3
        figureState = 1;
        onlyRun = 1;
        mazeSess = 1;
    elseif nargin == 4
        onlyRun = 1;
        mazeSess = 1;
    elseif nargin == 5
        mazeSess = 1;
    elseif nargin > 6
        disp('Too many input arguments.');
        return;
    end
    
    %%%%%%%%% initialize constants
    GlobalConst;
    
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    
    if(methodTheta == 0)
        th = 'H';
    else
        th = 'L';
    end
    fileNameBurst = [fileName '_burstAllAlignedRunCtrl_TH' th '_msess' num2str(mazeSess) '_Run' num2str(onlyRun) ...
                     '.mat'];

    if(methodTheta == 0)
        fileNameThetaPhase = [fileName '_ThetaPhaseHAlignedCtrl_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    else
        fileNameThetaPhase = [fileName '_ThetaPhaseLAlignedCtrl_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    end
    fileNameOrig = [fileName '.mat'];
    
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
    
    fullPath = [path fileName '_Info.mat']; 
    if(exist(fullPath) == 0)
        disp('The behavioral information file does not exist');
        return;
    end
    load(fullPath,'rec');
   
    % burst parameters
    paramBurst = struct('numNeurons', rec.numNeurons,...
                        'burstIsi', burstIsi,...
                        'burstIsi1st', burstIsi1st);
    
    % extract fields from trials structure
    
    spikes = trialsRunSpikes.Time;
    dist = trialsRunSpikes.Dist;
    if(methodTheta == 0)
        theta = trialsRunSpikes.thPhaseHilbSpike;
    else
        theta = trialsRunSpikes.thPhaseInterpSpike;
    end
   
    %%
    disp('Calculate bursts for all the non-stimulated laps');
    %%%%%%%%%%% calculate burst per neuron per trial
    disp('Calculate bursts per neuron per trial')
    burstPerNeuPerTrNonStim = ...
        BurstPerNeuronPerTrialAligned(spikes,theta,dist,trialNoNonStim,paramBurst);     
                                    % struct including the start
                                    % time of each burst, each cell of a
                                    % field containing the bursts from one
                                    % particular neuron and one particular trial
        
    %%%%%%%%%%% calculate bursts per neuron
    disp('Calculate bursts per neuron')
    burstIsiPerNeuronNonStim = BurstPerNeuron(burstPerNeuPerTrNonStim,trialNoNonStim,...
                                       rec.numNeurons);  
                                  
    %%%%%%%%% save data
    fullPath = [path fileNameBurst];
    save(fullPath, 'burstPerNeuPerTrNonStim','burstIsiPerNeuronNonStim',...
         'paramBurst');

end
