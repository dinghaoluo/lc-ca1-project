function [Track,Laps] = SortTrials_3armMaze_smTr(FileNameBase, lfpFreq, sampleFreq, isInt)

% sorts maze Tracks according to different bevaioral parameters and
% calculate various behavioral parameters
%
% output: matlab file '*_Tracks.mat'
%

% cd /groups/pastalkova/home/pastalkovae/data/data/Yingxue_recordings/A632/A632_110725_01_e
% SortTrials_3armMaze_v1('A632_110725_01');

    % added 08/14/2018, removed _Track_Laps.mat file to reduce data duplication
    if exist([FileNameBase '_Behav2PDataLFP.mat'], 'file') == 2
        fprintf(...
            '\nCheck whether %s file already contains Track and Laps.\n', ...
            [FileNameBase '_Behav2PDataLFP.mat']);
        load([FileNameBase '_Behav2PDataLFP.mat'], 'Processing');
        if(sum(Processing == 1) > 0)
            return;
        end
    end
    
    %% Load tracking file
    load([FileNameBase '-whl.mat']) ;            
    whlm = length(whlDataLfp(:,2));
   
     %% load Arduino behavioral data with TDT time stamps
    load([FileNameBase 'BTDT.mat']);                       
    noLfpTrials = behEventsTdt.trialDescr(:,4)>whlm;
    behEventsTdt.trialDescr(noLfpTrials,:) = [];

    totNLaps = behEventsTdt.trialDescr(end,5);           % total N of laps
    
    %% Initialize output structure
    Track.xMM = whlDataLfp(:,2);
    Track.xMMAll = whlDataLfp(:,3); % including the stopping period
    Track.startLfpInd = zeros(whlm,1);
    Track.lapID = zeros(whlm,1);                % order of Tracks
    Track.corrChoice = zeros(whlm,1);       % error[0]/correct[1]/not assigned[-1]
    Track.behavType = zeros(whlm,1);        % good-1, late decision-2, exploration-3, bad=not running-4
    Track.speed_MMsec = zeros(whlm,1);
    Track.accel_MMsecSq = zeros(whlm,1);
    Track.speed_MMsecAll = zeros(whlm,1); % including the stopping period
    Track.accel_MMsecSqAll = zeros(whlm,1); % including the stopping period
    Track.mazeType = zeros(whlm,1);
    Track.mazeSess = zeros(whlm,1);
    
    Track.dFF = interpTspdFF';
    Track.F = interpTspF';
    Track.Fneu = interpTspFneu';
    Track.spks = interpTspSpks';

    Laps.behavType = behType;
    Laps.mazeType = behEventsTdt.mazeType';
    Laps.mazeSess = behEventsTdt.mazeSess';
    Laps.trackLen = trackLenArr;
    Laps.lapID = zeros(totNLaps,1);
    Laps.startLfpInd = zeros(totNLaps,1);
    Laps.endLfpInd = zeros(totNLaps,1);
    Laps.startT = zeros(totNLaps,1);
    Laps.endT = zeros(totNLaps,1);
    Laps.corrChoice = zeros(totNLaps,1);
    Laps.lickT = cell(totNLaps,1);
    Laps.lickLfpInd = cell(totNLaps,1);
        % lick events per trial
    Laps.pumpT = cell(totNLaps,1);
    Laps.pumpLfpInd = cell(totNLaps,1);
        % pump events per trial
    Laps.movieOnT = cell(totNLaps,1);
    Laps.movieOnLfpInd = cell(totNLaps,1);
    Laps.movieLocation = zeros(totNLaps,1);
        % movie on time per trial
    Laps.lickPeriodT = cell(totNLaps,1);
    Laps.lickPeriodInd = cell(totNLaps,1);
        % lick period per trial
    Laps.airpuffT = cell(totNLaps,1);
    Laps.airpuffLfpInd = cell(totNLaps,1);
        % airpuff events per trial
        
    %%%%%% added by Yingxue 11/2/2019    
    Laps.stimOnT = cell(totNLaps,1);
    Laps.stimOnLfpInd = cell(totNLaps,1);
        % stimulation per trial
    Laps.stimPulseT = cell(totNLaps,1);
    Laps.stimPulseLfpInd = cell(totNLaps,1);
    Laps.stimPulseWidth = cell(totNLaps,1);
    Laps.stimPulseWidthLfp = cell(totNLaps,1);
    Laps.stimDiode = cell(totNLaps,1);
    Laps.stimDiodeCurr = cell(totNLaps,1);
    Laps.stimPulseMethod = cell(totNLaps,1);
    Laps.stimPulseLoc = cell(totNLaps,1);
        % stimulation pulses
        
    Laps.movieOnPulseT = cell(totNLaps,1);
    Laps.movieOnPulseLfpInd = cell(totNLaps,1);
    Laps.movieOffPulseT = cell(totNLaps,1);
    Laps.movieOffPulseLfpInd = cell(totNLaps,1);
    
    load([FileNameBase '_corrFluo.mat'],'Clu');
    
    %% threshold spikes and dF/F based on dF/F noise level
    param.stdsm = 0.3*lfpFreq/6;
    param.threStd = 3;
    param.prctile = 90;
    h = gaussFilter2P(6*param.stdsm,param.stdsm);
    lenGaussKernel = length(h);
    
    Track.dFFSM = zeros(size(Track.dFF,1),size(Track.dFF,2));
    Track.indNoise = zeros(size(Track.dFF,1),size(Track.dFF,2));
    Track.spikesSM = Track.spks;
            
    for i = 1:length(Clu.localClu)
        dFFSM = gauss_filter(Track.dFF(:,i),h,lenGaussKernel);
        Track.dFFGF(:,i) = dFFSM; % changed on 9/6/2023, saving smoothed dFF
        if(isInt == 0) % pyramidal neurons
            stdDFF = mad(dFFSM,1)/0.6745; % robust SD, David Tank, nature, 2021
            indNoise = dFFSM <= param.threStd * stdDFF;
    %         stdDFF = prctile(dFFSM,param.prctile);
    %         indNoise = dFFSM <= stdDFF;
            Track.indNoise(:,i) = indNoise;
            dFFSM(indNoise) = 0;
            Track.dFFSM(:,i) = dFFSM;
            Track.spikesSM(indNoise,i) = 0;
        else % interneurons
            Track.dFFSM(:,i) = Track.dFF(:,i); % changed on 9/6/2023, using raw dFF instead of smoothed dFF
        end
    end 
       
    %%%%%%%%%
        
    %% Analyze individual trials
    % use behEventsTdt.TrackDescr to get the beginning of each Track
    % behEventsTdt.TrackDescr: total time, session time, Track ID, ????, correct/incorrect
    % correct for repeated trail IDs
    for nTr = 1 : totNLaps
        if(nTr == 135)
            a = 1;
        end
        trStartLfpInd = trStartLfpIndArr(nTr);
        trEndLfpInd = trEndLfpIndArr(nTr);
        if(nTr ~= totNLaps)
            trNextStartLfpInd = trStartLfpIndArr(nTr+1);
        else
            if(length(behEventsTdt.taskDescr(:,4)) == totNLaps + 1)
                trNextStartLfpInd = behEventsTdt.taskDescr(end,4);
            else
                trNextStartLfpInd = trEndLfpIndArr(nTr);
            end
        end
%       figure(10); clf;

        % adjust the start index for each trial
        if(behType(nTr) == 1)
            % check whether the beginning of each trial is at ~0
            % distance
            if(nTr ~= 1 & behType(nTr-1) == 1)
                shift = 0;
                while(Track.xMM(trStartLfpInd) > 50)
                    trStartLfpInd = trStartLfpInd+1;
                    shift = shift+1;
                end
                if(shift ~= 0)
                    disp(['Trial no. ' num2str(nTr) ...
                        ': the start lfp index is shifted forward by ' ...
                        num2str(shift) ' points.']);
                end
            end

            % check whether the end of each trial is at the max
            % distance
            shift = 0;
            while(Track.xMM(trEndLfpInd) < trackLenArr(nTr) - 50)
                trEndLfpInd = trEndLfpInd-1;
                shift = shift+1;
            end
            if(shift ~= 0)
                disp(['Trial no. ' num2str(nTr)...
                    ': the stop lfp index is shifted backward by '...
                    num2str(shift) ' points.']);
            else
                while(Track.xMM(trEndLfpInd+1) > trackLenArr(nTr) - 1)
                    trEndLfpInd = trEndLfpInd+1;
                    shift = shift+1;
                end
                if(shift ~= 0)
                    disp(['Trial no. ' num2str(nTr) ...
                        ': the stop lfp index is shifted forward by '...
                        num2str(shift) ' points.']);
                end
            end
        end         

        Track.startLfpInd(trStartLfpInd:trEndLfpInd,1) = trStartLfpInd;
        Laps.startLfpInd(nTr,1) = trStartLfpInd;
        Laps.endLfpInd(nTr,1) = trEndLfpInd;
        Laps.startT(nTr,1) = trStartLfpInd/lfpFreq;
        Laps.endT(nTr,1) = trEndLfpInd/lfpFreq;
        
        % pump on time
        ind = behEventsTdt.pump(:,4) >= trStartLfpInd &...
            behEventsTdt.pump(:,4) <= trEndLfpInd;
        Laps.pumpT{nTr} = whlDataLfp(behEventsTdt.pump(ind,4),1);
        Laps.pumpLfpInd{nTr} = behEventsTdt.pump(ind,4);

        % lick on time
        ind = behEventsTdt.lick(:,4) >= trStartLfpInd &...
            behEventsTdt.lick(:,4) <= trEndLfpInd;
        Laps.lickT{nTr} = whlDataLfp(behEventsTdt.lick(ind,4),1);
        Laps.lickLfpInd{nTr} = behEventsTdt.lick(ind,4);
        
        % airpuff on time
        if(isfield(behEventsTdt,'airpuff'))
            ind = behEventsTdt.airpuff(:,4) >= trStartLfpInd &...
                behEventsTdt.airpuff(:,4) <= trEndLfpInd;
            Laps.airpuffT{nTr} = whlDataLfp(behEventsTdt.airpuff(ind,4),1);
            Laps.airpuffLfpInd{nTr} = behEventsTdt.airpuff(ind,4);
        end
        
        % stimulation on time
        % added on 11/2/2019 by Yingxue
        if(isfield(behEventsTdt,'stimOn')) % need to check this session once starting optogenetics 11/1/2021
            pulseMethod = behEventsTdt.movieTDescr{nTr}(12);
            Laps.stimPulseMethod{nTr} = pulseMethod;
            ind = behEventsTdt.stimOn(:,4) >= trStartLfpInd &...
                behEventsTdt.stimOn(:,4) <= trEndLfpInd;
            if(sum(ind) > 0)
                Laps.stimOnT{nTr} = whlDataLfp(behEventsTdt.stimOn(ind,4),1);
                Laps.stimOnLfpInd{nTr} = behEventsTdt.stimOn(ind,4);
                
                if(pulseMethod > 0)
                    indS = behEventsTdt.stimPulse(:,4) >= trStartLfpInd &...
                        behEventsTdt.stimPulse(:,4) <= trEndLfpInd;
                    if(sum(indS) ~= behEventsTdt.stimOn(ind,6))
                        indTmp = find(ind == 1);
                        disp(['Please check stimulation no. ' num2str(indTmp) ...
                            ', the no. of pulses does not match']);
                    end
                    Laps.stimPulseT{nTr} = whlDataLfp(behEventsTdt.stimPulse(indS,4),1);
                    Laps.stimPulseLfpInd{nTr} = behEventsTdt.stimPulse(indS,4);
                    Laps.stimPulseWidth{nTr} = behEventsTdt.stimPulse(indS,8);
                    Laps.stimPulseWidthLfp{nTr} = behEventsTdt.stimPulse(indS,9);
                    Laps.stimDiode{nTr} = [];
                    Laps.stimDiodeCurr{nTr} = [];
                    indS1 = find(indS == 1);
                    for m = 1:sum(indS)
                        Laps.stimDiode{nTr} = [Laps.stimDiode{nTr}; ...
                            behEventsTdt.stimOnDiode{indS1(m)}.indDiode];
                        Laps.stimDiodeCurr{nTr} = [Laps.stimDiodeCurr{nTr}; ...
                            behEventsTdt.stimOnDiode{indS1(m)}.diodeCurr]; 
                    end
                    
                    if(pulseMethod == 3)
                        Laps.stimPulseLoc{nTr} = behEventsTdt.movieTDescr{nTr}(9);
                    else
                        Laps.stimPulseLoc{nTr} = -1;
                    end
                else
                    behEventsTdt.stimOn(ind,7) = -1; % tagging stimulations
                    disp(['Trial ' num2str(nTr) ' has stimulation, but the '...
                        'stimulation is not locked to the trial']);
                    save([FileNameBase 'BTDT.mat'],'behEventsTdt');
                end
            end
            
            if(nTr == totNLaps)
               ind =  behEventsTdt.stimOn(:,4) >= trEndLfpInd;
               if(sum(ind) > 0)
                    behEventsTdt.stimOn(ind,7) = -1; % tagging stimulations
                    disp(['Trial ' num2str(nTr) ' has stimulation, but the '...
                        'stimulation is not locked to the trial']);
                    save([FileNameBase 'BTDT.mat'],'behEventsTdt');
               end
            end
        end
        
        % movie on time
        if(isfield(behEventsTdt,'movieOn'))
            ind = behEventsTdt.movieOn(:,4) >= trStartLfpInd &...
                behEventsTdt.movieOn(:,4) <= trEndLfpInd;
            Laps.movieOnT{nTr} = whlDataLfp(behEventsTdt.movieOn(ind,4),1);
            Laps.movieOnLfpInd{nTr} = behEventsTdt.movieOn(ind,4);
        else
            Laps.movieOnT{nTr} = [];
            Laps.movieOnLfpInd{nTr} = [];
        end

        %%        
        % lick period (all the licks between the start of the current trial 
        % and the start of the next trial)
        if(isfield(behEventsTdt,'lick'))
            ind = behEventsTdt.lick(:,4) >= trStartLfpInd &...
                behEventsTdt.lick(:,4) < trNextStartLfpInd;
            Laps.lickPeriodT{nTr} = whlDataLfp(behEventsTdt.lick(ind,4),1);
            Laps.lickPeriodLfpInd{nTr} = behEventsTdt.lick(ind,4);
        end
        
        %%%%%%%%%%%%%%%%%%%% N of Track (zero for the initial Track)
        Track.lapID(trStartLfpInd:trEndLfpInd,1) = nTr;
        Laps.lapID(nTr,1) = nTr;

        %%%%%%%%%%%%%%%%%%%% corrChoice - error[0]/correct[1]/not assigned[-1]
        Track.behavType(trStartLfpInd:trEndLfpInd,1) = behType(nTr);
        Track.mazeType(trStartLfpInd:trEndLfpInd,1) = behEventsTdt.mazeType(nTr);
        Track.mazeSess(trStartLfpInd:trEndLfpInd,1) = behEventsTdt.mazeSess(nTr);
        if(behType(nTr) == 1)
            Track.corrChoice(trStartLfpInd:trEndLfpInd,1) = 1;
            Laps.corrChoice(nTr,1) = 1;
        end 
    end
    
    
    %% calculate speed and acceleration for each frame (smooth first)
    fprintf('\nCalculating running speed and acceleration....\n');
    [Speed,Accel] = MazeSpeedAccel_e(whlDataLfp, lfpFreq);
    Track.speed_MMsec = Speed;
    Track.accel_MMsecSq = Accel;
    
    whlDataLfpTmp = whlDataLfp(:,[1 3]);
    [Speed1,Accel1] = MazeSpeedAccel_e(whlDataLfpTmp, lfpFreq);
    Track.speed_MMsecAll = Speed1;
    Track.accel_MMsecSqAll = Accel1;
    
    Processing = 1; % processing stage one, getting Track and Laps
    fprintf('\nTracking data saved into the structure file: %s....\n',...
            [FileNameBase '_Behav2PDataLFP.mat']);
    save([FileNameBase '_Behav2PDataLFP.mat'], ...
            'Track', 'Laps', 'Clu', 'Processing', 'sampleFreq', 'lfpFreq', 'param', '-v7.3');
end


%% Calculate speed and acceleration

function [sp,ac] = MazeSpeedAccel_e(whldata, varargin)

    % function [speed, accel] = MazeSpeedAccel(whldata,smoothlen)
    % Calculates running speed and acceleration from XY position data. 

%     [SamplRate, maxSpeed, minAccel, maxAccel, viewSpeedMap, viewAccelMap] = DefaultArgs(varargin,{1250, 100, -200, 200, 0, 0});

%     smoothlen = round(SamplRate / 10);
% 
%     % filter with a hanning window
%     hanfilter = hanning(smoothlen);
%     hanfilter = hanfilter./sum(hanfilter);

%     whldata(:,2) = Filter0(hanfilter,whldata(:,2));

    %calculate speed for values that aren't -1 or distorded by filtering
    %---------------------------------
    speeddata = diff(whldata(:,2))./diff(whldata(:,1))*1000; % mm/s
    acdata = [diff(speeddata(1:2,1));diff(speeddata(:,1))]./diff(whldata(:,1));
    sp = [speeddata(1); speeddata];
    ac = [acdata(1); acdata];
end
