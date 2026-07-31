function [Track,Laps] = SortTrials_3armMaze_smTr_pix(FileNameBase, lfpFreq)

% sorts maze Tracks according to different bevaioral parameters and
% calculate various behavioral parameters
%
% output: matlab file '*_Tracks.mat'
%

% cd /groups/pastalkova/home/pastalkovae/data/data/Yingxue_recordings/A632/A632_110725_01_e
% SortTrials_3armMaze_v1('A632_110725_01');

    % added 08/14/2018, removed _Track_Laps.mat file to reduce data duplication
    if exist([FileNameBase '_BehavElectrDataLFP.mat'], 'file') == 2
        fprintf(...
            '\nCheck whether %s file already contains Track and Laps.\n', ...
            [FileNameBase '_BehavElectrDataLFP.mat']);
        load([FileNameBase '_BehavElectrDataLFP.mat'], 'Processing');
        if(sum(Processing == 1) > 0)
            return;
        end
    end
    
    %% Load tracking file
    load([FileNameBase '-whl.mat']);            
    whlm = length(whlDataLfp(:,2));
   
     %% load Arduino behavioral data with TDT time stamps
    load([FileNameBase 'BTDT.mat']);                       
    noLfpTrials = behEventsTdt.trialDescr(:,2)>whlm;
    behEventsTdt.trialDescr(noLfpTrials,:) = [];

    totNLaps = behEventsTdt.trialDescr(end,3);           % total N of laps
    
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
        
    %%%%%%% added by Yingxue 3/30/2019
    if(isfield(behEventsTdt,'trStartCueLR')) 
        LRstatus = 1;
        Laps.trStartCueLR = zeros(totNLaps,1);
        Laps.movieID = cell(totNLaps,1);
        Laps.movieOnLPulseT = cell(totNLaps,1);
        Laps.movieOnLPulseLfpInd = cell(totNLaps,1);
        Laps.movieOffLPulseT = cell(totNLaps,1);
        Laps.movieOffLPulseLfpInd = cell(totNLaps,1);
        
        Laps.movieOnRPulseT = cell(totNLaps,1);
        Laps.movieOnRPulseLfpInd = cell(totNLaps,1);
        Laps.movieOffRPulseT = cell(totNLaps,1);
        Laps.movieOffRPulseLfpInd = cell(totNLaps,1);
    else
        LRstatus = 0;
        Laps.movieOnPulseT = cell(totNLaps,1);
        Laps.movieOnPulseLfpInd = cell(totNLaps,1);
        Laps.movieOffPulseT = cell(totNLaps,1);
        Laps.movieOffPulseLfpInd = cell(totNLaps,1);
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
            if(length(behEventsTdt.taskDescr(:,2)) == totNLaps + 1)
                trNextStartLfpInd = behEventsTdt.taskDescr(end,2);
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
        
        %%%%%% added by Yingxue on 03/30/3019
        % left or right trial
        if(LRstatus == 1)
            Laps.trStartCueLR(nTr) = behEventsTdt.trStartCueLR(nTr);
        end
        %%%%%%

        % pump on time
        ind = behEventsTdt.pump(:,2) >= trStartLfpInd &...
            behEventsTdt.pump(:,2) <= trEndLfpInd;
        Laps.pumpT{nTr} = whlDataLfp(behEventsTdt.pump(ind,2),1);
        Laps.pumpLfpInd{nTr} = behEventsTdt.pump(ind,2);

        % lick on time
        ind = behEventsTdt.lick(:,2) >= trStartLfpInd &...
            behEventsTdt.lick(:,2) <= trEndLfpInd;
        Laps.lickT{nTr} = whlDataLfp(behEventsTdt.lick(ind,2),1);
        Laps.lickLfpInd{nTr} = behEventsTdt.lick(ind,2);
        
        % airpuff on time
        if(isfield(behEventsTdt,'airpuff'))
            ind = behEventsTdt.airpuff(:,2) >= trStartLfpInd &...
                behEventsTdt.airpuff(:,2) <= trEndLfpInd;
            Laps.airpuffT{nTr} = whlDataLfp(behEventsTdt.airpuff(ind,2),1);
            Laps.airpuffLfpInd{nTr} = behEventsTdt.airpuff(ind,2);
        end
        
        % stimulation on time
        % added on 11/2/2019 by Yingxue
        if(isfield(behEventsTdt,'stimOn'))
            pulseMethod = behEventsTdt.movieTDescr{nTr}(12);
            Laps.stimPulseMethod{nTr} = pulseMethod;
            ind = behEventsTdt.stimOn(:,2) >= trStartLfpInd &...
                behEventsTdt.stimOn(:,2) <= trEndLfpInd;
            if(sum(ind) > 0)
                Laps.stimOnT{nTr} = whlDataLfp(behEventsTdt.stimOn(ind,2),1);
                Laps.stimOnLfpInd{nTr} = behEventsTdt.stimOn(ind,2);
                
                if(pulseMethod > 0)
                    indS = behEventsTdt.stimPulse(:,2) >= trStartLfpInd &...
                        behEventsTdt.stimPulse(:,2) <= trEndLfpInd;
                    if(sum(indS) ~= behEventsTdt.stimOn(ind,4))
                        indTmp = find(ind == 1);
                        disp(['Please check stimulation no. ' num2str(indTmp) ...
                            ', the no. of pulses does not match']);
                    end
                    Laps.stimPulseT{nTr} = whlDataLfp(behEventsTdt.stimPulse(indS,2),1);
                    Laps.stimPulseLfpInd{nTr} = behEventsTdt.stimPulse(indS,2);
                    Laps.stimPulseWidth{nTr} = behEventsTdt.stimPulse(indS,4);
                    Laps.stimPulseWidthLfp{nTr} = behEventsTdt.stimPulse(indS,6);
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
                    behEventsTdt.stimOn(ind,5) = -1; % tagging stimulations
                    disp(['Trial ' num2str(nTr) ' has stimulation, but the '...
                        'stimulation is not locked to the trial']);
                    save([FileNameBase 'BTDT.mat'],'behEventsTdt');
                end
            end
            
            if(nTr == totNLaps)
               ind =  behEventsTdt.stimOn(:,2) >= trEndLfpInd;
               if(sum(ind) > 0)
                    behEventsTdt.stimOn(ind,5) = -1; % tagging stimulations
                    disp(['Trial ' num2str(nTr) ' has stimulation, but the '...
                        'stimulation is not locked to the trial']);
                    save([FileNameBase 'BTDT.mat'],'behEventsTdt');
               end
            end
        end
        
        % movie on time
        if(isfield(behEventsTdt,'movieOn'))
            ind = behEventsTdt.movieOn(:,2) >= trStartLfpInd &...
                behEventsTdt.movieOn(:,2) <= trEndLfpInd;
            Laps.movieOnT{nTr} = whlDataLfp(behEventsTdt.movieOn(ind,2),1);
            Laps.movieOnLfpInd{nTr} = behEventsTdt.movieOn(ind,2);
            %%%%%% added by Yingxue on 03/30/2019
            if(LRstatus == 1)
                Laps.movieID{nTr} = behEventsTdt.movieID(ind);
            end
            %%%%%%
            %%%%%% added by Yingxue on 03/11/2020
            if(behEventsTdt.mazeType(nTr) == 23)
                Laps.movieLocation(nTr) = behEventsTdt.taskDescr(nTr,5);
            end
            %%%%%%
        else
            Laps.movieOnT{nTr} = [];
            Laps.movieOnLfpInd{nTr} = [];
            %%%%%% added by Yingxue on 03/30/2019
            if(LRstatus == 1)
                Laps.movieID{nTr} = [];
            end
            %%%%%%
            %%%%%% added by Yingxue on 03/11/2020
            if(behEventsTdt.mazeType(nTr) == 23)
                Laps.movieLocation(nTr) = behEventsTdt.taskDescr(nTr,5);
            end
            %%%%%%
        end
        
        if(LRstatus == 0)
            % movie pulse on time
            if(isfield(behEventsTdt,'TDTMovieOn'))
                ind = behEventsTdt.TDTMovieOn(:,2) >= trStartLfpInd &...
                    behEventsTdt.TDTMovieOn(:,2) <= trEndLfpInd;
                Laps.movieOnPulseT{nTr} = ...
                    whlDataLfp(behEventsTdt.TDTMovieOn(ind,2),1);
                Laps.movieOnPulseLfpInd{nTr} = behEventsTdt.TDTMovieOn(ind,2);
            else
                Laps.movieOnPulseT{nTr} = [];
                Laps.movieOnPulseLfpInd{nTr} = [];
            end

            % movie pulse off time
            if(isfield(behEventsTdt,'TDTMovieOff'))
                ind = behEventsTdt.TDTMovieOff(:,2) >= trStartLfpInd &...
                    behEventsTdt.TDTMovieOff(:,2) <= trEndLfpInd;
                Laps.movieOffPulseT{nTr} = ...
                    whlDataLfp(behEventsTdt.TDTMovieOff(ind,2),1);
                Laps.movieOffPulseLfpInd{nTr} = behEventsTdt.TDTMovieOff(ind,2);
            else
                Laps.movieOffPulseT{nTr} = [];
                Laps.movieOffPulseLfpInd{nTr} = [];
            end
        else
            %%
            %%%%% added by Yingxue on 03/30/3019
            % movie pulse on time left
            if(isfield(behEventsTdt,'TDTMovieOnL'))
                ind = behEventsTdt.TDTMovieOnL(:,2) >= trStartLfpInd &...
                    behEventsTdt.TDTMovieOnL(:,2) <= trEndLfpInd;
                Laps.movieOnLPulseT{nTr} = ...
                    whlDataLfp(behEventsTdt.TDTMovieOnL(ind,2),1);
                Laps.movieOnLPulseLfpInd{nTr} = behEventsTdt.TDTMovieOnL(ind,2);
            else
                Laps.movieOnLPulseT{nTr} = [];
                Laps.movieOnLPulseLfpInd{nTr} = [];
            end

            % movie pulse off time left
            if(isfield(behEventsTdt,'TDTMovieOffL'))
                ind = behEventsTdt.TDTMovieOffL(:,2) >= trStartLfpInd &...
                    behEventsTdt.TDTMovieOffL(:,2) <= trEndLfpInd;
                Laps.movieOffLPulseT{nTr} = ...
                    whlDataLfp(behEventsTdt.TDTMovieOffL(ind,2),1);
                Laps.movieOffLPulseLfpInd{nTr} = behEventsTdt.TDTMovieOffL(ind,2);
            else
                Laps.movieOffLPulseT{nTr} = [];
                Laps.movieOffLPulseLfpInd{nTr} = [];
            end

            % movie pulse on time left
            if(isfield(behEventsTdt,'TDTMovieOnR'))
                ind = behEventsTdt.TDTMovieOnR(:,2) >= trStartLfpInd &...
                    behEventsTdt.TDTMovieOnR(:,2) <= trEndLfpInd;
                Laps.movieOnRPulseT{nTr} = ...
                    whlDataLfp(behEventsTdt.TDTMovieOnR(ind,2),1);
                Laps.movieOnRPulseLfpInd{nTr} = behEventsTdt.TDTMovieOnR(ind,2);
            else
                Laps.movieOnRPulseT{nTr} = [];
                Laps.movieOnRPulseLfpInd{nTr} = [];
            end
            
            % movie pulse off time left
            if(isfield(behEventsTdt,'TDTMovieOffR'))
                ind = behEventsTdt.TDTMovieOffR(:,2) >= trStartLfpInd &...
                    behEventsTdt.TDTMovieOffR(:,2) <= trEndLfpInd;
                Laps.movieOffRPulseT{nTr} = ...
                    whlDataLfp(behEventsTdt.TDTMovieOffR(ind,2),1);
                Laps.movieOffRPulseLfpInd{nTr} = behEventsTdt.TDTMovieOffR(ind,2);
            else
                Laps.movieOffRPulseT{nTr} = [];
                Laps.movieOffRPulseLfpInd{nTr} = [];
            end
            %%%%% added by Yingxue on 03/30/3019
        end
        %%
        
        % lick period (all the licks between the start of the current trial 
        % and the start of the next trial)
        if(isfield(behEventsTdt,'lick'))
            ind = behEventsTdt.lick(:,2) >= trStartLfpInd &...
                behEventsTdt.lick(:,2) < trNextStartLfpInd;
            Laps.lickPeriodT{nTr} = whlDataLfp(behEventsTdt.lick(ind,2),1);
            Laps.lickPeriodLfpInd{nTr} = behEventsTdt.lick(ind,2);
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
            [FileNameBase '_BehavElectrDataLFP.mat']);
    save([FileNameBase '_BehavElectrDataLFP.mat'], ...
            'Track', 'Laps', 'Processing', '-v7.3');
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
