function [Track,Laps] = SortTrials_Beh_smTr(FileNameBase, lfpFreq)

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
    load([FileNameBase '-whl.mat']) ;            
    whlm = length(whlDataLfp(:,2));
   
     %% load Arduino behavioral data with TDT time stamps
    load([FileNameBase 'BLFP.mat']);                       
    noLfpTrials = behEventsTdt.trialDescr(:,1)>whlm+behEventsTdt.trialDescr(1,1);
    behEventsTdt.trialDescr(noLfpTrials,:) = [];

    totNLaps = behEventsTdt.trialDescr(end,end);           % total N of laps
    
    %% Initialize output structure
    Track.xMMAll = whlDataLfp(:,2); % including the stopping period
    Track.startLfpInd = zeros(whlm,1);
    Track.lapID = zeros(whlm,1);                % order of Tracks
    Track.corrChoice = zeros(whlm,1);       % error[0]/correct[1]/not assigned[-1]
    Track.behavType = zeros(whlm,1);        % good-1, late decision-2, exploration-3, bad=not running-4
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
        % movie on time per trial
    Laps.lickPeriodT = cell(totNLaps,1);
    Laps.lickPeriodInd = cell(totNLaps,1);
        % lick period per trial
    Laps.airpuffT = cell(totNLaps,1);
    Laps.airpuffLfpInd = cell(totNLaps,1);
        % airpuff events per trial
        
    %%%%%% added by Yingxue 11/2/2019    
%     Laps.stimOnT = cell(totNLaps,1);
%     Laps.stimOnLfpInd = cell(totNLaps,1);
%         % stimulation per trial
%     Laps.stimPulseT = cell(totNLaps,1);
%     Laps.stimPulseLfpInd = cell(totNLaps,1);
%     Laps.stimPulseWidth = cell(totNLaps,1);
%     Laps.stimPulseWidthLfp = cell(totNLaps,1);
%     Laps.stimDiode = cell(totNLaps,1);
%     Laps.stimDiodeCurr = cell(totNLaps,1);
%     Laps.stimPulseMethod = cell(totNLaps,1);
%     Laps.stimPulseLoc = cell(totNLaps,1);
%         % stimulation pulses
        
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
        trStartLfpInd = trStartLfpIndArr(nTr);
        trStartLfpIndTmp = trStartLfpInd + behEventsTdt.taskDescr(1,1) -1;
        trEndLfpInd = trEndLfpIndArr(nTr);
        trEndLfpIndTmp = trEndLfpInd + behEventsTdt.taskDescr(1,1) -1;
%       figure(10); clf;      

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
        if(~isempty(behEventsTdt.pump))
            ind = behEventsTdt.pump(:,1) >= trStartLfpIndTmp &...
                behEventsTdt.pump(:,1) <= trEndLfpIndTmp;
            pumpLfpInd = behEventsTdt.pump(ind,1) - behEventsTdt.taskDescr(1,1) + 1;
            Laps.pumpT{nTr} = whlDataLfp(pumpLfpInd,1);
            Laps.pumpLfpInd{nTr} = pumpLfpInd;
        else
            Laps.pumpT{nTr} = [];
            Laps.pumpLfpInd{nTr} = [];
        end

        % lick on time
        if(~isempty(behEventsTdt.lick))
            ind = behEventsTdt.lick(:,1) >= trStartLfpIndTmp &...
                behEventsTdt.lick(:,1) <= trEndLfpIndTmp;
            lickLfpInd = behEventsTdt.lick(ind,1) - behEventsTdt.taskDescr(1,1) + 1;
            Laps.lickT{nTr} = whlDataLfp(lickLfpInd,1);
            Laps.lickLfpInd{nTr} = lickLfpInd;
        else
            Laps.lickT{nTr} = [];
            Laps.lickLfpInd{nTr} = [];
        end
       
%         % stimulation on time
%         % added on 11/2/2019 by Yingxue
%         if(isfield(behEventsTdt,'stimOn'))
%             pulseMethod = behEventsTdt.movieTDescr{nTr}(12);
%             Laps.stimPulseMethod{nTr} = pulseMethod;
%             ind = behEventsTdt.stimOn(:,1) >= trStartLfpInd &...
%                 behEventsTdt.stimOn(:,1) <= trEndLfpInd;
%             if(sum(ind) > 0)
%                 Laps.stimOnT{nTr} = whlDataLfp(behEventsTdt.stimOn(ind,1),1);
%                 Laps.stimOnLfpInd{nTr} = behEventsTdt.stimOn(ind,1);
%                 
%                 if(pulseMethod > 0)
%                     indS = behEventsTdt.stimPulse(:,1) >= trStartLfpInd &...
%                         behEventsTdt.stimPulse(:,1) <= trEndLfpInd;
%                     if(sum(indS) ~= behEventsTdt.stimOn(ind,4))
%                         indTmp = find(ind == 1);
%                         disp(['Please check stimulation no. ' num2str(indTmp) ...
%                             ', the no. of pulses does not match']);
%                     end
%                     Laps.stimPulseT{nTr} = whlDataLfp(behEventsTdt.stimPulse(indS,1),1);
%                     Laps.stimPulseLfpInd{nTr} = behEventsTdt.stimPulse(indS,1);
%                     Laps.stimPulseWidth{nTr} = behEventsTdt.stimPulse(indS,3);
%                     Laps.stimPulseWidthLfp{nTr} = behEventsTdt.stimPulse(indS,5);
%                     Laps.stimDiode{nTr} = behEventsTdt.stimOnDiode{indS}.indDiode;
%                     Laps.stimDiodeCurr{nTr} = behEventsTdt.stimOnDiode{indS}.diodeCurr; 
%                     
%                     if(pulseMethod == 3)
%                         Laps.stimPulseLoc{nTr} = behEventsTdt.movieTDescr{nTr}(9);
%                     else
%                         Laps.stimPulseLoc{nTr} = -1;
%                     end
%                 else
%                     behEventsTdt.stimOn(ind,5) = -1; % tagging stimulations
%                     disp(['Trial ' num2str(nTr) ' has stimulation, but the '...
%                         'stimulation is not locked to the trial']);
%                     save([FileNameBase 'BTDT.mat'],'behEventsTdt');
%                 end
%             end
%             
%             if(nTr == totNLaps)
%                ind =  behEventsTdt.stimOn(:,1) >= trEndLfpInd;
%                if(sum(ind) > 0)
%                     behEventsTdt.stimOn(ind,5) = -1; % tagging stimulations
%                     disp(['Trial ' num2str(nTr) ' has stimulation, but the '...
%                         'stimulation is not locked to the trial']);
%                     save([FileNameBase 'BTDT.mat'],'behEventsTdt');
%                end
%             end
%         end
        
        % movie on time
        if(isfield(behEventsTdt,'movieOn'))
            ind = behEventsTdt.movieOn(:,1) >= trStartLfpIndTmp &...
                behEventsTdt.movieOn(:,1) <= trEndLfpIndTmp;
            movieLfpInd = behEventsTdt.movieOn(ind,1) - behEventsTdt.taskDescr(1,1) + 1; 
            Laps.movieOnT{nTr} = whlDataLfp(movieLfpInd,1);
            Laps.movieOnLfpInd{nTr} = movieLfpInd;
            %%%%%% added by Yingxue on 03/30/3019
            if(LRstatus == 1)
                Laps.movieID{nTr} = behEventsTdt.movieID(ind);
            end
            %%%%%%
        else
            Laps.movieOnT{nTr} = [];
            Laps.movieOnLfpInd{nTr} = [];
            %%%%%% added by Yingxue on 03/30/3019
            if(LRstatus == 1)
                Laps.movieID{nTr} = [];
            end
            %%%%%%
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
    Track.speed_MMsecAll = Speed;
    Track.accel_MMsecSqAll = Accel;
    
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
