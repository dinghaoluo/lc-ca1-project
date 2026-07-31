function [Track,Laps] = SortTrials_3armMaze_mpfi(FileNameBase, lfpFreq, purpose)

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
    Track.xMM = whlDataLfp(:,2);
    Laps.behavType = behType;
    whlm = length(Track.xMM);
    
    load([FileNameBase '-segLen.mat']);
    totalTrackLen = seg1Len + seg2Len + seg3Len;
    
    % cut the Arduino trial description table to the length of dat file (= length of whlDataLfp)
    % behEventsTdt.trialDescr
    if(purpose == 1 || purpose == 4) % treadmill running
        load([FileNameBase 'BTDT.mat']);                       % Arduino behavioral data with TDT time stamps
        noLfpTrials = behEventsTdt.trialDescr(:,2)>whlm;
        behEventsTdt.trialDescr(noLfpTrials,:) = [];

        totNLaps = behEventsTdt.trialDescr(end,3);           % total N of laps
        trialList = behEventsTdt.trialDescr(:,3);
        
    end
    
    %% Initialize output structure
    Track.startLfpInd = zeros(whlm,1);
    Track.lapID = zeros(whlm,1);                % order of Tracks
    Track.corrChoice = zeros(whlm,1);       % error[0]/correct[1]/not assigned[-1]
    Track.behavType = zeros(whlm,1);        % good-1, late decision-2, exploration-3, bad=not running-4
    Track.speed_MMsec = zeros(whlm,1);
    Track.accel_MMsecSq = zeros(whlm,1);
        
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
    
  
    %% Analyze individual trials
    % use behEventsTdt.TrackDescr to get the beginning of each Track
    % behEventsTdt.TrackDescr: total time, session time, Track ID, ????, correct/incorrect
    % correct for repeated trail IDs
    
    if (purpose == 1 || purpose == 4) % not opened field
        for nTr = 1 : totNLaps
            nLine = find(trialList == nTr);
            trStartLfpInd = trStartLfpIndArr(nTr);
            trEndLfpInd = trEndLfpIndArr(nTr);
%             figure(10); clf;
            
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
                while(Track.xMM(trEndLfpInd) < totalTrackLen - 50)
                    trEndLfpInd = trEndLfpInd-1;
                    shift = shift+1;
                end
                if(shift ~= 0)
                    disp(['Trial no. ' num2str(nTr)...
                        ': the stop lfp index is shifted backward by '...
                        num2str(shift) ' points.']);
                else
                    while(Track.xMM(trEndLfpInd+1) > totalTrackLen - 1)
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
            ind = behEventsTdt.pump(:,2) >= trStartLfpInd &...
                behEventsTdt.pump(:,2) <= trEndLfpInd;
            Laps.pumpT{nTr} = whlDataLfp(behEventsTdt.pump(ind,2),1);
            Laps.pumpLfpInd{nTr} = behEventsTdt.pump(ind,2);

            % lick on time
            ind = behEventsTdt.lick(:,2) >= trStartLfpInd &...
                behEventsTdt.lick(:,2) <= trEndLfpInd;
            Laps.lickT{nTr} = whlDataLfp(behEventsTdt.lick(ind,2),1);
            Laps.lickLfpInd{nTr} = behEventsTdt.lick(ind,2);
            
            %%%%%%%%%%%%%%%%%%%% N of Track (zero for the initial Track)
            Track.lapID(trStartLfpInd:trEndLfpInd,1) = nTr;
            Laps.lapID(nTr,1) = nTr;
            
            %%%%%%%%%%%%%%%%%%%% corrChoice - error[0]/correct[1]/not assigned[-1]
            if(purpose == 4)
                Track.behavType(trStartLfpInd:trEndLfpInd,1) = behType(nTr);
                if(behType(nTr) == 1)
                    Track.corrChoice(trStartLfpInd:trEndLfpInd,1) = 1;
                    Laps.corrChoice(nTr,1) = 1;
                end
            end  
        end
    end
       
    %% calculate speed and acceleration for each frame (smooth first)
    fprintf('\nCalculating running speed and acceleration....\n');
    [Speed,Accel] = MazeSpeedAccel_e(whlDataLfp, lfpFreq);
    Track.speed_MMsec = Speed;
    Track.accel_MMsecSq = Accel;
    
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
