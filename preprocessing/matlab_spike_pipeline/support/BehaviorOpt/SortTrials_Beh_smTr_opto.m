function [Track,Laps] = SortTrials_Beh_smTr_opto(FileNameBase, lfpFreq)

% sorts maze Tracks according to different bevaioral parameters and
% calculate various behavioral parameters
%
% output: matlab file '*_Tracks.mat'
%

% cd /groups/pastalkova/home/pastalkovae/data/data/Yingxue_recordings/A632/A632_110725_01_e
% SortTrials_3armMaze_v1('A632_110725_01');

    % added 06/27/2020
    Thre_speed = 5;
    Thre_time = 4;

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

    totNLaps = behEventsTdt.trialDescr(end,3);           % total N of laps
    
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
    Track.stimOn = zeros(whlm,1);
        
    Laps.behavType = behType;
    Laps.mazeType = behEventsTdt.mazeType';
    Laps.mazeTypeMod = behEventsTdt.mazeTypeMod';
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
    Laps.stimOnT = cell(totNLaps,1);
    Laps.stimOnLfpInd = cell(totNLaps,1);
    Laps.stimOffT = cell(totNLaps,1);
    Laps.stimOffLfpInd = cell(totNLaps,1);
        % stimulation per trial
    Laps.stimTotNumPulse = cell(totNLaps,1);
    Laps.stimNumPulse = cell(totNLaps,1);
    Laps.stimPulseWidth = cell(totNLaps,1);
    Laps.stimPulsePeriod = cell(totNLaps,1);
    Laps.stimOffTime = cell(totNLaps,1);
    Laps.stimRepeats = cell(totNLaps,1);
    Laps.stimPulseMethod = cell(totNLaps,1);
    Laps.stimPulseLoc = cell(totNLaps,1);
    Laps.stimOffLfpExp = cell(totNLaps,1);
        
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
        if(nTr == totNLaps)
            a = 1;
        end
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
       
        % stimulation on time
        % added on 06/27/2020 by Yingxue
        if(isfield(behEventsTdt,'stimOn'))
            pulseMethod = behEventsTdt.movieTDescr{nTr}(12);
            Laps.stimPulseMethod{nTr} = pulseMethod;
            ind = behEventsTdt.stimOn(:,1) >= trStartLfpIndTmp &...
                behEventsTdt.stimOn(:,1) <= trEndLfpIndTmp;
            if(sum(ind) > 0)
                if(nTr > 70)
                    a = 1;
                end
                Laps.stimOnT{nTr} = whlDataLfp(behEventsTdt.stimOn(ind,1) - behEventsTdt.taskDescr(1,1) + 1,1);
                Laps.stimOnLfpInd{nTr} = behEventsTdt.stimOn(ind,1) - behEventsTdt.taskDescr(1,1) + 1;
                
                % changed on 1/1/2021 by Yingxue 
                % for stim type 7 and 8, the stimOff can be one element longer than stimOn
                if(isfield(behEventsTdt,'stimOff'))
%                         length(behEventsTdt.stimOff(:,1)) == length(behEventsTdt.stimOn(:,1))
                    Laps.stimOffT{nTr} = whlDataLfp(behEventsTdt.stimOff(ind,1) - behEventsTdt.taskDescr(1,1) + 1,1);
                    Laps.stimOffLfpInd{nTr} = behEventsTdt.stimOff(ind,1) - behEventsTdt.taskDescr(1,1) + 1;
                end
                                   
                Laps.stimTotNumPulse{nTr} = behEventsTdt.stimOn(ind,2);
                Laps.stimNumPulse{nTr} = behEventsTdt.stimOn(ind,3);
                Laps.stimPulseWidth{nTr} = behEventsTdt.stimOn(ind,4);
                Laps.stimPulsePeriod{nTr} = behEventsTdt.stimOn(ind,5);
                Laps.stimOffTime{nTr} = behEventsTdt.stimOn(ind,6);
                Laps.stimRepeats{nTr} = behEventsTdt.stimOn(ind,7);
                Laps.stimPulseMethod{nTr} = pulseMethod;

                if(pulseMethod == 3 || pulseMethod == 8) % changed on 1/1/2021 by Yingxue for stim type 7 and 8
                    Laps.stimPulseLoc{nTr} = behEventsTdt.movieTDescr{nTr}(9);
                else
                    Laps.stimPulseLoc{nTr} = -1;
                end

                indStim = find(ind == 1);
                for n = 1:length(indStim)
                    if(behEventsTdt.stimOn(indStim(n),3) == 1 && behEventsTdt.stimOn(indStim(n),7) == 1)
                        totPulseWidth = resamp(behEventsTdt.stimOn(indStim(n),4),lfpFreq);
                    elseif(behEventsTdt.stimOn(indStim(n),3) > 1 && ...
                            behEventsTdt.stimOn(indStim(n),7) <= 1) % changed on 1/1/2021 by Yingxue from == 1 to <= 1
                        totPulseWidth = behEventsTdt.stimOn(indStim(n),5)*...
                            (behEventsTdt.stimOn(indStim(n),3) - 1) + ...
                            behEventsTdt.stimOn(indStim(n),4);
                        totPulseWidth = resamp(totPulseWidth,lfpFreq);
                    else
                        totPulseWidth = (behEventsTdt.stimOn(indStim(n),5)*...
                            behEventsTdt.stimOn(indStim(n),3) + behEventsTdt.stimOn(indStim(n),6)) * ...
                            behEventsTdt.stimOn(indStim(n),7) - behEventsTdt.stimOn(indStim(n),6) - ...
                            (behEventsTdt.stimOn(indStim(n),5) - behEventsTdt.stimOn(indStim(n),4));
                        totPulseWidth = resamp(totPulseWidth,lfpFreq);
                    end

                    if((pulseMethod >= 1 && pulseMethod <= 4) || pulseMethod == 7 || pulseMethod == 8) % start cue, run onset, mid-cue, and water delivery
                        Laps.stimOffLfpExp{nTr}(n) = behEventsTdt.stimOn(indStim(n),1) ...
                            - behEventsTdt.taskDescr(1,1) + 1 ...
                            + totPulseWidth - 1;  
                    elseif(pulseMethod == 5)
                        Laps.stimOffLfpExp{nTr}(n) = behEventsTdt.stimOn(indStim(n),1) ...
                            - behEventsTdt.taskDescr(1,1) + 1 ...
                            + totPulseWidth - 1;
                        if(Laps.stimOffLfpExp{nTr}(n) > trEndLfpInd)
                            Laps.stimOffLfpExp{nTr}(n) = trEndLfpInd;
                        end
                    elseif(pulseMethod == 6)
                        Laps.stimOffLfpExp{nTr}(n) = behEventsTdt.stimOn(indStim(n),1) ...
                            - behEventsTdt.taskDescr(1,1) + 1 ...
                            + totPulseWidth - 1;
                        if(nTr < totNLaps)
                            trStartLfpIndNext = trStartLfpIndArr(nTr+1) + behEventsTdt.taskDescr(1,1) -1;
                            trEndLfpIndNext = trEndLfpIndArr(nTr+1)+ behEventsTdt.taskDescr(1,1) -1;
                            indWheelNext = behEventsTdt.wheel(:,1) >= trStartLfpIndNext & ...
                                behEventsTdt.wheel(:,1) < trEndLfpIndNext;
                            wheelEvents = behEventsTdt.wheel(indWheelNext,:);
                            
                            indWheelCur = behEventsTdt.wheel(:,1) >= trStartLfpIndTmp & ...
                                behEventsTdt.wheel(:,1) < trStartLfpIndNext;
                            wheelEventsCur = behEventsTdt.wheel(indWheelCur,:);

                            if(sum(wheelEventsCur(end-Thre_time:end,2) > Thre_speed) > Thre_time)
                                Laps.stimOffLfpExp{nTr}(n) = trStartLfpIndArr(nTr+1);
                            else
                                wheelEventsTmp = [wheelEventsCur(end-Thre_time:end,:); wheelEvents(:,:)];
                                indSpeedHigh = wheelEventsTmp(:,2) > Thre_speed;
                                % find blocks of wheel run speed that are above threshold
                                f = find(diff([0 indSpeedHigh' 0] == 1));
                                p = f(1:2:end-1); % start indices
                                y = f(2:2:end) - p; % consecutive ones' count

                                x = find(y > Thre_time,1);
                                indTmp = p(x) + Thre_time;
                                stopTime = wheelEventsTmp(indTmp,1);
                                if((stopTime - behEventsTdt.stimOn(indStim(n),1) + 1) ...
                                        < totPulseWidth)
                                    Laps.stimOffLfpExp{nTr}(n) = stopTime - behEventsTdt.taskDescr(1,1) + 1;
                                end
                            end
                        end
                    end
                    
                    stimOnLfp = behEventsTdt.stimOn(indStim(n),1)...
                            - behEventsTdt.taskDescr(1,1) + 1;
                    if(isfield(behEventsTdt,'stimOff') && ...
                            length(behEventsTdt.stimOff(:,1)) == length(behEventsTdt.stimOn(:,1)))
                        stimOffLfp = behEventsTdt.stimOff(indStim(n),1)...
                            - behEventsTdt.taskDescr(1,1) + 1;
                        Track.stimOn(stimOnLfp:stimOffLfp) = 1;
                    elseif(isfield(behEventsTdt,'stimOff') && ... % added by Yingxue on 1/1/2021, to adjust for stim condition 7 and 8
                            length(behEventsTdt.stimOff(:,1)) > length(behEventsTdt.stimOn(:,1))+1)
                        stimOffLfp = behEventsTdt.stimOff(indStim(n),1)...
                            - behEventsTdt.taskDescr(1,1) + 1;
                        Track.stimOn(stimOnLfp:stimOffLfp) = 1;
                    else
                        disp(['Tr = ' num2str(nTr) ' OnT = ' num2str(stimOnLfp) ' OffT = ' ...
                            num2str(Laps.stimOffLfpExp{nTr}(n))]);
                        Track.stimOn(stimOnLfp:....
                            Laps.stimOffLfpExp{nTr}(n)) = 1;
                    end
                end
            end
        else
            Laps.stimPulseMethod{nTr} = 0;
        end
        
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
    
    Processing = 1;
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
