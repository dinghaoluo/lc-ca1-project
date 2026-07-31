function Arduino2TDTtime_smTr(baseFileName,sampleFreq,lfpFreq)
% using the sync signal to align the arduino time and the TDT time

%%%%%%%%% check arguments
if nargin<3
    disp('At least three arguments are needed for this function.');
    return;
end

%%%%%%%%% load recording file
fullPathB = [baseFileName 'B.mat'];

if(exist(fullPathB))
    load(fullPathB);
else
    disp('Behavioral event file does not exist.')
    return;
end

if(~isfield(behEvents,'ArdSyncMsec'))
    ArdSyncMsec = [];
else
    ArdSyncMsec = behEvents.ArdSyncMsec;
end
ardSyncLen = size(ArdSyncMsec,1);
tdtSyncLen = length(behEvents.TDTsyncInd);
TDTSyncMsecTmp = behEvents.TDTsyncMsec;
if(ardSyncLen > tdtSyncLen)
    TDTSyncMsecTmp = [zeros(ardSyncLen-tdtSyncLen,1); TDTSyncMsecTmp];
end

%%%%%%%% convert the arduino time stamps to the TDT time stamps
behEventsTdt = [];

%%
% behEvents.taskDescr
% first column: time stamps sampled at sampleFreq
% second column: time stamps sampled at lfpFreq
if(isfield(behEvents,'taskDescr'))
    
    % check whether there is time reverse
    diffTm = diff(behEvents.taskDescr(:,1));
    indNeg = find(diffTm < 0);
    if(~isempty(indNeg))
        behEvents.taskDescr(1:indNeg,1) = 0;
    end
    
    indLast = 1;
    for i = 1:size(behEvents.taskDescr,1)
        % find the sync event immediately after the trial description event
        indTmp = find(ArdSyncMsec(indLast:end,1) ...
            >= behEvents.taskDescr(i,1));
        if(~isempty(indTmp))
            indLast = indTmp(1)+indLast-1;
            if(TDTSyncMsecTmp(indLast) == 0)
                behEventsTdt.taskDescr(i,1) = 0;
                continue;
            end
            if(indLast ~= 1) % convert the arduino time to TDT time
                behEventsTdt.taskDescr(i,1) = TDTSyncMsecTmp(indLast-1) ...
                    + (TDTSyncMsecTmp(indLast) - TDTSyncMsecTmp(indLast-1))...
                    /(ArdSyncMsec(indLast,1) - ArdSyncMsec(indLast-1,1))...
                    *(behEvents.taskDescr(i,1)-ArdSyncMsec(indLast-1,1));
                if(behEventsTdt.taskDescr(i,1) < 0)
                    behEventsTdt.taskDescr(i,1) = 0;
                end
            else
                behEventsTdt.taskDescr(i,1) = TDTSyncMsecTmp(indLast)...
                    - (TDTSyncMsecTmp(indLast+1) - TDTSyncMsecTmp(indLast))...
                    /(ArdSyncMsec(indLast+1,1) - ArdSyncMsec(indLast,1))...
                    *(ArdSyncMsec(indLast,1) - behEvents.taskDescr(i,1));
            end
        else
            behEventsTdt.taskDescr(i,1) = TDTSyncMsecTmp(end)...
                + (TDTSyncMsecTmp(end) - TDTSyncMsecTmp(end-1))...
                /(ArdSyncMsec(end,1) - ArdSyncMsec(end-1,1))...
                *(behEvents.taskDescr(i,1)-ArdSyncMsec(end,1));
        end
    end
    behEventsTdt.taskDescr(:,2) = resamp(behEventsTdt.taskDescr(:,1),lfpFreq);
    behEventsTdt.taskDescr(:,1) = resamp(behEventsTdt.taskDescr(:,1),sampleFreq);
    behEventsTdt.taskDescr(:,3:4) = behEvents.taskDescr(:,2:3);
end

% behEvents.movieTDescr
% first column: time stamps sampled at sampleFreq
% second column: time stamps sampled at lfpFreq
if(isfield(behEvents,'movieTDescr'))
    behEventsTdt.movieTDescr = behEvents.movieTDescr;
end

% behEvents.trialDescr 
% first column: time stamps sampled at sampleFreq
% second column: time stamps sampled at lfpFreq
if(isfield(behEvents,'trialDescr'))
     % check whether there is time reverse
    diffTm = diff(behEvents.trialDescr(:,1));
    indNeg = find(diffTm < 0);
    if(~isempty(indNeg))
        behEvents.trialDescr(1:indNeg,1) = 0;
    end

    indLast = 1;
    for i = 1:size(behEvents.trialDescr,1)
       % find the sync event immediately after the trial description event
       indTmp = find(ArdSyncMsec(indLast:end,1) ...
                        >= behEvents.trialDescr(i,1));
       if(~isempty(indTmp))
           indLast = indTmp(1)+indLast-1;
           if(TDTSyncMsecTmp(indLast) == 0)
               behEventsTdt.trialDescr(i,1) = 0;
               continue;
           end
           if(indLast ~= 1) % convert the arduino time to TDT time
                behEventsTdt.trialDescr(i,1) = TDTSyncMsecTmp(indLast-1) ...
                    + (TDTSyncMsecTmp(indLast) - TDTSyncMsecTmp(indLast-1))...
                    /(ArdSyncMsec(indLast,1) - ArdSyncMsec(indLast-1,1))...
                    *(behEvents.trialDescr(i,1)-ArdSyncMsec(indLast-1,1));
                if(behEventsTdt.trialDescr(i,1) < 0)
                    behEventsTdt.trialDescr(i,1) = 0;
                end               
           else
               behEventsTdt.trialDescr(i,1) = TDTSyncMsecTmp(indLast)...
                   - (TDTSyncMsecTmp(indLast+1) - TDTSyncMsecTmp(indLast))...
                   /(ArdSyncMsec(indLast+1,1) - ArdSyncMsec(indLast,1))...
                   *(ArdSyncMsec(indLast,1) - behEvents.trialDescr(i,1));
           end
       else
           behEventsTdt.trialDescr(i,1) = TDTSyncMsecTmp(end)...
               + (TDTSyncMsecTmp(end) - TDTSyncMsecTmp(end-1))...
               /(ArdSyncMsec(end,1) - ArdSyncMsec(end-1,1))...
               *(behEvents.trialDescr(i,1)-ArdSyncMsec(end,1));
       end
    end
    behEventsTdt.trialDescr(:,2) = resamp(behEventsTdt.trialDescr(:,1),lfpFreq);
    behEventsTdt.trialDescr(:,1) = resamp(behEventsTdt.trialDescr(:,1),sampleFreq);
    behEventsTdt.trialDescr(:,3) = behEvents.trialDescr(:,3);
    behEventsTdt.trialDescr(:,4) = behEvents.trialDescr(:,2); % performance
end

%%
% behEvents.lick
% first column: time stamps sampled at sampleFreq
% second column: time stamps sampled at lfpFreq
if(isfield(behEvents,'lick'))
    % check whether there is time reverse, if so, delete the events
    % before the time reverse
    diffTm = diff(behEvents.lick(:,1));
    indNeg = find(diffTm < 0);
    if(~isempty(indNeg))
        behEvents.lick = behEvents.lick(indNeg(end)+1:end,:);
    end
    
    indLast = 1;
    for i = 1:size(behEvents.lick,1)
        indTmp = find(ArdSyncMsec(indLast:end,1) >= behEvents.lick(i,1));
        if(~isempty(indTmp))
            indLast = indTmp(1)+indLast-1;
            if(TDTSyncMsecTmp(indLast) == 0)
                behEventsTdt.lick(i,1) = 0;
                continue;
            end
            if(indLast ~= 1) % convert the arduino time to TDT time
                behEventsTdt.lick(i,1) = TDTSyncMsecTmp(indLast-1)...
                    + (TDTSyncMsecTmp(indLast) - TDTSyncMsecTmp(indLast-1))...
                    /(ArdSyncMsec(indLast,1) - ArdSyncMsec(indLast-1,1))...
                    *(behEvents.lick(i,1)-ArdSyncMsec(indLast-1,1));
                if(behEventsTdt.lick(i,1) < 0)
                    behEventsTdt.lick(i,1) = 0;
                end
            else
                behEventsTdt.lick(i,1) = TDTSyncMsecTmp(indLast) ...
                    - (TDTSyncMsecTmp(indLast+1) - TDTSyncMsecTmp(indLast))...
                    /(ArdSyncMsec(indLast+1,1) - ArdSyncMsec(indLast,1))...
                    *(ArdSyncMsec(indLast,1) - behEvents.lick(i,1));
            end
        else
            behEventsTdt.lick(i,1) = TDTSyncMsecTmp(end) ...
                + (TDTSyncMsecTmp(end) - TDTSyncMsecTmp(end-1))...
                /(ArdSyncMsec(end,1) - ArdSyncMsec(end-1,1))...
                *(behEvents.lick(i,1)-ArdSyncMsec(end,1));
        end
    end
    behEventsTdt.lick(:,2) = resamp(behEventsTdt.lick(:,1),lfpFreq);
    behEventsTdt.lick(:,1) = resamp(behEventsTdt.lick(:,1),sampleFreq);
    behEventsTdt.lick(:,3:4) = behEvents.lick(:,2:3);
end

%%
%%%%%% added by Yingxue 11/2/2019  
% behEvents.stimOn
% first column: time stamps sampled at sampleFreq
% second column: time stamps sampled at lfpFreq
if(isfield(behEvents,'stimOn'))
    % exclude stimulations that are not recorded in the behavioral file
    % which are very likely tagging stimulations
    ind = behEvents.stimOn(:,1) ~= -1;
    % check whether there is time reverse
    diffTm = diff(behEvents.stimOn(ind,1));
    indNeg = find(diffTm < 0);
    if(~isempty(indNeg))
        behEvents.stimOn = behEvents.stimOn(indNeg(end)+1:end,:);
        behEvents.stimOnPar = behEvents.stimOnPar(indNeg(end)+1:end,:);
    end
    
    indLast = 1;
    indStimP = 0;
    for i = 1:size(behEvents.stimOn,1)
        if(behEvents.stimOn(i,1) ~= -1)
            % find the sync event immediately after the trial description event
            indTmp = find(ArdSyncMsec(indLast:end,1) ...
                >= behEvents.stimOn(i,1));
            if(~isempty(indTmp))
                indLast = indTmp(1)+indLast-1;
                if(TDTSyncMsecTmp(indLast) == 0)
                    behEventsTdt.stimOn(i,1) = 0;
                    continue;
                end
                if(indLast ~= 1) % convert the arduino time to TDT time
                    behEventsTdt.stimOn(i,1) = TDTSyncMsecTmp(indLast-1) ...
                        + (TDTSyncMsecTmp(indLast) - TDTSyncMsecTmp(indLast-1))...
                        /(ArdSyncMsec(indLast,1) - ArdSyncMsec(indLast-1,1))...
                        *(behEvents.stimOn(i,1)-ArdSyncMsec(indLast-1,1));
                    if(behEventsTdt.stimOn(i,1) < 0)
                        behEventsTdt.stimOn(i,1) = 0;
                    end
                else
                    behEventsTdt.stimOn(i,1) = TDTSyncMsecTmp(indLast)...
                        - (TDTSyncMsecTmp(indLast+1) - TDTSyncMsecTmp(indLast))...
                        /(ArdSyncMsec(indLast+1,1) - ArdSyncMsec(indLast,1))...
                        *(ArdSyncMsec(indLast,1) - behEvents.stimOn(i,1));
                end
            else
                behEventsTdt.stimOn(i,1) = TDTSyncMsecTmp(end)...
                    + (TDTSyncMsecTmp(end) - TDTSyncMsecTmp(end-1))...
                    /(ArdSyncMsec(end,1) - ArdSyncMsec(end-1,1))...
                    *(behEvents.stimOn(i,1)-ArdSyncMsec(end,1));
            end            
        end
            
        indStim = behEvents.stimOn(i,3);      
        for n = 1:behEvents.stimOn(i,4)
            behEventsTdt.stimPulse(indStimP+n,3) = i; % stimulation no.
            behEventsTdt.stimPulse(indStimP+n,4) = behEvents.stimOnPar(i,2); % pulse width
            behEventsTdt.stimPulse(indStimP+n,5) = n; % num within stim

            behEventsTdt.stimPulse(indStimP+n,1) = behEvents.TDTstimMsec(indStimP+n);
            behEventsTdt.stimPulse(indStimP+n,6) = behEvents.TDTpulseWidthMsec(indStimP+n);
               
            behEventsTdt.stimOnDiode{indStimP+n} = behEvents.stimOnDiode{i};
            
            if(behEvents.stimOn(i,1) == -1 && n == 1)
                behEventsTdt.stimOn(i,1) = behEventsTdt.stimPulse(indStimP+n,1);
            end
            indStim = indStim + 1;
        end
        
        indStimP = indStimP + behEvents.stimOn(i,4);
    end
    behEventsTdt.stimOn(:,2) = resamp(behEventsTdt.stimOn(:,1),lfpFreq);
    behEventsTdt.stimOn(:,1) = resamp(behEventsTdt.stimOn(:,1),sampleFreq);
    behEventsTdt.stimOn(:,3:4) = behEvents.stimOn(:,3:4);  
    % changed on 2/17/2010 because there is an error in getStim_mpfi_smTr line 21 (no -1 in column 5)
    indTag = behEvents.stimOn(:,1) < 0;
    stimIsTag = ones(length(behEvents.stimOn(:,1)),1);
    stimIsTag(indTag) = -1;
    behEventsTdt.stimOn(:,5) = stimIsTag;  
    behEventsTdt.stimOn(:,6:9) = behEvents.stimOnPar(:,2:end);
    %
    
    behEventsTdt.stimPulse(:,2) = resamp(behEventsTdt.stimPulse(:,1),lfpFreq);
    behEventsTdt.stimPulse(:,1) = resamp(behEventsTdt.stimPulse(:,1),sampleFreq);
    behEventsTdt.stimPulse(:,7) = resamp(behEventsTdt.stimPulse(:,6),lfpFreq);
    behEventsTdt.stimPulse(:,6) = resamp(behEventsTdt.stimPulse(:,6),sampleFreq);    
end

%% added on 1/21/2019
% trial start time is defined by taskDescr, and end time is defined by trialDescr
if(length(behEventsTdt.taskDescr(:,1)) == length(behEventsTdt.trialDescr(:,1)))
    if(sum(behEventsTdt.trialDescr(:,1)-behEventsTdt.taskDescr(:,1) <= 0) > 0)
        behEventsTdt.trialDescr = behEventsTdt.trialDescr(2:end,:);
        behEventsTdt.trialDescr(:,3) = behEventsTdt.trialDescr(:,3)-1;
        disp('behEventsTdt.trialDescr: remove the first NT event');
    end
elseif(length(behEventsTdt.taskDescr(:,1))-1 == length(behEventsTdt.trialDescr(:,1)))
    if(sum(behEventsTdt.trialDescr(:,1)-behEventsTdt.taskDescr(1:end-1,1) <= 0) > 0)
        disp('Error: trial end time <= trial start time');
        return;
    end
else
    disp('Error: unequal number of trial start time and end time');
    return;
end

%%
% behEvents.trialT
% first column: time stamps sampled at sampleFreq
% second column: time stamps sampled at lfpFreq
if(isfield(behEvents,'trialT'))
    indLast = 1;
    indTmp = find(ArdSyncMsec(indLast:end,1) >= behEvents.trialT(1));
    if(~isempty(indTmp))
        indLast = indTmp(1)+indLast-1;
        if(TDTSyncMsecTmp(indLast) == 0)
            behEventsTdt.trialT(1) = 0;
        else
            if(indLast ~= 1)
                behEventsTdt.trialT(1) = TDTSyncMsecTmp(indLast-1)...
                    + (TDTSyncMsecTmp(indLast) - TDTSyncMsecTmp(indLast-1))...
                    /(ArdSyncMsec(indLast,1) - ArdSyncMsec(indLast-1,1))...
                    *(behEvents.trialT(i,1)-ArdSyncMsec(indLast-1,1));
                if(behEventsTdt.trialT(1) < 0)
                    behEventsTdt.trialT(1) = 0;
                end
            else
                behEventsTdt.trialT(1) = TDTSyncMsecTmp(indLast) ...
                    - (TDTSyncMsecTmp(indLast+1) - TDTSyncMsecTmp(indLast))...
                    /(ArdSyncMsec(indLast+1,1) - ArdSyncMsec(indLast,1))...
                    *(ArdSyncMsec(indLast,1) - behEvents.trialT(1));
            end
        end
    else
        behEventsTdt.trialT(1) = TDTSyncMsecTmp(end) ...
            + (TDTSyncMsecTmp(end) - TDTSyncMsecTmp(end-1))...
            /(ArdSyncMsec(end,1) - ArdSyncMsec(end-1,1))...
            *(behEvents.trialT(1)-ArdSyncMsec(end,1));
    end
    behEventsTdt.trialT(2) = resamp(behEventsTdt.trialT(1),lfpFreq);
    behEventsTdt.trialT(1) = resamp(behEventsTdt.trialT(1),sampleFreq);
    behEventsTdt.trialT(3) = behEvents.trialT(2);
end

%%
% behEvents.wheel
% first column: time stamps sampled at sampleFreq
% second column: time stamps sampled at lfpFreq
if(isfield(behEvents,'wheel'))
    % check whether there is time reverse, if so, delete the events
    % before the time reverse
    diffTm = diff(behEvents.wheel(:,1));
    indNeg = find(diffTm < 0);
    if(~isempty(indNeg))
        behEvents.wheel = behEvents.wheel(indNeg(end)+1:end,:);
    end
    
    indLast = 1;
    for i = 1:size(behEvents.wheel,1)
        indTmp = find(ArdSyncMsec(indLast:end,1) >= behEvents.wheel(i,1));
        if(~isempty(indTmp))
            indLast = indTmp(1)+indLast-1;
            if(TDTSyncMsecTmp(indLast) == 0)
                behEventsTdt.wheel(i,1) = 0;
                continue;
            end
            if(indLast ~= 1) % convert the arduino time to TDT time
                behEventsTdt.wheel(i,1) = TDTSyncMsecTmp(indLast-1) ...
                    + (TDTSyncMsecTmp(indLast) - TDTSyncMsecTmp(indLast-1))...
                    /(ArdSyncMsec(indLast,1) - ArdSyncMsec(indLast-1,1))...
                    *(behEvents.wheel(i,1)-ArdSyncMsec(indLast-1,1));
                if(behEventsTdt.wheel(i,1) < 0)
                    behEventsTdt.wheel(i,1) = 0;
                end
            else
                behEventsTdt.wheel(i,1) = TDTSyncMsecTmp(indLast) ...
                    - (TDTSyncMsecTmp(indLast+1) - TDTSyncMsecTmp(indLast))...
                    /(ArdSyncMsec(indLast+1,1) - ArdSyncMsec(indLast,1))...
                    *(ArdSyncMsec(indLast,1) - behEvents.wheel(i,1));
            end
        else
            behEventsTdt.wheel(i,1) = TDTSyncMsecTmp(end) ...
                + (TDTSyncMsecTmp(end) - TDTSyncMsecTmp(end-1))...
                /(ArdSyncMsec(end,1) - ArdSyncMsec(end-1,1))...
                *(behEvents.wheel(i,1)-ArdSyncMsec(end,1));
        end
    end
    behEventsTdt.wheel(:,2) = resamp(behEventsTdt.wheel(:,1),lfpFreq);
    behEventsTdt.wheel(:,1) = resamp(behEventsTdt.wheel(:,1),sampleFreq);
    behEventsTdt.wheel(:,3:4) = behEvents.wheel(:,2:3);
end

%%
% behEvents.pump
% first column: time stamps sampled at sampleFreq
% second column: time stamps sampled at lfpFreq
if(isfield(behEvents,'pump'))
    % check whether there is time reverse, if so, delete the events
    % before the time reverse
    diffTm = diff(behEvents.pump(:,1));
    indNeg = find(diffTm < 0);
    if(~isempty(indNeg))
        behEvents.pump = behEvents.pump(indNeg(end)+1:end,:);
    end
    
    indLast = 1;
    for i = 1:size(behEvents.pump,1)
        indTmp = find(ArdSyncMsec(indLast:end,1) >= behEvents.pump(i,1));
        if(~isempty(indTmp))
            indLast = indTmp(1)+indLast-1;
            if(TDTSyncMsecTmp(indLast) == 0)
                behEventsTdt.pump(i,1) = 0;
                continue;
            end
            if(indLast ~= 1) % convert the arduino time to TDT time
                behEventsTdt.pump(i,1) = TDTSyncMsecTmp(indLast-1) ...
                    + (TDTSyncMsecTmp(indLast) - TDTSyncMsecTmp(indLast-1))...
                    /(ArdSyncMsec(indLast,1) - ArdSyncMsec(indLast-1,1))...
                    *(behEvents.pump(i,1)-ArdSyncMsec(indLast-1,1));
                if(behEventsTdt.pump(i,1) < 0)
                    behEventsTdt.pump(i,1) = 0;
                end
            else
                behEventsTdt.pump(i,1) = TDTSyncMsecTmp(indLast) ...
                    - (TDTSyncMsecTmp(indLast+1) - TDTSyncMsecTmp(indLast))...
                    /(ArdSyncMsec(indLast+1,1) - ArdSyncMsec(indLast,1))...
                    *(ArdSyncMsec(indLast,1) - behEvents.pump(i,1));
            end
        else
            behEventsTdt.pump(i,1) = TDTSyncMsecTmp(end)...
                + (TDTSyncMsecTmp(end) - TDTSyncMsecTmp(end-1))...
                /(ArdSyncMsec(end,1) - ArdSyncMsec(end-1,1))...
                *(behEvents.pump(i,1)-ArdSyncMsec(end,1));
        end
    end
    behEventsTdt.pump(:,2) = resamp(behEventsTdt.pump(:,1),lfpFreq);
    behEventsTdt.pump(:,1) = resamp(behEventsTdt.pump(:,1),sampleFreq);
    behEventsTdt.pump(:,3:4) = behEvents.pump(:,2:3);
end

%%
% behEvents.airpuff
% first column: time stamps sampled at sampleFreq
% second column: time stamps sampled at lfpFreq
if(isfield(behEvents,'airpuff'))
    % check whether there is time reverse, if so, delete the events
    % before the time reverse
    diffTm = diff(behEvents.airpuff(:,1));
    indNeg = find(diffTm < 0);
    if(~isempty(indNeg))
        behEvents.airpuff = behEvents.airpuff(indNeg(end)+1:end,:);
    end
    
    indLast = 1;
    for i = 1:size(behEvents.airpuff,1)
        indTmp = find(ArdSyncMsec(indLast:end,1) >= behEvents.airpuff(i,1));
        if(~isempty(indTmp))
            indLast = indTmp(1)+indLast-1;
            if(TDTSyncMsecTmp(indLast) == 0)
                behEventsTdt.airpuff(i,1) = 0;
                continue;
            end
            if(indLast ~= 1) % convert the arduino time to TDT time
                behEventsTdt.airpuff(i,1) = TDTSyncMsecTmp(indLast-1)...
                    + (TDTSyncMsecTmp(indLast) - TDTSyncMsecTmp(indLast-1))...
                    /(ArdSyncMsec(indLast,1) - ArdSyncMsec(indLast-1,1))...
                    *(behEvents.airpuff(i,1)-ArdSyncMsec(indLast-1,1));
                if(behEventsTdt.airpuff(i,1) < 0)
                    behEventsTdt.airpuff(i,1) = 0;
                end
            else
                behEventsTdt.airpuff(i,1) = TDTSyncMsecTmp(indLast) ...
                    - (TDTSyncMsecTmp(indLast+1) - TDTSyncMsecTmp(indLast))...
                    /(ArdSyncMsec(indLast+1,1) - ArdSyncMsec(indLast,1))...
                    *(ArdSyncMsec(indLast,1) - behEvents.airpuff(i,1));
            end
        else
            behEventsTdt.airpuff(i,1) = TDTSyncMsecTmp(end) ...
                + (TDTSyncMsecTmp(end) - TDTSyncMsecTmp(end-1))...
                /(ArdSyncMsec(end,1) - ArdSyncMsec(end-1,1))...
                *(behEvents.airpuff(i,1)-ArdSyncMsec(end,1));
        end
    end
    behEventsTdt.airpuff(:,2) = resamp(behEventsTdt.airpuff(:,1),lfpFreq);
    behEventsTdt.airpuff(:,1) = resamp(behEventsTdt.airpuff(:,1),sampleFreq);
    behEventsTdt.airpuff(:,3) = behEvents.airpuff(:,2);
end


%%
% behEvents.movieOn
% first column: time stamps sampled at sampleFreq
% second column: time stamps sampled at lfpFreq
if(isfield(behEvents,'movieOn'))
    % check whether there is time reverse, if so, delete the events
    % before the time reverse
    diffTm = diff(behEvents.movieOn(:,1));
    indNeg = find(diffTm < 0);
    if(~isempty(indNeg))
        behEvents.movieOn = behEvents.movieOn(indNeg(end)+1:end,:);
    end
    
    indLast = 1;
    for i = 1:size(behEvents.movieOn,1)
        indTmp = find(ArdSyncMsec(indLast:end,1) >= behEvents.movieOn(i,1));
        if(~isempty(indTmp))
            indLast = indTmp(1)+indLast-1;
            if(TDTSyncMsecTmp(indLast) == 0)
                behEventsTdt.movieOn(i,1) = 0;
                continue;
            end
            if(indLast ~= 1) % convert the arduino time to TDT time
                behEventsTdt.movieOn(i,1) = TDTSyncMsecTmp(indLast-1)...
                    + (TDTSyncMsecTmp(indLast) - TDTSyncMsecTmp(indLast-1))...
                    /(ArdSyncMsec(indLast,1) - ArdSyncMsec(indLast-1,1))...
                    *(behEvents.movieOn(i,1)-ArdSyncMsec(indLast-1,1));
                if(behEventsTdt.movieOn(i,1) < 0)
                    behEventsTdt.movieOn(i,1) = 0;
                end
            else
                behEventsTdt.movieOn(i,1) = TDTSyncMsecTmp(indLast) ...
                    - (TDTSyncMsecTmp(indLast+1) - TDTSyncMsecTmp(indLast))...
                    /(ArdSyncMsec(indLast+1,1) - ArdSyncMsec(indLast,1))...
                    *(ArdSyncMsec(indLast,1) - behEvents.movieOn(i,1));
            end
        else
            behEventsTdt.movieOn(i,1) = TDTSyncMsecTmp(end) ...
                + (TDTSyncMsecTmp(end) - TDTSyncMsecTmp(end-1))...
                /(ArdSyncMsec(end,1) - ArdSyncMsec(end-1,1))...
                *(behEvents.movieOn(i,1)-ArdSyncMsec(end,1));
        end
    end
    behEventsTdt.movieOn(:,2) = resamp(behEventsTdt.movieOn(:,1),lfpFreq);
    behEventsTdt.movieOn(:,1) = resamp(behEventsTdt.movieOn(:,1),sampleFreq);
end

%%%%%% added by Yingxue on 03/30/3019
% behEvents.movieId
if(isfield(behEvents,'movieID'))
    behEventsTdt.movieID = behEvents.movieID;
end
%%%%%% added by Yingxue on 03/30/3019

%%
% behEvents.TDTMovieOn
% first column: time stamps sampled at sampleFreq
% second column: time stamps sampled at lfpFreq
if(isfield(behEvents,'TDTMovieOn'))
    % check whether there is time reverse, if so, delete the events
    % before the time reverse
    diffTm = diff(behEvents.TDTMovieOn(:,1));
    indNeg = find(diffTm < 0);
    if(~isempty(indNeg))
        behEvents.TDTMovieOn = behEvents.TDTMovieOn(indNeg(end)+1:end,:);
    end
    
    behEventsTdt.TDTMovieOn(:,2) = resamp(behEvents.TDTMovieOn(:,1),lfpFreq);
    behEventsTdt.TDTMovieOn(:,1) = resamp(behEvents.TDTMovieOn(:,1),sampleFreq);
end

% behEvents.TDTMovieOff
% first column: time stamps sampled at sampleFreq
% second column: time stamps sampled at lfpFreq
if(isfield(behEvents,'TDTMovieOff'))
    % check whether there is time reverse, if so, delete the events
    % before the time reverse
    diffTm = diff(behEvents.TDTMovieOff(:,1));
    indNeg = find(diffTm < 0);
    if(~isempty(indNeg))
        behEvents.TDTMovieOff = behEvents.TDTMovieOff(indNeg(end)+1:end,:);
    end
    
    behEventsTdt.TDTMovieOff(:,2) = resamp(behEvents.TDTMovieOff(:,1),lfpFreq);
    behEventsTdt.TDTMovieOff(:,1) = resamp(behEvents.TDTMovieOff(:,1),sampleFreq);
end

%%
%%%%%% added by Yingxue on 03/30/3019
% behEvents.TDTMovieOnL
% first column: time stamps sampled at sampleFreq
% second column: time stamps sampled at lfpFreq
if(isfield(behEvents,'TDTMovieOnL'))
    % check whether there is time reverse, if so, delete the events
    % before the time reverse
    diffTm = diff(behEvents.TDTMovieOnL(:,1));
    indNeg = find(diffTm < 0);
    if(~isempty(indNeg))
        behEvents.TDTMovieOnL = behEvents.TDTMovieOnL(indNeg(end)+1:end,:);
    end
    
    behEventsTdt.TDTMovieOnL(:,2) = resamp(behEvents.TDTMovieOnL(:,1),lfpFreq);
    behEventsTdt.TDTMovieOnL(:,1) = resamp(behEvents.TDTMovieOnL(:,1),sampleFreq);
end

% behEvents.TDTMovieOffL
% first column: time stamps sampled at sampleFreq
% second column: time stamps sampled at lfpFreq
if(isfield(behEvents,'TDTMovieOffL'))
    % check whether there is time reverse, if so, delete the events
    % before the time reverse
    diffTm = diff(behEvents.TDTMovieOffL(:,1));
    indNeg = find(diffTm < 0);
    if(~isempty(indNeg))
        behEvents.TDTMovieOffL = behEvents.TDTMovieOffL(indNeg(end)+1:end,:);
    end
    
    behEventsTdt.TDTMovieOffL(:,2) = resamp(behEvents.TDTMovieOffL(:,1),lfpFreq);
    behEventsTdt.TDTMovieOffL(:,1) = resamp(behEvents.TDTMovieOffL(:,1),sampleFreq);
end

% behEvents.TDTMovieOnR
% first column: time stamps sampled at sampleFreq
% second column: time stamps sampled at lfpFreq
if(isfield(behEvents,'TDTMovieOnR'))
    % check whether there is time reverse, if so, delete the events
    % before the time reverse
    diffTm = diff(behEvents.TDTMovieOnR(:,1));
    indNeg = find(diffTm < 0);
    if(~isempty(indNeg))
        behEvents.TDTMovieOnR = behEvents.TDTMovieOnR(indNeg(end)+1:end,:);
    end
    
    behEventsTdt.TDTMovieOnR(:,2) = resamp(behEvents.TDTMovieOnR(:,1),lfpFreq);
    behEventsTdt.TDTMovieOnR(:,1) = resamp(behEvents.TDTMovieOnR(:,1),sampleFreq);
end

% behEvents.TDTMovieOffR
% first column: time stamps sampled at sampleFreq
% second column: time stamps sampled at lfpFreq
if(isfield(behEvents,'TDTMovieOffR'))
    % check whether there is time reverse, if so, delete the events
    % before the time reverse
    diffTm = diff(behEvents.TDTMovieOffR(:,1));
    indNeg = find(diffTm < 0);
    if(~isempty(indNeg))
        behEvents.TDTMovieOffR = behEvents.TDTMovieOffR(indNeg(end)+1:end,:);
    end
    
    behEventsTdt.TDTMovieOffR(:,2) = resamp(behEvents.TDTMovieOffR(:,1),lfpFreq);
    behEventsTdt.TDTMovieOffR(:,1) = resamp(behEvents.TDTMovieOffR(:,1),sampleFreq);
end
%%%%% added by Yingxue on 03/30/3019
%%

% behEvents.lickPeriod
% first column: time stamps sampled at sampleFreq
% second column: time stamps sampled at lfpFreq
if(isfield(behEvents,'lickPeriod'))
    % check whether there is time reverse, if so, delete the events
    % before the time reverse
    diffTm = diff(behEvents.lickPeriod(:,1));
    indNeg = find(diffTm < 0);
    if(~isempty(indNeg))
        behEvents.lickPeriod = behEvents.lickPeriod(indNeg(end)+1:end,:);
    end
    
    indLast = 1;
    for i = 1:size(behEvents.lickPeriod,1)
        indTmp = find(ArdSyncMsec(indLast:end,1) >= behEvents.lickPeriod(i,1));
        if(~isempty(indTmp))
            indLast = indTmp(1)+indLast-1;
            if(TDTSyncMsecTmp(indLast) == 0)
                behEventsTdt.lickPeriod(i,1) = 0;
                continue;
            end
            if(indLast ~= 1) % convert the arduino time to TDT time
                behEventsTdt.lickPeriod(i,1) = TDTSyncMsecTmp(indLast-1)...
                    + (TDTSyncMsecTmp(indLast) - TDTSyncMsecTmp(indLast-1))...
                    /(ArdSyncMsec(indLast,1) - ArdSyncMsec(indLast-1,1))...
                    *(behEvents.lickPeriod(i,1)-ArdSyncMsec(indLast-1,1));
                if(behEventsTdt.lickPeriod(i,1) < 0)
                    behEventsTdt.lickPeriod(i,1) = 0;
                end
            else
                behEventsTdt.lickPeriod(i,1) = TDTSyncMsecTmp(indLast) ...
                    - (TDTSyncMsecTmp(indLast+1) - TDTSyncMsecTmp(indLast))...
                    /(ArdSyncMsec(indLast+1,1) - ArdSyncMsec(indLast,1))...
                    *(ArdSyncMsec(indLast,1) - behEvents.lick(i,1));
            end
        else
            behEventsTdt.lickPeriod(i,1) = TDTSyncMsecTmp(end) ...
                + (TDTSyncMsecTmp(end) - TDTSyncMsecTmp(end-1))...
                /(ArdSyncMsec(end,1) - ArdSyncMsec(end-1,1))...
                *(behEvents.lickPeriod(i,1)-ArdSyncMsec(end,1));
        end
    end
    behEventsTdt.lickPeriod(:,2) = resamp(behEventsTdt.lickPeriod(:,1),lfpFreq);
    behEventsTdt.lickPeriod(:,1) = resamp(behEventsTdt.lickPeriod(:,1),sampleFreq);
    behEventsTdt.lickPeriod(:,3) = behEvents.lickPeriod(:,2);
end

%%
%%%%%% added by Yingxue on 03/30/3019
% left vs right trial
if(isfield(behEvents,'trStartCueLR'))
    behEventsTdt.trStartCueLR = behEvents.trStartCueLR;
    behEventsTdt.startCueEndDist = behEvents.startCueEndDist;
else
    behEventsTdt.trStartCueLR = [];
end

%% analyze maze type and maze sessions
% added on 1/22/2019 by Yingxue Wang
% changed by Yingxue on 3/30/2019
behEventsTdt.mazeType = zeros(1,length(behEventsTdt.trialDescr(:,1)));
behEventsTdt.mazeSess = zeros(1,length(behEventsTdt.trialDescr(:,1)));
mazeSess = 1;

for i = 1:length(behEventsTdt.trialDescr(:,1))
    if(behEventsTdt.taskDescr(i,4) == 8)
        if(behEventsTdt.movieTDescr{i}(end-1) == 1) % auto water delivery
            if(behEventsTdt.movieTDescr{i}(end) == 0)
                % passive water delivery with start cue and not clamp the belt
                behEventsTdt.mazeType(i) = 9;
            else % passive water delivery with start cue, but the belt is clamped at the end of each trial
                behEventsTdt.mazeType(i) = 8;
            end
        else
            if(behEventsTdt.movieTDescr{i}(end) == 0)
                % active water delivery with start cue and not clamp the belt
                behEventsTdt.mazeType(i) = 10;
            else % active water delivery with start cue, but the belt is clamped at the end of each trial
                behEventsTdt.mazeType(i) = 11;
            end
        end
    elseif(behEventsTdt.taskDescr(i,4) == 9)
        % go/no go task, cue on the left screen -- go; cue on the right screen -- no go
        %%%% added by Yingxue on 04/11/2019
        if(behEventsTdt.startCueEndDist(i) > 180 && ...
                behEventsTdt.startCueEndDist(i) < 210)
            %%% if the start cue ends at 180 cm
            behEventsTdt.mazeType(i) = 15;
        elseif(behEventsTdt.startCueEndDist(i) > 120 && ...
                behEventsTdt.startCueEndDist(i) < 150)
            %%% if the start cue ends at 120 cm
            behEventsTdt.mazeType(i) = 14;
        elseif(behEventsTdt.startCueEndDist(i) > 60 && ...
                behEventsTdt.startCueEndDist(i) < 90)
            %%% if the start cue ends at 60 cm
            behEventsTdt.mazeType(i) = 13;
        elseif(behEventsTdt.startCueEndDist(i) > 0 && ...
                behEventsTdt.startCueEndDist(i) < 30)
            behEventsTdt.mazeType(i) = 12;
        else
            behEventsTdt.mazeType(i) = 0;
        end
    elseif(behEventsTdt.taskDescr(i,4) == 22)
        % double lick port
        behEventsTdt.mazeType(i) = 22;        
    end
    if(behEventsTdt.mazeType(i) == 0)
        behEventsTdt.mazeSess(i) = 0;
    else
        if(i > 1 && behEventsTdt.mazeSess(i-1) ~= 0 && ...
                (behEventsTdt.mazeType(i) ~= behEventsTdt.mazeType(i-1) || ...
                behEventsTdt.taskDescr(i,3) ~= behEventsTdt.taskDescr(i-1,3)))
            mazeSess = mazeSess + 1;
        end
        behEventsTdt.mazeSess(i) = mazeSess;
    end
end


%%
fullPathB = [baseFileName 'BTDT.mat'];
save(fullPathB,'behEventsTdt');
