function check = Arduino2TDTtime_smTr_opto2P(baseFileName,sampleFreq,lfpFreq)
% using the sync signal to align the arduino time and the TDT time

check = 0;

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

if(~isfield(behEvents,'TwoPSyncMsec'))
    ArdSyncMsec = [];
else
    ArdSyncMsec = behEvents.TwoPSyncMsec;
end
ardSyncLen = size(ArdSyncMsec,1);
tdtSyncLen = length(behEvents.TDTsyncMsec);
TDTSyncMsecTmp = behEvents.TDTsyncMsec;
if(ardSyncLen > tdtSyncLen)
    if(ardSyncLen ~= tdtSyncLen + 2 && ardSyncLen ~= tdtSyncLen + 3)
        disp('The number of frames recorded in the behavior file is not equal to the number recorded in the imaging + 2 or + 3.');
        ArdSyncMsec = ArdSyncMsec(1:tdtSyncLen);
        check = 1;
    else
        ArdSyncMsec = ArdSyncMsec(1:tdtSyncLen);
    end
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
    taskDescr = behEventsTdt.taskDescr(:,1);
    behEventsTdt.taskDescr(:,2) = resamp2P(taskDescr,lfpFreq);
    behEventsTdt.taskDescr(:,1) = resamp2P(taskDescr,sampleFreq);
    taskDescr = behEvents.taskDescr(:,1) - ArdSyncMsec(1);
    behEventsTdt.taskDescr(:,3) = resamp2P(taskDescr,sampleFreq); % time based on behavioral time
    behEventsTdt.taskDescr(:,4) = resamp2P(taskDescr,lfpFreq);
    behEventsTdt.taskDescr(:,5:6) = behEvents.taskDescr(:,2:3);
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
    trialDescr = behEventsTdt.trialDescr(:,1);
    behEventsTdt.trialDescr(:,1) = resamp2P(trialDescr,sampleFreq);
    behEventsTdt.trialDescr(:,2) = resamp2P(trialDescr,lfpFreq);    
    trialDescr = behEvents.trialDescr(:,1) - ArdSyncMsec(1);
    behEventsTdt.trialDescr(:,3) = resamp2P(trialDescr,sampleFreq); % time based on behavioral time
    behEventsTdt.trialDescr(:,4) = resamp2P(trialDescr,lfpFreq);
    behEventsTdt.trialDescr(:,5) = behEvents.trialDescr(:,3);
    behEventsTdt.trialDescr(:,6) = behEvents.trialDescr(:,2); % performance
end

%% added on 1/21/2019
% trial start time is defined by taskDescr, and end time is defined by trialDescr
if(length(behEventsTdt.taskDescr(:,1)) == length(behEventsTdt.trialDescr(:,1)))
    if(sum(behEventsTdt.trialDescr(:,1)-behEventsTdt.taskDescr(:,1) <= 0) > 0)
        behEventsTdt.trialDescr = behEventsTdt.trialDescr(2:end,:);
        behEventsTdt.trialDescr(:,5) = behEventsTdt.trialDescr(:,5)-1;
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

% % disp('scatter plot...');
% % disp(length(behEventsTdt.taskDescr(:,1)))
% % disp(length(behEventsTdt.taskDescr(:,2)))
% %disp(transpose(behEventsTdt.taskDescr))
% 
% y = [transpose(behEventsTdt.taskDescr(:,1)); transpose(behEventsTdt.taskDescr(:,2))];
% x = [(1:length(behEvents.taskDescr(:,1))); (1:length(behEvents.taskDescr(:,2)))];
% %x = [(1:length(behEvents.taskDescr(:,1))); 2:120];
% 
% %disp((x))
% %disp((y))
% 
% c1 = [1, 0, 0]
% c2 = [0, 1, 0]
% 
% % scatter(x(1,:), y(1,:), [], c1); hold on;
% % scatter(x(2,:), y(2,:), [], c2); 

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
    lick = behEventsTdt.lick(:,1);
    behEventsTdt.lick(:,1) = resamp2P(lick,sampleFreq);
    behEventsTdt.lick(:,2) = resamp2P(lick,lfpFreq);    
    lick = behEvents.lick(:,1) - ArdSyncMsec(1);
    behEventsTdt.lick(:,3) = resamp2P(lick,sampleFreq); % time based on behavioral time
    behEventsTdt.lick(:,4) = resamp2P(lick,lfpFreq);
    behEventsTdt.lick(:,5:6) = behEvents.lick(:,2:3);
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
            behEventsTdt.stimPulse(indStimP+n,5) = i; % stimulation no.
            behEventsTdt.stimPulse(indStimP+n,6) = behEvents.stimOnPar(i,2); % pulse width
            behEventsTdt.stimPulse(indStimP+n,7) = n; % num within stim

            behEventsTdt.stimPulse(indStimP+n,1) = behEvents.TDTstimMsec(indStimP+n);
            behEventsTdt.stimPulse(indStimP+n,8) = behEvents.TDTpulseWidthMsec(indStimP+n);
               
            behEventsTdt.stimOnDiode{indStimP+n} = behEvents.stimOnDiode{i};
            
            if(behEvents.stimOn(i,1) == -1 && n == 1)
                behEventsTdt.stimOn(i,1) = behEventsTdt.stimPulse(indStimP+n,1);
            end
            indStim = indStim + 1;
        end
        
        indStimP = indStimP + behEvents.stimOn(i,4);
    end
    stimOnTmp = behEventsTdt.stimOn(:,1);
    behEventsTdt.stimOn(:,1) = resamp2P(stimOnTmp,sampleFreq);
    behEventsTdt.stimOn(:,2) = resamp2P(stimOnTmp,lfpFreq);    
    stimOnTmp = behEvents.stimOn(:,1) - ArdSyncMsec(1);
    behEventsTdt.stimOn(:,3) = resamp2P(stimOnTmp,sampleFreq); % time based on behavioral time
    behEventsTdt.stimOn(:,4) = resamp2P(stimOnTmp,lfpFreq);
    behEventsTdt.stimOn(:,5:6) = behEvents.stimOn(:,3:4);  
    % changed on 2/17/2010 because there is an error in getStim_mpfi_smTr line 21 (no -1 in column 5)
    indTag = behEvents.stimOn(:,1) < 0;
    stimIsTag = ones(length(behEvents.stimOn(:,1)),1);
    stimIsTag(indTag) = -1;
    behEventsTdt.stimOn(:,7) = stimIsTag;  
    behEventsTdt.stimOn(:,8:11) = behEvents.stimOnPar(:,2:end);
    %
    stimPulseTmp = behEventsTdt.stimPulse(:,1);
    behEventsTdt.stimPulse(:,1) = resamp2P(stimPulseTmp,sampleFreq);
    behEventsTdt.stimPulse(:,2) = resamp2P(stimPulseTmp,lfpFreq);
    stimPulseTmp = behEvents.stimPulse(:,1) - ArdSyncMsec(1);
    behEventsTdt.stimPulse(:,3) = resamp2P(stimPulseTmp,sampleFreq); % time based on behavioral time
    behEventsTdt.stimPulse(:,4) = resamp2P(stimPulseTmp,lfpFreq);
    behEventsTdt.stimPulse(:,8) = resamp2P(behEventsTdt.stimPulse(:,8),sampleFreq);
    behEventsTdt.stimPulse(:,9) = resamp2P(behEventsTdt.stimPulse(:,8),lfpFreq);
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
    trialT = behEventsTdt.trialT(1);
    behEventsTdt.trialT(1) = resamp2P(trialT,sampleFreq);
    behEventsTdt.trialT(2) = resamp2P(trialT,lfpFreq);
    trialT = behEvents.trialT(1) - ArdSyncMsec(1);
    behEventsTdt.trialT(:,3) = resamp2P(trialT,sampleFreq); % time based on behavioral time
    behEventsTdt.trialT(:,4) = resamp2P(trialT,lfpFreq);
    behEventsTdt.trialT(5) = behEvents.trialT(2);
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
                
%                    disp((TDTSyncMsecTmp(indLast) - TDTSyncMsecTmp(indLast-1))...
%                     /(ArdSyncMsec(indLast,1) - ArdSyncMsec(indLast-1,1))...
%                     *(behEvents.wheel(i,1)-ArdSyncMsec(indLast-1,1)))
                if(behEventsTdt.wheel(i,1) < 0)
                    behEventsTdt.wheel(i,1) = 0;
                end
            else
                behEventsTdt.wheel(i,1) = (TDTSyncMsecTmp(indLast) ...
                    - (TDTSyncMsecTmp(indLast+1) - TDTSyncMsecTmp(indLast))...
                    /(ArdSyncMsec(indLast+1,1) - ArdSyncMsec(indLast,1))...
                    *(ArdSyncMsec(indLast,1) - behEvents.wheel(i,1))) ;
            end
        else
            behEventsTdt.wheel(i,1) = (TDTSyncMsecTmp(end) ...
                + (TDTSyncMsecTmp(end) - TDTSyncMsecTmp(end-1))...
                /(ArdSyncMsec(end,1) - ArdSyncMsec(end-1,1))...
                *(behEvents.wheel(i,1)-ArdSyncMsec(end,1)));
        end

    end
    wheelTmp = behEventsTdt.wheel(:,1);
    behEventsTdt.wheel(:,1) = resamp2P(wheelTmp,sampleFreq);
    behEventsTdt.wheel(:,2) = resamp2P(wheelTmp,lfpFreq);
    wheelTmp = behEvents.wheel(:,1) - ArdSyncMsec(1);
    behEventsTdt.wheel(:,3) = resamp2P(wheelTmp,sampleFreq); % time based on behavioral time
    behEventsTdt.wheel(:,4) = resamp2P(wheelTmp,lfpFreq);
    behEventsTdt.wheel(:,5:6) = behEvents.wheel(:,2:3);
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
    pumpTmp = behEventsTdt.pump(:,1);
    behEventsTdt.pump(:,1) = resamp2P(pumpTmp,sampleFreq);
    behEventsTdt.pump(:,2) = resamp2P(pumpTmp,lfpFreq);
    pumpTmp = behEvents.pump(:,1) - ArdSyncMsec(1);
    behEventsTdt.pump(:,3) = resamp2P(pumpTmp,sampleFreq); % time based on behavioral time
    behEventsTdt.pump(:,4) = resamp2P(pumpTmp,lfpFreq);
    behEventsTdt.pump(:,5:6) = behEvents.pump(:,2:3);
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
    airpuffTmp = behEventsTdt.airpuff(:,1);
    behEventsTdt.airpuff(:,1) = resamp2P(airpuffTmp,sampleFreq);
    behEventsTdt.airpuff(:,2) = resamp2P(airpuffTmp,lfpFreq);
    airpuffTmp = behEvents.airpuff(:,1) - ArdSyncMsec(1);
    behEventsTdt.airpuff(:,3) = resamp2P(airpuffTmp,sampleFreq); % time based on behavioral time
    behEventsTdt.airpuff(:,4) = resamp2P(airpuffTmp,lfpFreq);
    behEventsTdt.airpuff(:,5) = behEvents.airpuff(:,2);
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
    movieOnTmp = behEventsTdt.movieOn(:,1);
    behEventsTdt.movieOn(:,1) = resamp2P(movieOnTmp,sampleFreq);
    behEventsTdt.movieOn(:,2) = resamp2P(movieOnTmp,lfpFreq);
    movieOnTmp = behEvents.movieOn(:,1) - ArdSyncMsec(1);
    behEventsTdt.movieOn(:,3) = resamp2P(movieOnTmp,sampleFreq); % time based on behavioral time
    behEventsTdt.movieOn(:,4) = resamp2P(movieOnTmp,lfpFreq);
end

%% analyze maze type and maze sessions
% added on 1/22/2019 by Yingxue Wang
% changed by Yingxue on 3/30/2019
behEventsTdt.mazeType = ones(1,length(behEventsTdt.trialDescr(:,1)));
behEventsTdt.mazeSess = zeros(1,length(behEventsTdt.trialDescr(:,1)));
mazeSess = 1;

for i = 1:length(behEventsTdt.trialDescr(:,1))
%     if(behEventsTdt.taskDescr(i,6) == 8)
%         if(behEventsTdt.movieTDescr{i}(end-1) == 1) % auto water delivery
%             if(behEventsTdt.movieTDescr{i}(end) == 0)
%                 % passive water delivery with start cue and not clamp the belt
%                 behEventsTdt.mazeType(i) = 9;
%             else % passive water delivery with start cue, but the belt is clamped at the end of each trial
%                 behEventsTdt.mazeType(i) = 8;
%             end
%         else
%             if(behEventsTdt.movieTDescr{i}(end) == 0)
%                 % active water delivery with start cue and not clamp the belt
%                 behEventsTdt.mazeType(i) = 10;
%             else % active water delivery with start cue, but the belt is clamped at the end of each trial
%                 behEventsTdt.mazeType(i) = 11;
%             end
%         end
%     elseif(behEventsTdt.taskDescr(i,6) == 22)
%         % double lick port
%         behEventsTdt.mazeType(i) = 22;        
%     end
    if(behEventsTdt.mazeType(i) == 0)
        behEventsTdt.mazeSess(i) = 0;
    else
        if(i > 1 && behEventsTdt.mazeSess(i-1) ~= 0 && ...
                (behEventsTdt.mazeType(i) ~= behEventsTdt.mazeType(i-1) || ...
                behEventsTdt.taskDescr(i,5) ~= behEventsTdt.taskDescr(i-1,5)))
            mazeSess = mazeSess + 1;
        end
        behEventsTdt.mazeSess(i) = mazeSess;
    end
end

behEventsTdt.TDTsyncMsec = behEvents.TDTsyncMsec;
behEventsTdt.TDTsyncMsecLfp = round(behEvents.TDTsyncMsec/1000*lfpFreq);

behEventsTdt.ArdSyncMsec = ArdSyncMsec - ArdSyncMsec(1);
behEventsTdt.ArdSyncMsecLfp = round((ArdSyncMsec-ArdSyncMsec(1))/1000*lfpFreq);
behEventsTdt.ArdSyncMsecSamp = round((ArdSyncMsec-ArdSyncMsec(1))/1000*sampleFreq);

%%  removes all events that occur after the final FM frame. this
%   means all behavior events are between the first and last FM
%   frames
cutoffIdx = find(behEventsTdt.wheel(:,4) > behEventsTdt.ArdSyncMsecLfp(end));
if(~isempty(cutoffIdx))
    behEventsTdt.wheel(cutoffIdx(1):end, :) = [];
end
cutoffIdx = find(behEventsTdt.lick(:,4) > behEventsTdt.ArdSyncMsecLfp(end));
if(~isempty(cutoffIdx))
    behEventsTdt.lick(cutoffIdx(1):end, :) = [];
end
cutoffIdx = find(behEventsTdt.movieOn(:,4) > behEventsTdt.ArdSyncMsecLfp(end));
if(~isempty(cutoffIdx))
    behEventsTdt.movieOn(cutoffIdx(1):end, :) = [];
end
cutoffIdx = find(behEventsTdt.pump(:,4) > behEventsTdt.ArdSyncMsecLfp(end));
if(~isempty(cutoffIdx))
    behEventsTdt.pump(cutoffIdx(1):end, :) = [];
end
if(isfield(behEvents,'airpuff'))
    cutoffIdx = find(behEventsTdt.airpuff(:,4) > behEventsTdt.ArdSyncMsecLfp(end));
    if(~isempty(cutoffIdx))
        behEventsTdt.airpuff(cutoffIdx(1):end, :) = [];
    end
end
cutoffIdx = find(behEventsTdt.trialDescr(:,4) > behEventsTdt.ArdSyncMsecLfp(end));
if(~isempty(cutoffIdx))
    behEventsTdt.trialDescr(cutoffIdx(1):end, :) = [];
end
cutoffIdx = find(behEventsTdt.taskDescr(:,4) > behEventsTdt.ArdSyncMsecLfp(end));
if(~isempty(cutoffIdx))
    behEventsTdt.taskDescr(cutoffIdx(1):end, :) = [];
end
behEventsTdt.movieTDescr = behEventsTdt.movieTDescr(1:size(behEventsTdt.taskDescr,1));
behEventsTdt.mazeSess = behEventsTdt.mazeSess(1:size(behEventsTdt.trialDescr,1));
behEventsTdt.mazeType = behEventsTdt.mazeType(1:size(behEventsTdt.trialDescr,1));

%%
fullPathB = [baseFileName 'BTDT.mat'];
save(fullPathB,'behEventsTdt');
