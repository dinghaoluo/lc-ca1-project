function Arduino2LFPtime_smTr(fileName,lfpFreq)
% convert the arduino time to lfp time
    
    %%%%%%%%% check arguments
    if nargin<2
        disp('At least two arguments are needed for this function.');
        return;
    end

    %%%%%%%%% load recording file
    fullPathB = [fileName 'B.mat'];

    if(exist(fullPathB))
        load(fullPathB);
    else
        disp('Behavioral event file does not exist.')
        return;
    end

    behEventsTdt = [];

    %%
    % behEvents.taskDescr
    % first column: time stamps sampled at lfpFreq
    if(isfield(behEvents,'taskDescr'))
        
        % check whether there is time reverse
        diffTm = diff(behEvents.taskDescr(:,1));
        indNeg = find(diffTm < 0);
        if(~isempty(indNeg))
            behEvents.taskDescr(1:indNeg,1) = 0;
        end
        
        behEventsTdt.taskDescr(:,1) = resamp(behEvents.taskDescr(:,1),lfpFreq);
        behEventsTdt.taskDescr(:,2:3) = behEvents.taskDescr(:,2:3);
    end
    
    % behEvents.movieTDescr
    if(isfield(behEvents,'movieTDescr'))
        behEventsTdt.movieTDescr = behEvents.movieTDescr;
    end
    
    % behEvents.trialDescr 
    % first column: time stamps sampled at lfpFreq
    if(isfield(behEvents,'trialDescr'))
        behEventsTdt.trialDescr(:,1) = resamp(behEvents.trialDescr(:,1),lfpFreq);
        behEventsTdt.trialDescr(:,2) = behEvents.trialDescr(:,2);
    end
    
    %%
    % behEvents.lick
    % first column: time stamps sampled at lfpFreq
    if(isfield(behEvents,'lick'))
        behEventsTdt.lick(:,1) = resamp(behEvents.lick(:,1),lfpFreq);
        behEventsTdt.lick(:,2:3) = behEvents.lick(:,2:3);
    else
        behEventsTdt.lick = [];
    end
    
    %%
    % behEvents.stimOn
    % first column: time stamps sampled at lfpFreq
    if(isfield(behEvents,'stimOn'))
        behEventsTdt.stimOn(:,1) = resamp(behEvents.stimOn(:,1),lfpFreq);
%         behEventsTdt.stimOn(:,2:4) = behEvents.stimOn(:,2:4);
    end
    
    %%
    % behEvents.wheel
    % first column: time stamps sampled at lfpFreq
    if(isfield(behEvents,'wheel'))
        behEventsTdt.wheel(:,1) = resamp(behEvents.wheel(:,1),lfpFreq);
        behEventsTdt.wheel(:,2:3) = behEvents.wheel(:,2:3);
    end
    
    %%
    % behEvents.pump
    % first column: time stamps sampled at lfpFreq
    if(isfield(behEvents,'pump'))
        behEventsTdt.pump(:,1) = resamp(behEvents.pump(:,1),lfpFreq);
        behEventsTdt.pump(:,2:3) = behEvents.pump(:,2:3);
    else
        behEventsTdt.pump = [];
    end
    
    %%
    % behEvents.movieOn
    % first column: time stamps sampled at lfpFreq
    if(isfield(behEvents,'movieOn'))
        behEventsTdt.movieOn(:,1) = resamp(behEvents.movieOn(:,1),lfpFreq);
        behEventsTdt.movieOn(:,2) = behEvents.movieOn(:,2);
    end
    
    %% analyze maze type and maze sessions
    % added on 1/22/2019 by Yingxue Wang
    % changed by Yingxue on 3/30/2019
    behEventsTdt.mazeType = zeros(1,length(behEvents.trialDescr(:,1)));
    behEventsTdt.mazeSess = zeros(1,length(behEventsTdt.trialDescr(:,1)));
    mazeSess = 1;

    for i = 1:length(behEventsTdt.trialDescr(:,1))
        if(behEventsTdt.taskDescr(i,3) == 8)
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
        elseif(behEventsTdt.taskDescr(i,3) == 9)
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
    fullPathB = [fileName 'BLFP.mat'];
    save(fullPathB,'behEventsTdt');
end
