function Arduino2LFPtime_smTr_opto(fileName,lfpFreq)
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
        behEventsTdt.trialDescr(:,2:3) = behEvents.trialDescr(:,2:3);
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
        numPulses = 1;

        numPC = length(behEvents.stimOn);
        for i = 1:numPC
            indTmp = find(behEvents.pulsePar(:,1) <= behEvents.stimOn(i,1),1,'last');
            if(~isempty(indTmp))
                numPulsePerStim =  behEvents.pulsePar(indTmp,5) * ...
                    behEvents.pulsePar(indTmp,7);
                behEvents.stimOn(i,2) = numPulses; % starting pulse index of the current stimulation
                behEvents.stimOn(i,3) = behEvents.pulsePar(indTmp,5); % number of pulses
                behEvents.stimOn(i,4) = behEvents.pulsePar(indTmp,3)/1000; % pulse width in ms
                behEvents.stimOn(i,5) = behEvents.pulsePar(indTmp,4)/1000; % pulse period in ms
                behEvents.stimOn(i,6) = behEvents.pulsePar(indTmp,6)/1000; % off time between pulse trains in ms
                behEvents.stimOn(i,7) = behEvents.pulsePar(indTmp,7); % repeats
            else
                disp('no pulse parameter detected')
            end
            numPulses = numPulses+numPulsePerStim;
        end
        behEventsTdt.stimOn(:,1) = resamp(behEvents.stimOn(:,1),lfpFreq);
        behEventsTdt.stimOn(:,2:7) = behEvents.stimOn(:,2:7);
    end
    save(fullPathB,'behEvents');
    
    %%
    % behEvents.stimOff
    % first column: time stamps sampled at lfpFreq
    if(isfield(behEvents,'stimOff'))
        numPO = length(behEvents.stimOff);
        if(numPO ~= numPC)
            disp('non-equal number of stim on and off pulses');
        end
        behEventsTdt.stimOff(:,1) = resamp(behEvents.stimOff(:,1),lfpFreq);
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

    for i = 1:length(behEventsTdt.trialDescr(:,1))
        if(behEventsTdt.taskDescr(i,3) == 22)
            % active licking control
            behEventsTdt.mazeType(i) = behEventsTdt.movieTDescr{i}(12);
        end
    end
    
    behEventsTdt.mazeTypeMod = behEventsTdt.mazeType;
    
    % count the number of control trials
    indZ = behEventsTdt.mazeType == 0;
    f = find(diff([0 indZ 0] == 1));
    p = f(1:2:end-1); % start indices
    y = f(2:2:end) - p; % consecutive ones' count
    x = find(y < 5);
    for i = 1:length(x)
        if(p(x(i)) > 1)
            behEventsTdt.mazeTypeMod(p(x(i)):p(x(i))+y(x(i))-1) = behEventsTdt.mazeTypeMod(p(x(i))-1);
        else
            behEventsTdt.mazeTypeMod(p(x(i)):p(x(i))+y(x(i))-1) = behEventsTdt.mazeTypeMod(p(x(i+1))-1);
        end
    end
    mazeSess = 1;
    indZ = find(diff([behEventsTdt.mazeTypeMod(1) behEventsTdt.mazeTypeMod]) ~= 0);
    if(~isempty(indZ))
        for i = 1:length(indZ)
            if(i == 1)
                behEventsTdt.mazeSess(1:indZ(1)-1) = mazeSess;
            else
                behEventsTdt.mazeSess(indZ(i-1):indZ(i)-1) = mazeSess;
            end
            mazeSess = mazeSess + 1;
        end
        behEventsTdt.mazeSess(indZ(end):end) = mazeSess;
    else
        behEventsTdt.mazeSess = ones(1,length(behEventsTdt.mazeTypeMod));
    end
    
    %%
    fullPathB = [fileName 'BLFP.mat'];
    save(fullPathB,'behEventsTdt');
end
