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
        behEventsTdt.trialDescr(:,3) = behEvents.trialDescr(:,2);
    end
    
    %% added on 1/21/2019 
    % trial start time is defined by taskDescr, and end time is defined by trialDescr
    if(length(behEventsTdt.taskDescr(:,1)) == length(behEventsTdt.trialDescr(:,1)))
        if(sum(behEventsTdt.trialDescr(:,1)-behEventsTdt.taskDescr(:,1) <= 0) > 0)
            disp('Error: trial end time <= trial start time');
            return;
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
    
    %% analyze maze type and maze sessions
    % added on 1/22/2019 by Yingxue Wang
    behEventsTdt.mazeType = zeros(1,length(behEventsTdt.trialDescr(:,1)));
    behEventsTdt.mazeSess = zeros(1,length(behEventsTdt.trialDescr(:,1)));
    mazeSess = 1;
    
    % label the free running session (the first session in most cases)
    indTrFirstMovie = ...
        find(behEvents.taskDescr(:,1) < behEvents.movieOn(1,1),1,'last');
    if(~isempty(indTrFirstMovie) && indTrFirstMovie > 1)
        behEventsTdt.mazeType(1:indTrFirstMovie-1) = 1;
        behEventsTdt.mazeSess(1:indTrFirstMovie-1) = mazeSess;
        mazeSess = mazeSess + 1;
    else
        indTrFirstMovie = 0;
    end
    
    % label the start cue session (likely the second session)
    % in both the start cue and the free running sessions, the second 
    % parameter in the taskDescr is 0, but the start cue session has MV
    % command will the free running does not
    indTrStartCues = ...
        find(behEvents.taskDescr(:,3) == 0);
    indTrStartCues = indTrStartCues(indTrStartCues > indTrFirstMovie);
    if(~isempty(indTrStartCues))
        diffTr = diff(indTrStartCues);
        sessGap = find(diffTr > 1);
        if(indTrStartCues(1) == indTrFirstMovie+1)
            if(~isempty(sessGap))
                indTrStartCuesLeft = indTrStartCues(sessGap(1)+1:end);
                indTrStartCues = indTrStartCues(1:sessGap(1));
            else
                indTrStartCuesLeft = [];
            end
            behEventsTdt.mazeType(indTrStartCues(1:end-1))= 2;
            behEventsTdt.mazeSess(indTrStartCues(1:end-1)) = mazeSess;   
            mazeSess = mazeSess + 1;
        end
    end
    
    % label the three cue session (likely the second or third session)
    indTrThreeCues = ...
        find(behEvents.taskDescr(:,3) > 1);
    if(~isempty(indTrThreeCues)) 
        if(indTrThreeCues(1) == indTrStartCues(end)+1 || ...
                indTrThreeCues(1) == indTrFirstMovie+1)
            behEventsTdt.mazeType(indTrThreeCues(1:end-1))= 3;
            behEventsTdt.mazeSess(indTrThreeCues(1:end-1)) = mazeSess;
            mazeSess = mazeSess + 1;   
        end
    end
    
    % label the free running session after the three cue session if there
    % is one, and label the double reward session followed by the free
    % running session if there is one
    if(~isempty(indTrStartCuesLeft))   
        for tr = indTrStartCuesLeft'
            if(behEventsTdt.movieTDescr{tr}(end)>1) % swith to double reward
                indTrFreeRun = indTrStartCuesLeft(1:tr-1);
                indTrDoubleReward = indTrStartCuesLeft(tr:end);
                break;
            end
        end
        if(tr == indTrStartCuesLeft(end))
            indTrFreeRun = indTrStartCuesLeft;
            indTrDoubleReward = [];
        end
        if(~isempty(indTrFreeRun) && ...
                indTrFreeRun(1) == indTrThreeCues(end)+1)
            if(indTrFreeRun(end) == length(behEventsTdt.mazeType))
                behEventsTdt.mazeType(indTrFreeRun(1:end))= 1;
                behEventsTdt.mazeSess(indTrFreeRun(1:end)) = mazeSess;
            else
                behEventsTdt.mazeType(indTrFreeRun(1:end-1))= 1;
                behEventsTdt.mazeSess(indTrFreeRun(1:end-1)) = mazeSess;
            end
            mazeSess = mazeSess + 1;   
        end
        if(~isempty(indTrDoubleReward) && ...
                indTrDoubleReward(1) == indTrFreeRun(end)+1)
            behEventsTdt.mazeType(indTrDoubleReward(1:end))= 4;
            behEventsTdt.mazeSess(indTrDoubleReward(1:end)) = mazeSess;
        end
    end
    
    
    %%
    fullPathB = [baseFileName 'BTDT.mat'];
    save(fullPathB,'behEventsTdt');
