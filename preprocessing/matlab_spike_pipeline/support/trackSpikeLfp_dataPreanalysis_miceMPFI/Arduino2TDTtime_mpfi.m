function Arduino2TDTtime_mpfi(baseFileName,sampleFreq,lfpFreq)
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
        behEventsTdt.trialDescr(:,3:6) = behEvents.trialDescr(:,2:5);
    end
    
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
    
    % behEvents.beam
    % first column: time stamps sampled at sampleFreq
    % second column: time stamps sampled at lfpFreq
    if(isfield(behEvents,'beam'))
        % check whether there is time reverse, if so, delete the events
        % before the time reverse
        diffTm = diff(behEvents.beam(:,1));
        indNeg = find(diffTm < 0);
        if(~isempty(indNeg))
            behEvents.beam = behEvents.beam(indNeg(end)+1:end,:);
        end
        
        indLast = 1;
        for i = 1:size(behEvents.beam,1)
           indTmp = find(ArdSyncMsec(indLast:end,1) >= behEvents.beam(i,1));
           if(~isempty(indTmp))
               indLast = indTmp(1)+indLast-1;
               if(TDTSyncMsecTmp(indLast) == 0)
                   behEventsTdt.beam(i,1) = 0;
                   continue;
               end
               if(indLast ~= 1) % convert the arduino time to TDT time
                    behEventsTdt.beam(i,1) = TDTSyncMsecTmp(indLast-1) ...
                        +(TDTSyncMsecTmp(indLast) - TDTSyncMsecTmp(indLast-1))...
                        /(ArdSyncMsec(indLast,1) - ArdSyncMsec(indLast-1,1))...
                        *(behEvents.beam(i,1)-ArdSyncMsec(indLast-1,1));
                    if(behEventsTdt.beam(i,1) < 0)
                        behEventsTdt.beam(i,1) = 0;
                    end
               else
                   behEventsTdt.beam(i,1) = TDTSyncMsecTmp(indLast) ...
                        - (TDTSyncMsecTmp(indLast+1) - TDTSyncMsecTmp(indLast))...
                        /(ArdSyncMsec(indLast+1,1) - ArdSyncMsec(indLast,1))...
                        *(ArdSyncMsec(indLast,1) - behEvents.beam(i,1));
               end
           else
               behEventsTdt.beam(i,1) = TDTSyncMsecTmp(end) ...
                   + (TDTSyncMsecTmp(end) - TDTSyncMsecTmp(end-1))...
                   /(ArdSyncMsec(end,1) - ArdSyncMsec(end-1,1))...
                   *(behEvents.beam(i,1)-ArdSyncMsec(end,1));
           end
        end
        behEventsTdt.beam(:,2) = resamp(behEventsTdt.beam(:,1),lfpFreq);
        behEventsTdt.beam(:,1) = resamp(behEventsTdt.beam(:,1),sampleFreq);
        behEventsTdt.beam(:,3:4) = behEvents.beam(:,2:3);
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
    
    fullPathB = [baseFileName 'BTDT.mat'];
    save(fullPathB,'behEventsTdt');
    
    % check whether the beam breaks are triggered in the correct order 
    % (meaning 1 2 3 instead of 1 3 2)
    beamOrder=input('Are the beam breaks triggered in the correct order (1 2 3)? Y/N [Y]: ','s');
    if(isempty(beamOrder))
        beamOrder = 'N';
    end
    if(~contains(beamOrder,'Y') && ~contains(beamOrder,'y'))
       swapBeam2and3InBTDT(baseFileName); 
    end
