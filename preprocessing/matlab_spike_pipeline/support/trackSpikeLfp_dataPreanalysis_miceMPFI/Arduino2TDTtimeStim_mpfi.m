function Arduino2TDTtimeStim_mpfi(baseFileName,sampleFreq,lfpFreq)
% using the sync signal to align the arduino time and the TDT time for
% stimulation pulses

    %%%%%%%%% check arguments
    if nargin<3
        disp('At least three arguments are needed for this function.');
        return;
    end
    
    %%%%%%%%% load recording file
    fullPathP = [baseFileName 'P.mat'];
        
    if(exist(fullPathP,'file'))
        load(fullPathP);
    else
        disp('stimulation event file does not exist.')
        return;
    end
    
    %%%%%%%% convert the arduino time stamps to the TDT time stamps
    stimEventsTdt = [];
    
    %% check whether there is any stimulation by mistake
    totPulse = 1;
    totCorrPulse = 0; 
    for n = 1:size(stimEvents.pulse,1)
        if(n < size(stimEvents.pulse,1))
            indStimTmp = find(stimEvents.stim(:,1) > stimEvents.pulse(n,1) ...
                & stimEvents.stim(:,1) < stimEvents.pulse(n+1,1));
        else
            indStimTmp = find(stimEvents.stim(:,1) > stimEvents.pulse(n,1));
        end
        totCorrPulse = totCorrPulse + length(indStimTmp)...
                        *stimEvents.pulse(n,4)*stimEvents.pulse(n,6);
    end
    indStart = 1;
    if(totCorrPulse < length(stimEvents.stimInd))
        % only consider the simplest case, where the error happens before
        % the pulse definition
        indWrongStim = find(stimEvents.stim(:,1) < stimEvents.pulse(1,1));
        if(~isempty(indWrongStim))
            totPulse = length(stimEvents.stimInd) - totCorrPulse + 1;
            indStart = length(indWrongStim)+1;
        end
    end
    
    for i = indStart:size(stimEvents.stim,1)
        %% extract information about which diode is activated, here assume
        % only one diode is activated at a time 
        indDiode = find(stimEvents.diode(:,1) < stimEvents.stim(i,1),1,'last');
        diodeAct = find(stimEvents.diode(indDiode,2:7) == 1, 1);
        if(isempty(diodeAct))
            diodeAct = 0;
        end
        
        %% extract information about the stimulation pulse train
        indPulse = find(stimEvents.pulse(:,1) < stimEvents.stim(i,1),1,'last');
        
        numPulse = stimEvents.pulse(indPulse,4);
        pulseWidth = (stimEvents.pulse(indPulse,2)...
                        *stimEvents.pulse(indPulse,3)/100)/1000; % ms
        if(numPulse == 1)
            pulsePeriod = (stimEvents.pulse(indPulse,2)...
                        +stimEvents.pulse(indPulse,5))/1000; % ms
        else
            pulsePeriod = stimEvents.pulse(indPulse,2)/1000; % ms
        end
        
        m = 1;
        while m <= stimEvents.pulse(indPulse,6)
            if(m > 1)
                trainInterv = (stimEvents.stimMsec(totPulse) ...
                                - stimEvents.stimMsec(totPulse-1));
                if(trainInterv > stimEvents.pulse(indPulse,5)/1000*1.05 ...
                        | trainInterv < stimEvents.pulse(indPulse,5)/1000*0.95)
                    break;
                end
            end
                
            if(totPulse+numPulse-1 <= length(stimEvents.stimMsec))
                truePulsePeriod = ...
                    diff(stimEvents.stimMsec(totPulse:totPulse+numPulse-1));
                errPulse = find(truePulsePeriod > pulsePeriod*1.05 ...
                                | truePulsePeriod < pulsePeriod*0.95);
            else
                truePulsePeriod = diff(stimEvents.stimMsec(totPulse:end));
                errPulse = find(truePulsePeriod > pulsePeriod*1.05 ...
                                | truePulsePeriod < pulsePeriod*0.95);
                if(~isempty(errPulse))
                    errPulse = min(errPulse(1),length(truePulsePeriod)+1);
                else
                    errPulse = length(truePulsePeriod)+1;
                end
            end
        
            if(~isempty(errPulse))
                nPulses = errPulse(1)-1; 
                % if the number of valid pulses is < 80% of numPulses, 
                % then consider whether there were pulses from other
                % train being included into the current pulse train
                errPulse = find(truePulsePeriod > pulsePeriod*5);
                if(~isempty(errPulse))
                    totPulse = totPulse + errPulse(end);
                    continue;
                end
            else
                nPulses = numPulse;
            end
            for j = 1:nPulses
                stimEventsTdt.pulse(totPulse+j-1,1) = ...
                    resamp(stimEvents.stimMsec(totPulse+j-1),sampleFreq);
                stimEventsTdt.pulse(totPulse+j-1,2) = ...
                    resamp(stimEvents.stimMsec(totPulse+j-1),lfpFreq);
                stimEventsTdt.pulse(totPulse+j-1,3) = pulseWidth;
                stimEventsTdt.pulse(totPulse+j-1,4) = pulsePeriod;
                stimEventsTdt.pulse(totPulse+j-1,5) = j; % pulse number in a train
                stimEventsTdt.pulse(totPulse+j-1,6) = diodeAct;
                stimEventsTdt.pulse(totPulse+j-1,7) = i; % stimulation number 
                stimEventsTdt.pulse(totPulse+j-1,8) = m; % train no. within a stimulation
            end
            totPulse = totPulse+nPulses;
            m = m+1;
        end
    end
        
    fullPathP = [baseFileName 'PTDT.mat'];
    save(fullPathP,'stimEventsTdt');
