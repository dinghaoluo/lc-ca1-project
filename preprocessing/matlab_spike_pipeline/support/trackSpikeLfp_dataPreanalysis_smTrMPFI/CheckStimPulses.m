function [behEvents,status] = CheckStimPulses(behEvents,UpCrossings,DownCrossings,sampleFreq)
% check every stimulation pulse for their pulse width

    status = 0;

    % number of pulses per stimulation
    totNumPulses = 0;
    numPC = 0;
    numPulses = 1;
    if(isfield(behEvents,'stimOn'))
        numPC = length(behEvents.stimOn);
        numPulsePerStim = zeros(1,numPC);
        for i = 1:numPC
            % modified to accommodate small inaccuracies, but source
            % unfound, Dinghao 23 Jan 2023
            indTmp = find(behEvents.pulsePar(:,1) <= behEvents.stimOn(i,1)+1250,1,'last');
            if(~isempty(indTmp))
                numPulsePerStim(i) =  behEvents.pulsePar(indTmp,5) * ...
                    behEvents.pulsePar(indTmp,7);
                totNumPulses = totNumPulses + numPulsePerStim(i);
                [behEvents.stimOnPar(i,:),behEvents.stimOnDiode{i}] = recStimPar(behEvents.pulsePar(indTmp,:),...
                    behEvents.diodePar(indTmp,:));

            end
        end
    
        % check each optogenetic pulses
        if(totNumPulses <= length(UpCrossings))
            for i = 1:numPC
                index = numPulses:numPulses+numPulsePerStim(i)-1;           
                DC = DownCrossings(index);
                UC = UpCrossings(index);

                behEvents.stimOn(i,2) = ...
                    checkPulseW(DC,UC,behEvents.stimOnPar(i,:),...
                        sampleFreq,numPulsePerStim(i));
                behEvents.stimOn(i,3) = numPulses; % starting pulse index of the current stimulation
                behEvents.stimOn(i,4) = numPulsePerStim(i); % number of pulses within the current stimulation
                if(behEvents.stimOn(i,2) == 0)
                    disp(['The recorded stimulation pulse width does not'...
                        ' match its parameter. PC number : ' num2str(i)]);
                    status = 1;
                    return;
                end
                numPulses = numPulses+numPulsePerStim(i);
            end
        end
    end

    % check whether there are stimulation pulses not recorded in the
    % behavior file
    if(totNumPulses < length(UpCrossings))
        if(isfield(behEvents,'stimOn'))
            indTmp = find(behEvents.pulsePar(:,1) > behEvents.stimOn(end,1));
        else
            indTmp = 1:length(behEvents.pulsePar(:,1));
        end
        if(~isempty(indTmp))
            totExtNP = 0;
            numPC1 = length(indTmp);
            for i = 1:numPC1
                totExtNP = totExtNP + behEvents.pulsePar(indTmp(i),5);
            end
            if(totNumPulses+totExtNP <= length(UpCrossings))
                numPulsePerStim = zeros(1,numPC1);
                numPulses = totNumPulses+1;
                for i = 1:numPC1
                    [behEvents.stimOnPar(numPC+i,:),...
                        behEvents.stimOnDiode{numPC+i}]= ...
                        recStimPar(behEvents.pulsePar(indTmp(i),:),...
                                    behEvents.diodePar(indTmp(i),:));
                    numPulsePerStim(i) = behEvents.pulsePar(indTmp(i),5) * ...
                        behEvents.pulsePar(indTmp(i),7 );

                    index = numPulses:numPulses+numPulsePerStim(i)-1;
                    DC = DownCrossings(index);
                    UC = UpCrossings(index);
                    behEvents.stimOn(numPC+i,1) = -1;
                    behEvents.stimOn(numPC+i,2) = ...
                        checkPulseW(DC,UC,behEvents.stimOnPar(numPC+i,:),...
                            sampleFreq,numPulsePerStim(i));
                    behEvents.stimOn(numPC+i,3) = numPulses;
                    behEvents.stimOn(numPC+i,4) = numPulsePerStim(i);

                    if(behEvents.stimOn(i,2) == 0)
                        disp(['The recorded stimulation pulse width does not'...
                            ' match its parameter. PP number : ' num2str(indTmp(i))]);
                        status = 1;
                        return;
                    end

                    numPulses = numPulses+numPulsePerStim(i);
                    disp(['tagging pulse ' num2str(i) ', total number of pulses ' num2str(numPulses)]);
                end
            end
        end
    elseif(totNumPulses > length(UpCrossings))
        disp(['Unmatching numbers of stimulation pulses between the '...
            'behavior file and the recording. Please check.']);
        status = 1;
        return;
    end
end

function [stimOnPar,stimOnDiode] = recStimPar(pulsePar,diodePar)
    stimOnPar(1) = pulsePar(5); % number of pulses
    stimOnPar(2) = pulsePar(3)/1000; % pulse width in ms
    stimOnPar(3) = pulsePar(4)/1000; % pulse period in ms
    stimOnPar(4) = pulsePar(6)/1000; % off time between pulse trains in ms
    stimOnPar(5) = pulsePar(7); % repeats
    indDiode = find(diodePar(2:end) > 0);
    stimOnDiode.indDiode = indDiode; % which diode is selected
    stimOnDiode.diodeCurr = diodePar(indDiode+1); % diode current
end

function stimOn = checkPulseW(DownCrossings,UpCrossings,stimOnPar,sampleFreq,...
                                numPulsePerStim)
    stimOn = 0; % 0: bad stimulation; 1: good stimulation
    pulseW = (DownCrossings-UpCrossings)/sampleFreq*1000;
    checkPW = pulseW > stimOnPar(2)*0.5 & pulseW < stimOnPar(2)*1.5;
    if(sum(checkPW) == numPulsePerStim) % check pulse width
        if(length(UpCrossings) == 1)
            stimOn = 1;
            return;
        end
        pulsePer = diff(UpCrossings)/sampleFreq*1000;
        checkPP = pulsePer >= stimOnPar(3)*0.5 & pulsePer <= stimOnPar(3)*1.5; 
            % check pulse period
        if(stimOnPar(5) == 1) % repeats == 1
            if(sum(checkPP) == numPulsePerStim-1)
                stimOn = 1; % good stimulation
            end
        else
            interval = stimOnPar(4)+ stimOnPar(3);
            checkPP1 = pulsePer > interval*0.89 & pulsePer < interval*1.11; 
                % check pulse period
            if(sum(checkPP) == ...
                    ((numPulsePerStim-1)*stimOnPar(5)) &&...
                    sum(checkPP1) == (stimOnPar(5)-1))
                stimOn = 1; % good stimulation
            end
        end
        % tagging (60), Dinghao 15 Feb 2024
        if(numPulsePerStim == 60 || numPulsePerStim == 10)
            stimOn = 1;
        end
    end
end
