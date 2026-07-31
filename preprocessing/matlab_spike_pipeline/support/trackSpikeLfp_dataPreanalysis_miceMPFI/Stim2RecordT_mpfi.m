function Stim2RecordT_mpfi(baseFileName,sampleFreq,lfpFreq,nChannelsTot)
% Convert stimulation time to recording time

    if exist([baseFileName 'PTDT.mat'], 'file') == 2
        disp('BTDT file already exists.')
        return;
    else
        fullNameP = [baseFileName 'P.mat'];
        if(exist(fullNameP,'file') ~= 0)
            load(fullNameP);
        else
            disp('parsing the stimulation file.');
            stimEvents = LoadStimFile([baseFileName 'PulseE.txt']);
                      
            save([baseFileName 'P.mat'], 'stimEvents');
        end

        if(~isfield(stimEvents,'stimInd'))
            % get the recording time of the stimulation pulses
            
            [UpCrossings] = getPulseTime(baseFileName,nChannelsTot,sampleFreq);
            if(~isempty(UpCrossings))

                % store the sync pulse recording time into Arduino behavior file
                stimEvents.stimInd = UpCrossings;
                stimEvents.stimMsec = UpCrossings / sampleFreq * 1000;
                
                save(fullNameP,'stimEvents');
            end
        end

        % convert stimulation time to recording time        
        Arduino2TDTtimeStim_mpfi(baseFileName,sampleFreq,lfpFreq);
    end
end
    
function [UpCrossings] = getPulseTime(baseFileName,nChannelsTot,sampleRate)

 % sync:  are we trying to detect sync pulses? 1: yes
 % 21 - file length (sec)
 % 26 - file size (bytes)
 
    %syncCh = str2double(metaInfo(38));    % the last ch is SYNC
    % ask for ID of the SYNC pulse
    resp2 = 'n';
    while strcmp(resp2, 'y') ~= 1
        resp1 = input(['\nWhich ch is the stimulation pulse ch? Remember' ...
            'to report the neuronscope ch # plus 1. [228]']);
        if isempty(resp1)
            syncCh = 228;
        else syncCh = resp1;
        end
        resp2 = input(['\nIs stimulation ch # ' num2str(syncCh) '? [y/n]'], 's');
        if isempty(resp2)
            resp2 = 'y';
        end
    end
    
    % changed 3/2/2017, line 60-80, reading the channel directly from .dat,
    % removed the dependence on .sev file
    datFileName = [baseFileName '.dat'];
    if exist(datFileName, 'file') == 2
        listing = dir(datFileName);
        Nsamples = listing.bytes/2/nChannelsTot;  % sec
        
        fprintf('\n Loading stimulation pulses from: %s  ', datFileName);        
        fid = fopen(datFileName,'r');
         % filter before downsampling
        sync = LoadDatFile_FL(fid, syncCh, Nsamples, sampleRate, nChannelsTot);
        sync = sync';
    end
                   
    [b,a]=butter(2,[10 9000]/(sampleRate/2));
    lfSync = filtfilt(b,a,sync);
    
    lfSync = lfSync - lfSync(1);
    % sync pulse goes first down and then up!!! This is because the sync pulse
    % from the behavior box went through a optocoupler 
    % changed by Yingxue 20150219 
    amp = max(lfSync) - mean(lfSync(1:sampleRate*500));
    [up down] = TriggerUpDownMarked(lfSync,(mean(lfSync)+amp/2),(mean(lfSync)+amp)/2);
    
    if(length(down) ~= length(up))
        disp('There are non-equal number of rising and falling edges.')
        in = input('\nDo you want to continue? Y/N: ');
        if(~isempty(strfind(in,'Y')) || ~isempty(strfind(in,'y')))
            UpCrossings = up;
        else
            UpCrossings = [];
        end
    else
        UpCrossings = up;
    end
end
