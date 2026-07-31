function [UpCrossings] = getSyncPulseTimeSmTr(baseFileName,nChannelsTot,sampleRate)

 % sync:  are we trying to detect sync pulses? 1: yes
 % 21 - file length (sec)
 % 26 - file size (bytes)
 
%syncCh = str2double(metaInfo(38));    % the last ch is SYNC
% ask for ID of the SYNC pulse    
% default to 33, dinghao 23 Jun 2022

%     resp2 = 'n';
%     while strcmp(resp2, 'y') ~= 1
%         resp1 = input(['\nWhich ch is the SYNC pulse ch? Remember to',...
%                        ' report the neuronscope ch # plus 1. [65]']);
%         if isempty(resp1)
%             syncCh = 65;
%         else syncCh = resp1;
%         end
%         resp2 = input(['\nIs SYNC ch # ' num2str(syncCh) '? [y/n]'], 's');
%         if isempty(resp2)
%             resp2 = 'y';
%         end
%     end

    syncCh = 33;
    
    % in case this is a hippocampus recording, Dinghao, 20 Jan 2023
    if exist([baseFileName '.fet.7'], 'file')==2
        syncCh = 65;
    end
    
    % for pixel rig recording concatenated, Dinghao, 12 July 2023
    if exist([baseFileName '-concat.dat'], 'file')==2
        syncCh = 65;
    end
    
    % changed 3/2/2017, line 22-44, reading the channel directly from .dat,
    % removed the dependence on .sev file
    datFileName = [baseFileName '.dat'];
    
    % for pixel rig recording concatenated, Dinghao, 12 July 2023
    if exist([baseFileName '-concat.dat'], 'file')==2
        datFileName = [baseFileName '-concat.dat'];
    end
    
    if exist(datFileName, 'file') == 2
        listing = dir(datFileName);
        Nsamples = listing.bytes/2/nChannelsTot;  % sec
        
        fprintf('\n Loading SYNC pulses from: %s  ', datFileName);        
        fid = fopen(datFileName,'r');
         % filter before downsampling
        sync = LoadDatFile_FL(fid, syncCh, Nsamples, sampleRate, nChannelsTot);
        sync = sync';
    end
    
    [b,a]=butter(2,[10 9000]/(sampleRate/2));
    lfSync = filtfilt(b,a,sync);
    
    lfSync = lfSync - lfSync(1);
    % sync pulse goes first down and then up!!! This is because the sync 
    % pulse from the behavior box went through a optocoupler 
    % changed by Yingxue 20150219 
    amp = max(lfSync) - min(lfSync);
    [up down] = TriggerUpDownMarked(lfSync,(mean(lfSync)+amp/6),(mean(lfSync)+amp)/6);
    
    %%% changed by Yingxue 20150219, check the pulse width and 
    % guarantee that the down and up has equal length 
    % remove any pulses which do not have both up and down edges
    if(~isempty(down) && ~isempty(up))
        lenDown = length(down)-1;
        lenUp = length(up);
        if(lenUp - lenDown > 0)
            lenTmp = lenDown;
        else
            lenTmp = lenUp;
        end

        upNew = zeros(lenTmp,1);
        downNew = zeros(lenTmp,1);
        lengthUp = 1;
        for i = 1:lenTmp
            if(lenUp - lenDown > 0)
                tmpDiff = down(i) - up;
            else
                tmpDiff = down - up(i);
            end
            tmpMatch = find(tmpDiff > sampleRate*0.004 & tmpDiff < sampleRate*0.014);
            if length(tmpMatch) == 0
                disp('here');
            end
            if(length(tmpMatch) == 2 && diff(tmpMatch) == 1)
                tmpMatch = tmpMatch(1);
            end
            if(length(tmpMatch) == 1)
                if(lenUp - lenDown > 0)
                    upNew(lengthUp) = up(tmpMatch);
                    downNew(lengthUp) = down(i);
                else
                    upNew(lengthUp) = up(i);
                    downNew(lengthUp) = down(tmpMatch);
                end
                lengthUp = lengthUp + 1;
            else
                disp('sync pulse width has a problem.');
            end
        end
        upNew = upNew(1:lengthUp-1);
        downNew = downNew(1:lengthUp-1);
        UpCrossings = upNew([sampleRate;diff(upNew)]>sampleRate*0.9);
    else
        UpCrossings = [];
    end
    %%% changed by Yingxue 20150219 

end
