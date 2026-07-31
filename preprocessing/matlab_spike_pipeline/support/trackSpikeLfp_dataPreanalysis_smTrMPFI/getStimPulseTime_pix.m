function [UpCrossings, DownCrossings] = getStimPulseTime_pix(baseFileName,nChannelsTot,sampleRate)

 % sync:  are we trying to detect sync pulses? 1: yes
 % 21 - file length (sec)
 % 26 - file size (bytes)
 
%syncCh = str2double(metaInfo(38));    % the last ch is SYNC
% ask for ID of the SYNC pulse
% default to 34, dinghao 23 Jun 2022
    
%     resp2 = 'n';
%     while strcmp(resp2, 'y') ~= 1
%         resp1 = input(['\nWhich ch is the stimulation pulse ch? Remember' ...
%             'to report the neuronscope ch # plus 1. [66]']);
%         if isempty(resp1)
%             syncCh = 66;
%         else syncCh = resp1;
%         end
%         resp2 = input(['\nIs stimulation ch # ' num2str(syncCh) '? [y/n]'], 's');
%         if isempty(resp2)
%             resp2 = 'y';
%         end
%     end
    
    syncCh = 3;  % check pix rig intan board connection if wrong

    % changed 3/2/2017, line 60-80, reading the channel directly from .dat,
    % removed the dependence on .sev file
    datFileName = [baseFileName '.dat'];
    if exist(datFileName, 'file') == 2
        listing = dir(datFileName);
        Nsamples = listing.bytes/2/nChannelsTot;  % sec
        
        fprintf('\n Loading stimulation pulses from: %s  ', datFileName);        
        fid = fopen(datFileName,'r');
         % filter before downsampling
        sync = LoadDatFile_FL_uint16(fid, syncCh, Nsamples, sampleRate, nChannelsTot);
        sync = sync';
    end
               
%     sevFileName = [baseFileName '.sev'];
%     listing = dir(sevFileName);
%     Nsamples = listing.bytes/2/nChannelsTot;  % sec
%     
%     disp('Loading stimulation pulse.....');
%     
%     fidSev = fopen(sevFileName,'r');
%     startRead = (syncCh-1)*Nsamples*2;
%     sync = LoadSevFile(fidSev,startRead,Nsamples);
    
    %[b,a]=butter(2,[900 9999]/(sampleRate/2));
    %lfSync = filtfilt(b,a,sync);
    lfSync = diff(sync);
    
    lfSync = lfSync - lfSync(1);
    % sync pulse goes first down and then up!!! This is because the sync pulse from the behavior box went through a optocoupler // changed by Yingxue
    % 20150219 
    amp = max(lfSync) - mean(lfSync(1:min(size(lfSync),sampleRate*500)));
    up = TriggerUpDownMarked(lfSync,(mean(lfSync)+amp/2),(mean(lfSync)+amp)/2);
    down = TriggerUpDownMarked(-lfSync,(mean(lfSync)+amp/2),(mean(lfSync)+amp)/2);
    
    if(up(1) > down(1) && length(up) > 1)
        up = up(2:end);
    end
    
%     for A069r-20230914-02 since one of the stim upcrossing was not
%     detected
%     up = cat(1, up(1:122), 19124543, up(123:end));  
    
    if(length(down) ~= length(up))
        disp('There are non-equal number of rising and falling edges.')
%         in = input('\nDo you want to continue? Y/N: ');
%         if(~isempty(strfind(in,'Y')) || ~isempty(strfind(in,'y')))
%             UpCrossings = up;
%             DownCrossings = down;
%         else
            UpCrossings = [];
            DownCrossings = [];
%         end
    else
        UpCrossings = up;
        DownCrossings = down;
    end
end
