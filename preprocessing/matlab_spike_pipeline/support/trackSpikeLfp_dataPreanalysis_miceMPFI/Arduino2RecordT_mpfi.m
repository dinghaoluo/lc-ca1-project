function Arduino2RecordT_mpfi(baseFileName,sampleFreq,lfpFreq,nChannelsTot)
% Convert arduino time to recording time

    if exist([baseFileName 'BTDT.mat'], 'file') == 2
        disp('BTDT file already exists.')
        return;
    else
        fullNameB = [baseFileName 'B.mat'];
        if(exist(fullNameB,'file') ~= 0)
            load(fullNameB);
        else
            disp('parsing the behavioral file.');
            behEvents = LoadBehMazeFile_mice([baseFileName 'T.txt']);
                      
            save([baseFileName 'B.mat'], 'behEvents');
        end
        
        if(~isfield(behEvents,'TDTsyncInd'))
            % get the recording time of the sync pulses
            
            [UpCrossings] = getSyncPulseTime(baseFileName,nChannelsTot,sampleFreq); 
          
            % store the sync pulse recording time into Arduino behavior file
            behEvents.TDTsyncInd = UpCrossings;
            behEvents.TDTsyncMsec = UpCrossings / sampleFreq * 1000;

            save(fullNameB,'behEvents');
        end

        % convert arduino time to recording time
        Arduino2TDTtime_mpfi(baseFileName,sampleFreq,lfpFreq); 
    end
