function check = Arduino2RecordT_smTr2P(baseFileName,sampleFreq,lfpFreq, data_2p)
% Convert arduino time to recording time
    
    check = 0; 
    
    if exist([baseFileName 'BTDT.mat'], 'file') == 2
        disp('BTDT file already exists.')
        return;
    else
        fullNameB = [baseFileName 'B.mat'];
        if(exist(fullNameB,'file') ~= 0)
            load(fullNameB);
        else
            disp('parsing the behavioral file.');
            behEvents = LoadBehMazeFile_smTr2P([baseFileName 'T.txt'],sampleFreq);

            save([baseFileName 'B.mat'], 'behEvents');
        end
             
        % time        
        framesPerMs = 1000/sampleFreq;
        behEvents.TDTsyncMsec = [];
        for n = 1:length(data_2p.F)
            behEvents.TDTsyncMsec(n) = (n-1) * framesPerMs;
        end        
        save(fullNameB,'behEvents');
               
        %% convert arduino time to recording time
        check = Arduino2TDTtime_smTr_opto2P(baseFileName,sampleFreq,lfpFreq); 
    end
