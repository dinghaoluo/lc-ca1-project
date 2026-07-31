function ConcatenateSpikes_smTr_CtrlOnly(path,fileName,onlyRun,mazeSess)
% Calculate CCG for each subsession
% 
% by Yingxue 8/25/2017

    %%%%%%%%% check arguments
    if nargin<2
        disp('At least three arguments are needed for this function.');
        return;
    elseif nargin == 2
        onlyRun = 1;
        mazeSess = 1;
    elseif nargin == 3
        mazeSess = 1;
    elseif nargin > 4
        disp('Too many input arguments');        
        return;
    end
        
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    fileNameFR = [fileName '_FR_Ctrl_Run' num2str(onlyRun) '_mazeSess' num2str(mazeSess) '.mat'];
    fileNameConSp = [fileName '_Concatsp_Ctrl_Run' num2str(onlyRun) '_mazeSess' num2str(mazeSess) '.mat'];
    fileNameExt = [fileName '_ext.mat'];
    fileName = [fileName '.mat'];
    
    fullPath = [path fileName]; 
    if(exist(fullPath) == 0)
        disp('The file does not exist');
        return;
    end
    load(fullPath);
    totClu = length(cluList.all);
    
    if(onlyRun == 1)
        fullPath = [path fileNameExt];
        if(exist(fullPath) == 0)
            disp(['Extended file does not exist. Please run function',...
                  ' SpikeDuringRun first']);
            return;
        end
        load(fullPath);
    end
    
    fullPath = [path fileNameFR];
    if(exist(fullPath,'file') == 0)
        disp('_FR_Ctrl file does not exist. Please run "MeanFiringRateCtrlOnly" function first');
    end
    load(fullPath,'mFRStructSessCtrl');
    
    %%%%%%%%% initialize constants
    GlobalConst;
    fileNameInfo = [fileName(1:end-4) '_Info.mat'];
        
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo_smTr(path,fileName);
    end
    load(fullPath);
    
    lfpIndStart = getRecField(trials,'lfpIndStart',1:length(lapList));
    lfpIndEnd = getRecField(trials,'lfpIndEnd',1:length(lapList));
    if(onlyRun == 1)
        spikes = getRecField(trialsExt,'spikes',1:length(lapList));
    else
        spikes = getRecField(trials,'spikes',1:length(lapList));
    end
    clear trials trialsExt
    
    % calculate CCG for each subsession
    disp('Concatenate spikes for ctrl trials');
    indLaps = mFRStructSessCtrl.indLapList; 
    [spTrainCtrl,spCluCtrl,totLfpIndCtrl] = ...
            concatenateSpikeTrain(spikes,lfpIndStart,lfpIndEnd,indLaps,totClu);
        
    fullPath = [path fileNameConSp];
    save(fullPath, 'spTrainCtrl','spCluCtrl','totLfpIndCtrl','-v7.3');
                   
    clear mydata;
    
end

function [spTrain,spClu,totLfpInd] = ...
        concatenateSpikeTrain(spikes,lfpIndStart,lfpIndEnd,indLaps,totClu)
    spTrain = [];
    spClu = [];
    if(~isempty(indLaps))
        totLfpInd = lfpIndEnd{max(indLaps)}...
                - lfpIndStart{min(indLaps)};
    else
        totLfpInd = 0;
        return;
    end
    
    for i = 1:length(indLaps)
        for j = 1:totClu
            tmp = spikes{indLaps(i)}{j};
            tmp = tmp + lfpIndStart{indLaps(i)} - lfpIndStart{indLaps(1)};
            spTrain = [spTrain;tmp];
            spClu = [spClu;j*ones(length(tmp),1)];
        end
    end
    
    for i = 1:totClu
        spNum = sum(spClu == i);
        if(spNum == 0)
            spTrain = [spTrain;1];
            spClu = [spClu;i];
        end
    end
    [spTrain,ind] = sort(spTrain);
    spClu = spClu(ind);
end
