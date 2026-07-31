function ConcatenateSpikesAlignedRunCtrlOnly(path,fileName,onlyRun,mazeSess)
% Calculate CCG for each subsession
% 
% by Yingxue 8/25/2017

    %%%%%%%%% check arguments
    if nargin<3
        disp('At least four arguments are needed for this function.');
        return;
    elseif nargin > 4
        disp('Too many input arguments');        
        return;
    end
        
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    fileNameConSp = [fileName '_ConcatspAlignedRunCtrl_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    fileNameOrig = [fileName '.mat'];
    
    fullPath = [path fileNameOrig]; 
    if(exist(fullPath) == 0)
        disp('The file does not exist');
        return;
    end
    load(fullPath,'cluList');
    totClu = length(cluList.all);
    
    fileNameRun = [fileName '_alignRun_msess' num2str(mazeSess) '.mat'];  
    fullPath = [path fileNameRun]; 
    if(exist(fullPath) == 0)
        disp('The _alignRun_msess file does not exist');
        return;
    end
    load(fullPath,'trialsRun');
    
    fullPath = [path fileName '_alignedSpikesPerNPerT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _alignedSpikesPerNPerT file does not exist');
        return;
    end
    load(fullPath,'trialsRunSpikes');
    
    %%%%%%%%% initialize constants
    GlobalConst;
    
    fullPath = [path fileName '_PeakFRAlignedCtrl_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];    
    if(exist(fullPath) == 0)
        disp('The _PeakFRAlignedCtrl file does not exist');
        return;
    end
    load(fullPath,'trialNoNonStim');
        
    lfpIndStart = trialsRun.startLfpInd;
    lfpIndEnd = trialsRun.endLfpInd;
    spikes = trialsRunSpikes.Time;
    clear trialsRun trialsRunSpikes
    
    % calculate CCG for each subsession
    disp('Concatenate spikes for nonstimulated trials')
    [spTrainNonStim,spCluNonStim,totLfpIndNonStim] = ...
            concatenateSpikeTrain(spikes,lfpIndStart,lfpIndEnd,trialNoNonStim,totClu);
       
    fullPath = [path fileNameConSp];
    save(fullPath, 'spTrainNonStim','spCluNonStim','totLfpIndNonStim','-v7.3');
                   
    clear mydata;
    
end

function [spTrain,spClu,totLfpInd] = ...
        concatenateSpikeTrain(spikes,lfpIndStart,lfpIndEnd,indLaps,totClu)
    spTrain = [];
    spClu = [];
    totLfpInd = lfpIndEnd(max(indLaps))...
                - lfpIndStart(min(indLaps));
    
    for i = 1:length(indLaps)
        for j = 1:totClu
            tmp = spikes{j,indLaps(i)};
            tmp = tmp + lfpIndStart(indLaps(i)) - lfpIndStart(indLaps(1));
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
