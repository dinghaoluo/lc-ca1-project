function ConcatenateSpikesAlignedRun(path,fileName,onlyRun,mazeSess)
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
    fileNameConSp = [fileName '_ConcatspAlignedRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
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
    
    fullPath = [path fileName '_PeakFRAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];    
    if(exist(fullPath) == 0)
        disp('The _PeakFRAligned file does not exist');
        return;
    end
    load(fullPath,'trialNoNonStimGood','trialNoNonStimBad','trialNoStim','trialNoStimCtrl','pulseMeth');
        
    lfpIndStart = trialsRun.startLfpInd;
    lfpIndEnd = trialsRun.endLfpInd;
    spikes = trialsRunSpikes.Time;
    clear trialsRun trialsRunSpikes
    
    % calculate CCG for each subsession
    disp('Concatenate spikes for nonstimulated good trials')
    [spTrainNonStimGood,spCluNonStimGood,totLfpIndNonStimGood] = ...
            concatenateSpikeTrain(spikes,lfpIndStart,lfpIndEnd,trialNoNonStimGood,totClu);
        
    disp('Concatenate spikes for nonstimulated Bad trials')
    [spTrainNonStimBad,spCluNonStimBad,totLfpIndNonStimBad] = ...
            concatenateSpikeTrain(spikes,lfpIndStart,lfpIndEnd,trialNoNonStimBad,totClu);
       
    disp('Concatenate spikes for stimulated trials')
    %% changed by Yingxue on 7/14/2021
    spTrainStim = [];
    spCluStim = [];
    totLfpIndStim = [];
    spTrainStimCtrl = [];
    spCluStimCtrl = [];
    totLfpIndStimCtrl = [];
    for i = 1:length(pulseMeth)
        [spTrainStim{i},spCluStim{i},totLfpIndStim{i}] = ...
            concatenateSpikeTrain(spikes,lfpIndStart,lfpIndEnd,trialNoStim{i},totClu);
        
        [spTrainStimCtrl{i},spCluStimCtrl{i},totLfpIndStimCtrl{i}] = ...
            concatenateSpikeTrain(spikes,lfpIndStart,lfpIndEnd,trialNoStimCtrl{i},totClu);
    end
        
    fullPath = [path fileNameConSp];
    save(fullPath, 'spTrainNonStimGood','spCluNonStimGood','totLfpIndNonStimGood',...
        'spTrainNonStimBad','spCluNonStimBad','totLfpIndNonStimBad',...
        'spTrainStim','spCluStim','totLfpIndStim',...
        'spTrainStimCtrl','spCluStimCtrl','totLfpIndStimCtrl',...
        '-v7.3');
                   
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
