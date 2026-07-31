function ConcatenateSpikes_smTr_GoodTr(path,fileName,onlyRun)
% Calculate CCG for each subsession
% 
% by Yingxue 8/25/2017

    %%%%%%%%% check arguments
    if nargin<2
        disp('At least three arguments are needed for this function.');
        return;
    elseif nargin == 2
        onlyRun = 1;
    elseif nargin > 3
        disp('Too many input arguments');        
        return;
    end
        
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    fileNameConSp = [fileName '_Concatsp_GoodTr_Run' num2str(onlyRun) '.mat'];
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
    
    %%%%%%%%% initialize constants
    GlobalConst;
    fileNameInfo = [fileName(1:end-4) '_Info.mat'];
        
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo_smTr(path,fileName);
    end
    load(fullPath);
    mazeSess = beh.mazeSessAll;
    
    lfpIndStart = getRecField(trials,'lfpIndStart',1:length(lapList));
    lfpIndEnd = getRecField(trials,'lfpIndEnd',1:length(lapList));
    if(onlyRun == 1)
        spikes = getRecField(trialsExt,'spikes',1:length(lapList));
    else
        spikes = getRecField(trials,'spikes',1:length(lapList));
    end
    clear trials trialsExt
    
    % calculate CCG for each subsession
    disp('Concatenate spikes for each subsession')
    spTrainGoodTr = cell(1,length(mazeSess));
    spCluGoodTr = cell(1,length(mazeSess));
    totLfpIndGoodTr = cell(1,length(mazeSess));
    spTrainBadTr = cell(1,length(mazeSess));
    spCluBadTr = cell(1,length(mazeSess));
    totLfpIndBadTr = cell(1,length(mazeSess));
    spTrainOKTr = cell(1,length(mazeSess));
    spCluOKTr = cell(1,length(mazeSess));
    totLfpIndOKTr = cell(1,length(mazeSess));
    if(length(mazeSess)>1)
        for i = 1:length(mazeSess)
            disp(['Session ' num2str(i)]);
            indLaps = find(beh.mazeSess == mazeSess(i));
            indLaps = intersect(indLaps,beh.indTrCtrl); 
            indLapsGoodTr = intersect(indLaps,beh.indGoodTrCtrl); 
            [spTrainGoodTr{i},spCluGoodTr{i},totLfpIndGoodTr{i}] = ...
                concatenateSpikeTrain(spikes,lfpIndStart,lfpIndEnd,indLapsGoodTr,totClu);
            indLapsBadTr = intersect(indLaps,beh.indBadTrCtrl); 
            [spTrainBadTr{i},spCluBadTr{i},totLfpIndBadTr{i}] = ...
                concatenateSpikeTrain(spikes,lfpIndStart,lfpIndEnd,indLapsBadTr,totClu);
            indLapsOKTr = setdiff(indLaps,[beh.indGoodTrCtrl,beh.indBadTrCtrl]); 
            [spTrainOKTr{i},spCluOKTr{i},totLfpIndOKTr{i}] = ...
                concatenateSpikeTrain(spikes,lfpIndStart,lfpIndEnd,indLapsOKTr,totClu);       
        end
    else
        indLapsGoodTr = beh.indGoodTrCtrl; 
        [spTrainGoodTr{1},spCluGoodTr{1},totLfpIndGoodTr{1}] = ...
            concatenateSpikeTrain(spikes,lfpIndStart,lfpIndEnd,indLapsGoodTr,totClu);
        indLapsBadTr = beh.indBadTrCtrl; 
        [spTrainBadTr{1},spCluBadTr{1},totLfpIndBadTr{1}] = ...
            concatenateSpikeTrain(spikes,lfpIndStart,lfpIndEnd,indLapsBadTr,totClu);
        indLapsOKTr = setdiff(beh.indTrCtrl,[beh.indGoodTrCtrl,beh.indBadTrCtrl]); 
        [spTrainOKTr{1},spCluOKTr{1},totLfpIndOKTr{1}] = ...
            concatenateSpikeTrain(spikes,lfpIndStart,lfpIndEnd,indLapsOKTr,totClu); 
    end

    fullPath = [path fileNameConSp];
    save(fullPath, 'spTrainGoodTr','spCluGoodTr','totLfpIndGoodTr',...
        'spTrainBadTr','spCluBadTr','totLfpIndBadTr',...
        'spTrainOKTr','spCluOKTr','totLfpIndOKTr','-v7.3');
                   
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
