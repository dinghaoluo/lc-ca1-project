function ConcatenateSpikes_smTr(path,pathsave,fileName,onlyRun)
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
    fileNameConSp = [fileName '_Concatsp_Run' num2str(onlyRun) '.mat'];
    fileNameExt = [fileName '_ext.mat'];
    fileName = [fileName '.mat'];
    
    fullPath = [path fileName]; 
    if(exist(fullPath) == 0)
        disp('The file does not exist');
        return;
    end
    load(fullPath);
    totClu = length(cluList.all);
    
    %% added on 7/22/2021 
    cd path
    path = pwd;
    indTmp = strfind(path,'\');
    path1 = path(indTmp(end-2)+1:end);
    pathNew = [pathsave '\' path1];
    if(~exist(pathNew,'dir') == 0)
        mkdir(pathNew);
    end
    
    if(onlyRun == 1)
        fullPath = [pathNew fileNameExt];
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
        
    fullPath = [pathNew fileNameInfo];
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
    spTrain = cell(1,length(mazeSess));
    spClu = cell(1,length(mazeSess));
    totLfpInd = cell(1,length(mazeSess));
    for i = 1:length(mazeSess)
        disp(['Session ' num2str(i)]);
        indLaps = find(beh.mazeSess == mazeSess(i)); 
        indLaps = intersect(indLaps, beh.indGoodLap);
        [spTrain{i},spClu{i},totLfpInd{i}] = ...
            concatenateSpikeTrain(spikes,lfpIndStart,lfpIndEnd,indLaps,totClu);
    end
    
    %% changed on 7/22/2021
    fullPath = [pathNew fileNameConSp];
    save(fullPath, 'spTrain','spClu','totLfpInd','-v7.3');
                   
    clear mydata;
    
end

function [spTrain,spClu,totLfpInd] = ...
        concatenateSpikeTrain(spikes,lfpIndStart,lfpIndEnd,indLaps,totClu)
    spTrain = [];
    spClu = [];
    totLfpInd = lfpIndEnd{max(indLaps)}...
                - lfpIndStart{min(indLaps)};
    
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
