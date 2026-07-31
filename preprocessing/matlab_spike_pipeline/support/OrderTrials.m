function OrderTrials(path,fileName,onlyRun)
% plot field according to licking time within a trial
% sortMethod: 1: based on the first lick position
%             2: based on the running distance of the first segment
%             3: based on the running speed of the first segment

    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    
    fileNameBeh = [fileName '_orderTrials' ...
                    '_Run' num2str(onlyRun) '.mat'];                    
    fileNameSpInfo = [fileName '_SpInfo_Run' num2str(onlyRun) '.mat'];                    
    fileNameInfo = [fileName '_Info.mat'];
    fileNameFR = [fileName '_FR_Run' num2str(onlyRun) '.mat'];
    fileNameDepth = [fileName '_Depth.mat'];
    fileNameSpeed = [fileName '_runSpeed.mat'];
    fileNameRec = [fileName '.mat'];
    
    fullPath = [path fileNameSpInfo];
    if(exist(fullPath) == 0)
        disp(['The spatial information file does not exist. Please call ',...
                'function "GetFRMapInfo" first.']);
        return;
    end
        
    fullPath = [path fileNameInfo];
    if(exist(fullPath,'file') == 0)
        BasicInfo(path,fileName);
    end
    load(fullPath);
    mazeSess = beh.mazeSessAll;
        
    fullPath = [path fileNameFR];
    if(exist(fullPath) == 0)
        disp(['The neuron depth file does not exist. Please call ',...
                'function "GetNeuRelativeDepth" first.']);
        return;
    end
    load(fullPath,'mFRStruct');
    
    fullPath = [path fileNameDepth];
    if(exist(fullPath) == 0)
        disp(['The mean firing rate file does not exist. Please call ',...
                'function "MeanFiringRateVR" first.']);
        return;
    end
    load(fullPath);
    
    fullPath = [path fileNameSpeed];
    if(exist(fullPath) == 0)
        disp(['The run speed file does not exist. Please call ',...
                'function "RunSpeed" first.']);
        return;
    end
    load(fullPath,'runSegments','speedSpectro');
    
    fullPath = [path fileNameRec];
    if(exist(fullPath) == 0)
        disp('The recording file does not exist.');
        disp(fullPath)
        return;
    end
    load(fullPath,'trials','cluList');
    
    indNeurons = 1:rec.numNeurons;
    indGoodLap = beh.indGoodLap;
    %% sort trials 
    mazeSessTr = beh.mazeSess;
    % 1. based on the first lick position
    % find first lick position
    orderByLick = struct(...
            'firstLickPos', [],...
            'firstLickPosOrdered',[],...
            'indGoodLapOrdered', []);
        
    orderByLick.firstLickPos = -1*ones(1,length(mazeSessTr));
    for i = 1:length(mazeSessTr)
        if(sum(indGoodLap == i)==0)
            continue;
        end
        lickDist = trials{i}.xMM(trials{i}.lickLfpInd);
        ind = find(lickDist > 300, 1, 'first'); 
        orderByLick.firstLickPos(i) = mean(lickDist(ind));
    end

    % order trials by where the lick occurs
    for i = mazeSess'
        ind = find(mazeSessTr == i);
        ind = intersect(ind,indGoodLap);
        [orderByLick.firstLickPosOrdered{i}, indGoodLapOrdered] = ...
            sort(orderByLick.firstLickPos(ind));
        orderByLick.indGoodLapOrdered{i} = ind(indGoodLapOrdered)';
    end
    orderByLick.firstLickPos  = orderByLick.firstLickPos(indGoodLap);
    
    % 2. based on the running distance of the longest segment
    orderByDistMaxRunSeg = struct(...
            'distMaxRunSegment', runSegments.distMaxRunSegment,...
            'distMaxRunSegOrdered',[],...
            'indGoodLapOrdered', []);
    for i = mazeSess'
        ind = find(mazeSessTr(indGoodLap) == i);
        [orderByDistMaxRunSeg.distMaxRunSegOrdered{i}, indGoodLapOrdered] = ...
            sort(runSegments.distMaxRunSegment(ind));
        orderByDistMaxRunSeg.indGoodLapOrdered{i} = ind(indGoodLapOrdered)';
    end
    
    % 3. based on the running speed of the longest segment
    orderBySpeedMaxRunSeg = struct(...
            'speedMaxRunSegment', runSegments.speedMaxRunSegment,...
            'speedMaxRunSegOrdered',[],...
            'indGoodLapOrdered', []);
    for i = mazeSess'
        ind = find(mazeSessTr(indGoodLap) == i);
        [orderBySpeedMaxRunSeg.speedMaxRunSegOrdered{i}, indGoodLapOrdered] = ...
            sort(runSegments.speedMaxRunSegment(ind));
        orderBySpeedMaxRunSeg.indGoodLapOrdered{i} = ind(indGoodLapOrdered)';
    end
    
    speedSpectrogram.mmeanAmpSpectro3_5 = speedSpectro.mmeanAmpSpectro3_5;
    speedSpectrogram.mmeanAmpSpectro6_10 = speedSpectro.mmeanAmpSpectro6_10;
    speedSpectrogram.percHighPower3_5 = speedSpectro.percHighPower3_5;
    speedSpectrogram.percHighPower6_10 = speedSpectro.percHighPower6_10;
        
    fullName = [path fileNameBeh];
    save(fullName, 'orderByLick','orderByDistMaxRunSeg',...
                   'orderBySpeedMaxRunSeg','speedSpectrogram');
end
