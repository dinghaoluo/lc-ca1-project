function GetFRMapInfo_CtrlOnly(path, fileName, onlyRun, mazeSess)
% calculate spatial information
%
% by Yingxue, 2017.08.24
    
    if nargin<2
        disp('At least three arguments are needed for this function.');
        return;
    elseif(nargin == 2)
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
    
    fileNameInfo = [fileName '_Info.mat'];
    fileNameExt = [fileName '_ext.mat'];
    fileNameSpInfo = [fileName '_SpInfo_Ctrl_Run' num2str(onlyRun) '_mazeSess' num2str(mazeSess) '.mat'];
    fileNameFR = [fileName '_FR_Ctrl_Run' num2str(onlyRun) '_mazeSess' num2str(mazeSess) '.mat'];
    fileName = [fileName '.mat'];
    
    fullPath = [path fileName];
    if(exist(fullPath) == 0)
        disp('The recording file does not exist');
        return;
    end
    load(fullPath);
    totClu = length(cluList.all);
    
    fullPath = [path fileNameFR];
    if(exist(fullPath,'file') == 0)
        disp('_FR_Ctrl file does not exist. Please run MeanFiringRateCtrlOnly function first');
    end
    load(fullPath,'mFRStructSessCtrl');
    
    if(onlyRun == 1)
        fullPath = [path fileNameExt];
        if(exist(fullPath) == 0)
            disp(['Extended file does not exist. Please run function',...
                  ' SpikeDuringRun first']);
            return;
        end
        load(fullPath);
    end
    
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo_smTr(path,fileName);
    end
    load(fullPath);
    
    GlobalConst; 
    
    param.divDist = 5; % mm
    param.smooth = 10; % mm
    param.sampleFq = sampleFq;

    drawFig = 0;
    if(drawFig == 1)
        figure;
    end
    
    if(onlyRun == 1)
        distSp = getRecField(trialsExt,'spikesMM',1:length(lapList));
    else
        distSp = getRecField(trials,'spikesMM',1:length(lapList));
    end
    
    spatialInfoSessCtrl = struct('meanFR',zeros(1,totClu),... % mean firing rate
                         'smoothedFR',[],... % binned and smoothed firing rate map
                         'smoothedBinProb',[],... % smoothed binProb
                         'normSmoothedFR',[],... % normalized smoothedFR
                         'binnedSpikes',[],... % amount of spikes occurred in each bin 
                         'binnedTime',[],... % amount of time spend in each bin
                         'binProb',[],... % percent of time spent in each bin
                         'spatialInfo',zeros(1,totClu),... % spatial information
                         'adaptSpatialInfo',zeros(1,totClu),... % spatial information after adaptive smoothing
                         'adaptSmoothedFR',[],... % FR after adaptive smoothing
                         'sparsity',zeros(1,totClu),... % sparsity
                         'SNR',zeros(1,totClu)); % SNR     
    
    disp('Calculate spatial information for ctrl tirals')
    indLaps = mFRStructSessCtrl.indLapList;
        
    % Collect the distance information over all the trials within a
    % subsession
    distTrCtrl = [];
    for j = 1:length(indLaps)
        if(onlyRun == 1)
            ind = trials{indLaps(j)}.speed > minSpeed;
        else
            ind = 1:trials{indLaps(j)}.Nsamples;
        end
        distTmp = trials{indLaps(j)}.xMM(ind);
        distTrCtrl = [distTrCtrl;distTmp];
    end

    % Collect the distance information over all the trials when spikes
    % occurs
    distSpikesCtrl = ...
        concatenateSpikes(distSp,indLaps,totClu);
    
    for j = 1:totClu
        if(~isempty(distSpikesCtrl{j}))
            spatialInfoSessCtrl.meanFR(j) = length(distSpikesCtrl{j}) / ...
                (length(distTrCtrl)/sampleFq);

            % Bin the spikes over space, and smooth the curve after binning
            [spatialInfoSessCtrl.smoothedFR(j,:), ...
                spatialInfoSessCtrl.smoothedBinProb(j,:), ...
                spatialInfoSessCtrl.normSmoothedFR(j,:),...
                spatialInfoSessCtrl.binnedSpikes(j,:),...
                spatialInfoSessCtrl.binnedTime(j,:),...
                spatialInfoSessCtrl.binProb(j,:)] =...
                binCoordiate(distSpikesCtrl{j},distTrCtrl,param);

            spatialInfoSessCtrl.spatialInfo(j) = ...
                getSpInfo1(spatialInfoSessCtrl.smoothedFR(j,:), ...
                    spatialInfoSessCtrl.smoothedBinProb(j,:),...
                    spatialInfoSessCtrl.meanFR(j));

            [spatialInfoSessCtrl.adaptSpatialInfo(j),...
                spatialInfoSessCtrl.adaptSmoothedFR(j,:)] = ...
                getSpInfo_addaptBin1(spatialInfoSessCtrl.binnedSpikes(j,:),...
                    spatialInfoSessCtrl.binnedTime(j,:), ...
                    spatialInfoSessCtrl.binProb(j,:),...
                    spatialInfoSessCtrl.meanFR(j));

            [spatialInfoSessCtrl.sparsity(j), spatialInfoSessCtrl.SNR(j)] = ...
                sparsityInfo(spatialInfoSessCtrl.smoothedBinProb(j,:), ...
                spatialInfoSessCtrl.smoothedFR(j,:),...
                spatialInfoSessCtrl.meanFR(j));
        else
            spatialInfoSessCtrl.adaptSpatialInfo(j) = nan;
            spatialInfoSessCtrl.spatialInfo(j) = nan;
            spatialInfoSessCtrl.sparsity(j) = nan;
            spatialInfoSessCtrl.SNR(j) = nan;
        end  
    end
    
    fullPath = [path fileNameSpInfo];
    save(fullPath, 'spatialInfoSessCtrl');
end

function distSpikes = ...
        concatenateSpikes(distSp,indLaps,totClu)
    
    distSpikes = cell(1,totClu);
    
    for i = 1:length(indLaps)
        for j = 1:totClu
            distSpikes{j} = [distSpikes{j};distSp{indLaps(i)}{j}];
        end
    end
end
