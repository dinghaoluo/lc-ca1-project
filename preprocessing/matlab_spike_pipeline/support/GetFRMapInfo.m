function GetFRMapInfo(path, fileName, onlyRun)
% calculate spatial information
%
% by Yingxue, 2017.08.24
    
    if nargin<2
        disp('At least three arguments are needed for this function.');
        return;
    elseif(nargin == 2)
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
    
    fileNameInfo = [fileName '_Info.mat'];
    fileNameExt = [fileName '_ext.mat'];
    fileNameSpInfo = [fileName '_SpInfo_Run' num2str(onlyRun) '.mat'];
    fileName = [fileName '.mat'];
    
    fullPath = [path fileName];
    if(exist(fullPath) == 0)
        disp('The recording file does not exist');
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
    
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo_smTr(path,fileName);
    end
    load(fullPath);
    mazeSess = beh.mazeSessAll;
    
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
    
    spatialInfo = struct('meanFR',zeros(1,totClu),... % mean firing rate
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
    
    disp('Calculate spatial information for each subsession')
    spatialInfoSess = cell(1,length(mazeSess));
    for i = 1:length(mazeSess)    
        disp(['Session ' num2str(i)]);
        spatialInfoSess{i} = spatialInfo;
        indLaps = find(beh.mazeSess == mazeSess(i)); 
        indLaps = intersect(indLaps, beh.indGoodLap);
        
        % Collect the distance information over all the trials within a
        % subsession
        distTr = [];
        for j = 1:length(indLaps)
            if(onlyRun == 1)
                ind = trials{indLaps(j)}.speed > minSpeed;
            else
                ind = 1:trials{indLaps(j)}.Nsamples;
            end
            distTmp = trials{indLaps(j)}.xMM(ind);
            distTr = [distTr;distTmp];
        end
        
        % Collect the distance information over all the trials when spikes
        % occurs
        distSpikes = ...
            concatenateSpikes(distSp,indLaps,totClu);
        
        for j = 1:totClu
            if(~isempty(distSpikes{j}))
                spatialInfoSess{i}.meanFR(j) = length(distSpikes{j}) / ...
                    (length(distTr)/sampleFq);

                % Bin the spikes over space, and smooth the curve after binning
                [spatialInfoSess{i}.smoothedFR(j,:), ...
                    spatialInfoSess{i}.smoothedBinProb(j,:), ...
                    spatialInfoSess{i}.normSmoothedFR(j,:),...
                    spatialInfoSess{i}.binnedSpikes(j,:),...
                    spatialInfoSess{i}.binnedTime(j,:),...
                    spatialInfoSess{i}.binProb(j,:)] =...
                    binCoordiate(distSpikes{j},distTr,param);

                spatialInfoSess{i}.spatialInfo(j) = ...
                    getSpInfo1(spatialInfoSess{i}.smoothedFR(j,:), ...
                        spatialInfoSess{i}.smoothedBinProb(j,:),...
                        spatialInfoSess{i}.meanFR(j));

                [spatialInfoSess{i}.adaptSpatialInfo(j),...
                    spatialInfoSess{i}.adaptSmoothedFR(j,:)] = ...
                    getSpInfo_addaptBin1(spatialInfoSess{i}.binnedSpikes(j,:),...
                        spatialInfoSess{i}.binnedTime(j,:), ...
                        spatialInfoSess{i}.binProb(j,:),...
                        spatialInfoSess{i}.meanFR(j));

                [spatialInfoSess{i}.sparsity(j), spatialInfoSess{i}.SNR(j)] = ...
                    sparsityInfo(spatialInfoSess{i}.smoothedBinProb(j,:), ...
                    spatialInfoSess{i}.smoothedFR(j,:),...
                    spatialInfoSess{i}.meanFR(j));

                if drawFig == 1
                    nl=1;
                    nc=2;
                    subplot(nl,nc,1); cla;
                    plot(1:distTr(1:2:end),ones(1,length(1:distTr(1:2:end))),...
                        '.', 'MarkerSize', 3,'Color', [0.6 0.6 0.6]);
                    hold on;
                    plot(distSpikes{j},ones(1,length(distSpikes{j})),...
                        '.','MarkerSize', 3,'Color', [1 0 0]);
                    xlim([0 max(distTr)]);

                    title(['totClu:' num2str(j) '     sh:' num2str(cluList.shank(j))...
                        '  locClu:' num2str(cluList.localClu(j))]);

                    % place fields
                    subplot(nl,nc,2); cla;
                    trBin = round(distTr/param.divDist)+1;
                    bins = unique(trBin);
                    maxC = max(spatialInfoSess{i}.smoothedFR(j,:));
                    imagesc(bins'*param.divDist, ones(1,length(bins)),...
                        spatialInfoSess{i}.smoothedFR(j,:),[0 maxC]); 
                    hold on;
                    set(gca,'YDir','normal','XLim',[0 max(bins)*param.divDist]);
                    title(num2str([spatialInfoSess{i}.spatialInfo(j) ...
                                spatialInfoSess{i}.sparsity(j)...
                                spatialInfoSess{i}.SNR(j)]));

                end
            else
                spatialInfoSess{i}.adaptSpatialInfo(j) = nan;
                spatialInfoSess{i}.spatialInfo(j) = nan;
                spatialInfoSess{i}.sparsity(j) = nan;
                spatialInfoSess{i}.SNR(j) = nan;
            end
        end
    end
    
    fullPath = [path fileNameSpInfo];
    save(fullPath, 'spatialInfoSess');
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
