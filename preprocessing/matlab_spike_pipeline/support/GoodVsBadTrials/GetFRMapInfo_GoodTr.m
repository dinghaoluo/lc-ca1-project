function GetFRMapInfo_GoodTr(path, fileName, onlyRun)
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
    fileNameSpInfo = [fileName '_SpInfo_GoodTr_Run' num2str(onlyRun) '.mat'];
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
    spatialInfoSessGoodTr = cell(1,length(mazeSess));
    spatialInfoSessBadTr = cell(1,length(mazeSess));
    spatialInfoSessOKTr = cell(1,length(mazeSess));
    for i = 1:length(mazeSess)    
        disp(['Session ' num2str(i)]);
        spatialInfoSessGoodTr{i} = spatialInfo;
        spatialInfoSessBadTr{i} = spatialInfo;
        spatialInfoSessOKTr{i} = spatialInfo;
        
        if(length(mazeSess) > 1)
            indLaps = find(beh.mazeSess == mazeSess(i)); 
            indLaps = intersect(indLaps, beh.indTrCtrl);
        else
            indLaps = beh.indTrCtrl;
        end
        indLapsGoodTr = intersect(indLaps,beh.indGoodTrCtrl); 
        indLapsBadTr = intersect(indLaps,beh.indBadTrCtrl); 
        indLapsOKTr = setdiff(indLaps,[beh.indGoodTrCtrl,beh.indBadTrCtrl]);
        
        % Collect the distance information over all the trials within a
        % subsession
        distTrGoodTr = [];
        for j = 1:length(indLapsGoodTr)
            if(onlyRun == 1)
                ind = trials{indLapsGoodTr(j)}.speed > minSpeed;
            else
                ind = 1:trials{indLapsGoodTr(j)}.Nsamples;
            end
            distTmp = trials{indLapsGoodTr(j)}.xMM(ind);
            distTrGoodTr = [distTrGoodTr;distTmp];
        end
        
        distTrBadTr = [];
        for j = 1:length(indLapsBadTr)
            if(onlyRun == 1)
                ind = trials{indLapsBadTr(j)}.speed > minSpeed;
            else
                ind = 1:trials{indLapsBadTr(j)}.Nsamples;
            end
            distTmp = trials{indLapsBadTr(j)}.xMM(ind);
            distTrBadTr = [distTrBadTr;distTmp];
        end
        
        distTrOKTr = [];
        for j = 1:length(indLapsOKTr)
            if(onlyRun == 1)
                ind = trials{indLapsOKTr(j)}.speed > minSpeed;
            else
                ind = 1:trials{indLapsOKTr(j)}.Nsamples;
            end
            distTmp = trials{indLapsOKTr(j)}.xMM(ind);
            distTrOKTr = [distTrOKTr;distTmp];
        end
        
        % Collect the distance information over all the trials when spikes
        % occurs
        distSpikesGoodTr = ...
            concatenateSpikes(distSp,indLapsGoodTr,totClu);
        distSpikesBadTr = ...
            concatenateSpikes(distSp,indLapsBadTr,totClu);
        distSpikesOKTr = ...
            concatenateSpikes(distSp,indLapsOKTr,totClu);        
        
        for j = 1:totClu
            if(~isempty(distSpikesGoodTr{j}))
                spatialInfoSessGoodTr{i}.meanFR(j) = length(distSpikesGoodTr{j}) / ...
                    (length(distTrGoodTr)/sampleFq);

                % Bin the spikes over space, and smooth the curve after binning
                [spatialInfoSessGoodTr{i}.smoothedFR(j,:), ...
                    spatialInfoSessGoodTr{i}.smoothedBinProb(j,:), ...
                    spatialInfoSessGoodTr{i}.normSmoothedFR(j,:),...
                    spatialInfoSessGoodTr{i}.binnedSpikes(j,:),...
                    spatialInfoSessGoodTr{i}.binnedTime(j,:),...
                    spatialInfoSessGoodTr{i}.binProb(j,:)] =...
                    binCoordiate(distSpikesGoodTr{j},distTrGoodTr,param);

                spatialInfoSessGoodTr{i}.spatialInfo(j) = ...
                    getSpInfo1(spatialInfoSessGoodTr{i}.smoothedFR(j,:), ...
                        spatialInfoSessGoodTr{i}.smoothedBinProb(j,:),...
                        spatialInfoSessGoodTr{i}.meanFR(j));

                [spatialInfoSessGoodTr{i}.adaptSpatialInfo(j),...
                    spatialInfoSessGoodTr{i}.adaptSmoothedFR(j,:)] = ...
                    getSpInfo_addaptBin1(spatialInfoSessGoodTr{i}.binnedSpikes(j,:),...
                        spatialInfoSessGoodTr{i}.binnedTime(j,:), ...
                        spatialInfoSessGoodTr{i}.binProb(j,:),...
                        spatialInfoSessGoodTr{i}.meanFR(j));

                [spatialInfoSessGoodTr{i}.sparsity(j), spatialInfoSessGoodTr{i}.SNR(j)] = ...
                    sparsityInfo(spatialInfoSessGoodTr{i}.smoothedBinProb(j,:), ...
                    spatialInfoSessGoodTr{i}.smoothedFR(j,:),...
                    spatialInfoSessGoodTr{i}.meanFR(j));

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
                    maxC = max(spatialInfoSessGoodTr{i}.smoothedFR(j,:));
                    imagesc(bins'*param.divDist, ones(1,length(bins)),...
                        spatialInfoSessGoodTr{i}.smoothedFR(j,:),[0 maxC]); 
                    hold on;
                    set(gca,'YDir','normal','XLim',[0 max(bins)*param.divDist]);
                    title(num2str([spatialInfoSessGoodTr{i}.spatialInfo(j) ...
                                spatialInfoSessGoodTr{i}.sparsity(j)...
                                spatialInfoSessGoodTr{i}.SNR(j)]));

                end
            else
                spatialInfoSessGoodTr{i}.adaptSpatialInfo(j) = nan;
                spatialInfoSessGoodTr{i}.spatialInfo(j) = nan;
                spatialInfoSessGoodTr{i}.sparsity(j) = nan;
                spatialInfoSessGoodTr{i}.SNR(j) = nan;
            end
            
            if(~isempty(distSpikesBadTr{j}))
                spatialInfoSessBadTr{i}.meanFR(j) = length(distSpikesBadTr{j}) / ...
                    (length(distTrBadTr)/sampleFq);

                % Bin the spikes over space, and smooth the curve after binning
                [spatialInfoSessBadTr{i}.smoothedFR(j,:), ...
                    spatialInfoSessBadTr{i}.smoothedBinProb(j,:), ...
                    spatialInfoSessBadTr{i}.normSmoothedFR(j,:),...
                    spatialInfoSessBadTr{i}.binnedSpikes(j,:),...
                    spatialInfoSessBadTr{i}.binnedTime(j,:),...
                    spatialInfoSessBadTr{i}.binProb(j,:)] =...
                    binCoordiate(distSpikesBadTr{j},distTrBadTr,param);

                spatialInfoSessBadTr{i}.spatialInfo(j) = ...
                    getSpInfo1(spatialInfoSessBadTr{i}.smoothedFR(j,:), ...
                        spatialInfoSessBadTr{i}.smoothedBinProb(j,:),...
                        spatialInfoSessBadTr{i}.meanFR(j));

                [spatialInfoSessBadTr{i}.adaptSpatialInfo(j),...
                    spatialInfoSessBadTr{i}.adaptSmoothedFR(j,:)] = ...
                    getSpInfo_addaptBin1(spatialInfoSessBadTr{i}.binnedSpikes(j,:),...
                        spatialInfoSessBadTr{i}.binnedTime(j,:), ...
                        spatialInfoSessBadTr{i}.binProb(j,:),...
                        spatialInfoSessBadTr{i}.meanFR(j));

                [spatialInfoSessBadTr{i}.sparsity(j), spatialInfoSessBadTr{i}.SNR(j)] = ...
                    sparsityInfo(spatialInfoSessBadTr{i}.smoothedBinProb(j,:), ...
                    spatialInfoSessBadTr{i}.smoothedFR(j,:),...
                    spatialInfoSessBadTr{i}.meanFR(j));
            else
                spatialInfoSessBadTr{i}.adaptSpatialInfo(j) = nan;
                spatialInfoSessBadTr{i}.spatialInfo(j) = nan;
                spatialInfoSessBadTr{i}.sparsity(j) = nan;
                spatialInfoSessBadTr{i}.SNR(j) = nan;
            end
            
            if(~isempty(distSpikesOKTr{j}))
                spatialInfoSessOKTr{i}.meanFR(j) = length(distSpikesOKTr{j}) / ...
                    (length(distTrOKTr)/sampleFq);

                % Bin the spikes over space, and smooth the curve after binning
                [spatialInfoSessOKTr{i}.smoothedFR(j,:), ...
                    spatialInfoSessOKTr{i}.smoothedBinProb(j,:), ...
                    spatialInfoSessOKTr{i}.normSmoothedFR(j,:),...
                    spatialInfoSessOKTr{i}.binnedSpikes(j,:),...
                    spatialInfoSessOKTr{i}.binnedTime(j,:),...
                    spatialInfoSessOKTr{i}.binProb(j,:)] =...
                    binCoordiate(distSpikesOKTr{j},distTrOKTr,param);

                spatialInfoSessOKTr{i}.spatialInfo(j) = ...
                    getSpInfo1(spatialInfoSessOKTr{i}.smoothedFR(j,:), ...
                        spatialInfoSessOKTr{i}.smoothedBinProb(j,:),...
                        spatialInfoSessOKTr{i}.meanFR(j));

                [spatialInfoSessOKTr{i}.adaptSpatialInfo(j),...
                    spatialInfoSessOKTr{i}.adaptSmoothedFR(j,:)] = ...
                    getSpInfo_addaptBin1(spatialInfoSessOKTr{i}.binnedSpikes(j,:),...
                        spatialInfoSessOKTr{i}.binnedTime(j,:), ...
                        spatialInfoSessOKTr{i}.binProb(j,:),...
                        spatialInfoSessOKTr{i}.meanFR(j));

                [spatialInfoSessOKTr{i}.sparsity(j), spatialInfoSessOKTr{i}.SNR(j)] = ...
                    sparsityInfo(spatialInfoSessOKTr{i}.smoothedBinProb(j,:), ...
                    spatialInfoSessOKTr{i}.smoothedFR(j,:),...
                    spatialInfoSessOKTr{i}.meanFR(j));
            else
                spatialInfoSessOKTr{i}.adaptSpatialInfo(j) = nan;
                spatialInfoSessOKTr{i}.spatialInfo(j) = nan;
                spatialInfoSessOKTr{i}.sparsity(j) = nan;
                spatialInfoSessOKTr{i}.SNR(j) = nan;
            end
        end
    end
    
    fullPath = [path fileNameSpInfo];
    save(fullPath, 'spatialInfoSessGoodTr','spatialInfoSessBadTr',...
        'spatialInfoSessOKTr');
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
