function GetFRMapInfo2P(path, fileName, onlyRun)
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
    
     %%%%%%%%% initialize constants
     GlobalConst2P;
    
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    
    fileNameInfo = [fileName '_Info.mat'];
    fileNameSpInfo = [fileName '_SpInfo_Run' num2str(onlyRun) '.mat'];
    fileName = [fileName '.mat'];
    
    fullPath = [path fileName];
    if(exist(fullPath) == 0)
        disp('The recording file does not exist');
        return;
    end
    load(fullPath);
    totClu = length(cluList.localClu);
    
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo_smTr(path,fileName);
    end
    load(fullPath);
    mazeSess = beh.mazeSessAll;
   
    param.divDist = 5; % mm
    param.smooth = 10; % mm
    param.sampleFq = sampleFq;

    drawFig = 0;
    if(drawFig == 1)
        figure;
    end
    
    spikes = getRecField2P(trials,'spikesSM',1:length(lapList));
    
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
        indLaps = setdiff(indLaps,1:startTrNo);
        
        % Collect the distance information over all the trials within a
        % subsession
        distTr = [];
        spikesTr = cell(1,totClu);
        for j = 1:length(indLaps)
            if(onlyRun == 1)
                ind = trials{indLaps(j)}.speed > minSpeed1;
            else
                ind = 1:trials{indLaps(j)}.Nsamples;
            end
            distTmp = trials{indLaps(j)}.xMM(ind);
            distTr = [distTr;distTmp];
            for n = 1:totClu
                spikesTr{n} = [spikesTr{n}; spikes{indLaps(j)}(ind,n)];
            end
        end
        
        for j = 1:totClu
            if(sum(spikesTr{j}) > 0)
                spatialInfoSess{i}.meanFR(j) = sum(spikesTr{j}) / ...
                    (length(distTr)/sampleFq);

                % Bin the spikes over space, and smooth the curve after binning
                [spatialInfoSess{i}.smoothedFR(j,:), ...
                    spatialInfoSess{i}.smoothedBinProb(j,:), ...
                    spatialInfoSess{i}.normSmoothedFR(j,:),...
                    spatialInfoSess{i}.binnedSpikes(j,:),...
                    spatialInfoSess{i}.binnedTime(j,:),...
                    spatialInfoSess{i}.binProb(j,:)] =...
                    binCoordiate2P(spikesTr{j},distTr,param);

                spatialInfoSess{i}.spatialInfo(j) = ...
                    getSpInfo12P(spatialInfoSess{i}.smoothedFR(j,:), ...
                        spatialInfoSess{i}.smoothedBinProb(j,:),...
                        spatialInfoSess{i}.meanFR(j));

                [spatialInfoSess{i}.adaptSpatialInfo(j),...
                    spatialInfoSess{i}.adaptSmoothedFR(j,:)] = ...
                    getSpInfo_addaptBin12P(spatialInfoSess{i}.binnedSpikes(j,:),...
                        spatialInfoSess{i}.binnedTime(j,:), ...
                        spatialInfoSess{i}.binProb(j,:),...
                        spatialInfoSess{i}.meanFR(j));

                [spatialInfoSess{i}.sparsity(j), spatialInfoSess{i}.SNR(j)] = ...
                    sparsityInfo2P(spatialInfoSess{i}.smoothedBinProb(j,:), ...
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
