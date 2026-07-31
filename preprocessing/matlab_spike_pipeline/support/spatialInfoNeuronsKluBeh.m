function spatialInfoNeuronsKluBeh(path,fileName,onlyRun)
% calculate the spatial information for the good and bad behavior trials
    
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
    fileNameSpInfo = [fileName '_SpInfoKluBeh_Run' num2str(onlyRun) '.mat'];
    fileNameCluster = [fileName '_kmeanBeh.mat'];
    fileName = [fileName '.mat'];
    
    fullPath = [path fileName];
    if(exist(fullPath) == 0)
        disp('The recording file does not exist');
        return;
    end
    load(fullPath);
        
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
        BasicInfo(path,fileName);
    end
    load(fullPath);
    mazeSess = beh.mazeSessAll;
    
    fullPath = [path fileNameCluster];
    if(exist(fullPath) == 0)
        disp('The behavior cluster file does not exist.');
        disp(fullPath)
        return;
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
    
    indPyrNeu = find(autoCorr.isPyrneuron == 1);
    totClu = length(indPyrNeu);
%     ind = abs(cluList.centerMax) > 500;

    [indCluGood,indCluBad] = HandSelectedCluFromBehKMean(fileName);
    
    indTrGood = indLaps(kmeanClu == indCluGood);
    indTrBad = indLaps(kmeanClu == indCluBad);
    
    spatialInfoKluBehGood = struct('indPyrNeu', indPyrNeu, ...
                         'indLaps', indTrGood,...
                         'meanFR',zeros(1,totClu),... % mean firing rate
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
                     
     spatialInfoKluBehBad = struct('indPyrNeu', indPyrNeu, ...
                         'indLaps', indTrBad,...
                         'meanFR',zeros(1,totClu),... % mean firing rate
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
    
    disp('Calculate spatial information for good behavior trials')
    % Collect the distance information over all the trials within the good
    % behavior kluster
    distTrGood = [];
    for j = 1:length(indTrGood)
        if(onlyRun == 1)
            ind = trials{indTrGood(j)}.speed > minSpeed;
        else
            ind = 1:trials{indTrGood(j)}.Nsamples;
        end
        distTmp = trials{indTrGood(j)}.xMM(ind);
        distTrGood = [distTrGood;distTmp];
    end
    
    % Collect the distance information over all the trials when spikes
    % occurs
    distSpikesGood = ...
        concatenateSpikes(distSp,indTrGood,indPyrNeu);
        
    spatialInfoKluBehGood = ...
            calSpatialInfo(distSpikesGood,distTrGood,param,totClu,sampleFq,0);
    
    disp('Calculate spatial information for bad behavior trials')
    % Collect the distance information over all the trials within the bad
    % behavior kluster
    distTrBad = [];
    for j = 1:length(indTrBad)
        if(onlyRun == 1)
            ind = trials{indTrBad(j)}.speed > minSpeed;
        else
            ind = 1:trials{indTrBad(j)}.Nsamples;
        end
        distTmp = trials{indTrBad(j)}.xMM(ind);
        distTrBad = [distTrBad;distTmp];
    end
    
    % Collect the distance information over all the trials when spikes
    % occurs
    distSpikesBad = ...
        concatenateSpikes(distSp,indTrBad,indPyrNeu);    
     
    spatialInfoKluBehBad = ...
        calSpatialInfo(distSpikesBad,distTrBad,param,totClu,sampleFq,0);
    
    indSelPyr = [];
    if(~isempty(strfind(fileName,'A011-20190218')))
        indSelPyr = [3 21 53 54 59 80 95 129 134 140 144 165 172 177 186 206 216 226];
    elseif(~isempty(strfind(fileName,'A011-20190219')))
        indSelPyr = [9 16 19 90 93 130 140 143 182];
    elseif(~isempty(strfind(fileName,'A012-20190224')))
        indSelPyr = [13 15 16 30 31 38 39 55 63 65 66 85 86 87 88 98 103 108 134];
    elseif(~isempty(strfind(fileName,'A012-20190221')))
        indSelPyr = [15 30 31 33 35 40 42 43 52 66 69 70];
    end
    
    if(~isempty(indSelPyr))
        [~,ia,ib] = intersect(indPyrNeu,indSelPyr);
        figure;
        plot(spatialInfoKluBehGood.spatialInfo(ia),...
                spatialInfoKluBehBad.spatialInfo(ia),'bo');
        hold on;
        plot(spatialInfoKluBehGood.sparsity(ia),...
                spatialInfoKluBehBad.sparsity(ia),'go');
        maxSparsity =  max([spatialInfoKluBehGood.sparsity(ia),...
                spatialInfoKluBehBad.sparsity(ia)]);
        plot([0,maxSparsity],[0, maxSparsity],'r');
        set(gca,'XLim',[0 maxSparsity]);
        title(fileName(1:13))
        xlabel('Good behavior')
        ylabel('Bad behavior')
    end
    
    fullPath = [path fileNameSpInfo];
    save(fullPath, 'spatialInfoKluBehGood','spatialInfoKluBehBad');
end

function distSpikes = ...
        concatenateSpikes(distSp,indLaps,indPyrNeu)
    
    distSpikes = cell(1,length(indPyrNeu));
    
    for i = 1:length(indLaps)
        for j = 1:length(indPyrNeu)
            distSpikes{j} = [distSpikes{j};distSp{indLaps(i)}{indPyrNeu(j)}];
        end
    end
end

function spatialInfo = calSpatialInfo(distSpikes,distTr,param,totClu,sampleFq,drawFig)
    for j = 1:totClu
        if(~isempty(distSpikes{j}))
            spatialInfo.meanFR(j) = length(distSpikes{j}) / ...
                (length(distTr)/sampleFq);

            % Bin the spikes over space, and smooth the curve after binning
            [spatialInfo.smoothedFR(j,:), ...
                spatialInfo.smoothedBinProb(j,:), ...
                spatialInfo.normSmoothedFR(j,:),...
                spatialInfo.binnedSpikes(j,:),...
                spatialInfo.binnedTime(j,:),...
                spatialInfo.binProb(j,:)] =...
                binCoordiate(distSpikes{j},distTr,param);

            spatialInfo.spatialInfo(j) = ...
                getSpInfo1(spatialInfo.smoothedFR(j,:), ...
                    spatialInfo.smoothedBinProb(j,:),...
                    spatialInfo.meanFR(j));

            [spatialInfo.adaptSpatialInfo(j),...
                spatialInfo.adaptSmoothedFR(j,:)] = ...
                getSpInfo_addaptBin1(spatialInfo.binnedSpikes(j,:),...
                    spatialInfo.binnedTime(j,:), ...
                    spatialInfo.binProb(j,:),...
                    spatialInfo.meanFR(j));

            [spatialInfo.sparsity(j), spatialInfo.SNR(j)] = ...
                sparsityInfo(spatialInfo.smoothedBinProb(j,:), ...
                spatialInfo.smoothedFR(j,:),...
                spatialInfo.meanFR(j));

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
                maxC = max(spatialInfo.smoothedFR(j,:));
                imagesc(bins'*param.divDist, ones(1,length(bins)),...
                    spatialInfo.smoothedFR(j,:),[0 maxC]); 
                hold on;
                set(gca,'YDir','normal','XLim',[0 max(bins)*param.divDist]);
                title(num2str([spatialInfo.spatialInfo(j) ...
                            spatialInfo.sparsity(j)...
                            spatialInfo.SNR(j)]));

            end
        else
            spatialInfo.adaptSpatialInfo(j) = nan;
            spatialInfo.spatialInfo(j) = nan;
            spatialInfo.sparsity(j) = nan;
            spatialInfo.SNR(j) = nan;
        end
    end
end
