function GetFRMapInfo_AlignedCue2P(path, fileName, mazeSess, onlyRun, intervalT)
% calculate spatial information
%
% by Yingxue, 2017.08.24
    
    if nargin<3
        disp('At least three arguments are needed for this function.');
        return;
    elseif(nargin == 3)
        onlyRun = 1;
        intervalT = 10; %sec.
    elseif(nargin == 4)
        intervalT = 10; %sec.
    elseif nargin > 5
        disp('Too many input arguments');        
        return;
    end
    
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    
    fileNameInfo = [fileName '_Info.mat'];
    fileNameRun = [fileName '_alignCue_msess' num2str(mazeSess) '.mat']; 
    fileNameSpInfo = [fileName '_SpInfoAlignedCue_msess' num2str(mazeSess) ...
                        '_Run' num2str(onlyRun) '.mat'];
    fileNamePeakFR = [fileName '_PeakFRAlignedRun_msess' num2str(mazeSess) ...
                        '_Run' num2str(onlyRun) '.mat'];
    fileName = [fileName '.mat'];
    
    fullPath = [path fileName];
    if(exist(fullPath) == 0)
        disp('The recording file does not exist');
        return;
    end
    load(fullPath,'cluList');
    totClu = length(cluList.localClu);
    
    fullPath = [path fileNameRun];
    if(exist(fullPath) == 0)
        disp(['align to run file does not exist. Please run function',...
              ' alignToCue first']);
        return;
    end
    load(fullPath);
        
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo_smTr(path,fileName);
    end
    load(fullPath);
    
    fullPath = [path fileNamePeakFR];  
    if(exist(fullPath) == 0)
        disp(['The peak firing rate file does not exist. Please call ',...
                'function "PeakFiringRate_Aligned" first.']);
        return;
    end
    load(fullPath,'trialNoNonStimGood');
    
    GlobalConst2P; 
    
    param.divDist = timeStep; % sec
    std = 0.03; % sec
    param.smooth = std*6; % sec

    drawFig = 0;
    if(drawFig == 1)
        figure;
    end
    
    spatialInfoSessNonStimGood = struct('meanFR',zeros(1,totClu),... % mean firing rate
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
    
    param.minSpeed = minSpeed;
    param.onlyRun = onlyRun;
    param.sampleFq = sampleFq;
    param.intervalT = intervalT;
    
    disp('Calculate spatial information for each subsession')
    indLaps = setdiff(trialNoNonStimGood,1:startTrNo); % changed on 2/19/2021

    % Collect the distance information over all the trials within a
    % subsession
    timeTr = [];
    spikesTr = cell(1,totClu);
    for j = 1:length(indLaps)
        if(onlyRun == 1)
            ind = find(trialsCue.speed_MMsec{indLaps(j)} > minSpeed);
        else
            ind = (1:trialsCue.numSamples(indLaps(j)))';
        end
        ind = ind(ind <= intervalT*sampleFq);
        timeTr = [timeTr; ind/sampleFq];
        for n = 1:totClu
            spikesTr{n} = [spikesTr{n}; trialsCue.spikes{indLaps(j)}(ind,n)];
        end
    end

    for j = 1:totClu
        if(sum(spikesTr{j}) > 0)
            spatialInfoSessNonStimGood.meanFR(j) = sum(spikesTr{j}) / ...
                (length(timeTr)/sampleFq);

            % Bin the spikes over space, and smooth the curve after binning
            [spatialInfoSessNonStimGood.smoothedFR(j,:), ...
                spatialInfoSessNonStimGood.smoothedBinProb(j,:), ...
                spatialInfoSessNonStimGood.normSmoothedFR(j,:),...
                spatialInfoSessNonStimGood.binnedSpikes(j,:),...
                spatialInfoSessNonStimGood.binnedTime(j,:),...
                spatialInfoSessNonStimGood.binProb(j,:)] =...
                binCoordiate2P(spikesTr{j},timeTr,param);

            spatialInfoSessNonStimGood.spatialInfo(j) = ...
                getSpInfo12P(spatialInfoSessNonStimGood.smoothedFR(j,:), ...
                    spatialInfoSessNonStimGood.smoothedBinProb(j,:),...
                    spatialInfoSessNonStimGood.meanFR(j));

            [spatialInfoSessNonStimGood.adaptSpatialInfo(j),...
                spatialInfoSessNonStimGood.adaptSmoothedFR(j,:)] = ...
                getSpInfo_addaptBin12P(spatialInfoSessNonStimGood.binnedSpikes(j,:),...
                    spatialInfoSessNonStimGood.binnedTime(j,:), ...
                    spatialInfoSessNonStimGood.binProb(j,:),...
                    spatialInfoSessNonStimGood.meanFR(j));

            [spatialInfoSessNonStimGood.sparsity(j), spatialInfoSessNonStimGood.SNR(j)] = ...
                sparsityInfo2P(spatialInfoSessNonStimGood.smoothedBinProb(j,:), ...
                spatialInfoSessNonStimGood.smoothedFR(j,:),...
                spatialInfoSessNonStimGood.meanFR(j));

            if drawFig == 1
                nl=1;
                nc=2;
                subplot(nl,nc,1); cla;
                plot(1:timeTr(1:2:end),ones(1,length(1:timeTr(1:2:end))),...
                    '.', 'MarkerSize', 3,'Color', [0.6 0.6 0.6]);
                hold on;
                plot(spikesTr{j},ones(1,length(spikesTr{j})),...
                    '.','MarkerSize', 3,'Color', [1 0 0]);
                xlim([0 max(timeTr)]);

                title(['totClu:' num2str(j) '     sh:' num2str(cluList.shank(j))...
                    '  locClu:' num2str(cluList.localClu(j))]);

                % place fields
                subplot(nl,nc,2); cla;
                trBin = round(timeTr/param.divDist)+1;
                bins = unique(trBin);
                maxC = max(spatialInfoSessNonStimGood.smoothedFR(j,:));
                imagesc(bins'*param.divDist, ones(1,length(bins)),...
                    spatialInfoSessNonStimGood.smoothedFR(j,:),[0 maxC]); 
                hold on;
                set(gca,'YDir','normal','XLim',[0 max(bins)*param.divDist]);
                title(num2str([spatialInfoSessNonStimGood.spatialInfo(j) ...
                            spatialInfoSessNonStimGood.sparsity(j)...
                            spatialInfoSessNonStimGood.SNR(j)]));

            end
        else
            spatialInfoSessNonStimGood.adaptSpatialInfo(j) = nan;
            spatialInfoSessNonStimGood.spatialInfo(j) = nan;
            spatialInfoSessNonStimGood.sparsity(j) = nan;
            spatialInfoSessNonStimGood.SNR(j) = nan;
        end
    end
    
    fullPath = [path fileNameSpInfo];
    save(fullPath, 'spatialInfoSessNonStimGood');
end
