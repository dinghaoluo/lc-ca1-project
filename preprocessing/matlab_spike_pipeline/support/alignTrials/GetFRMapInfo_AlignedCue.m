function GetFRMapInfo_AlignedCue(path, fileName, mazeSess, onlyRun, intervalT)
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
    fileNameSpike = [fileName '_alignedSpikesPerNPerT_msess' num2str(mazeSess)...
        '_Run' num2str(onlyRun) '.mat'];
    fileNameConv = [fileName '_convSpikesAligned_msess' num2str(mazeSess)...
        '_Run' num2str(onlyRun) '.mat'];
    fileNamePeakFR = [fileName '_PeakFRAligned_msess' num2str(mazeSess) ...
                        '_Run' num2str(onlyRun) '.mat'];
    fileName = [fileName '.mat'];
    
    fullPath = [path fileName];
    if(exist(fullPath) == 0)
        disp('The recording file does not exist');
        return;
    end
    load(fullPath,'cluList');
    totClu = length(cluList.all);
    
    fullPath = [path fileNameRun];
    if(exist(fullPath) == 0)
        disp(['align to run file does not exist. Please run function',...
              ' alignToCue first']);
        return;
    end
    load(fullPath);
    
    fullPath = [path fileNameSpike];
    if(exist(fullPath) == 0)
        disp(['spikes perNperT file does not exist. Please run function',...
              ' getSpikesPerNeuPerTr first']);
        return;
    end
    load(fullPath,'trialsCueSpikes');
    
    fullPath = [path fileNameConv];
    if(exist(fullPath) == 0)
        disp(['convSpikeAligned file does not exist. Please run function',...
              ' ConvSpikeTrain_Aligned first']);
        return;
    end
    load(fullPath,'paramC');
    
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
    
    GlobalConst; 
    
    param.divDist = 0.0008; % sec
    std = 0.03; % sec
    param.smooth = std*6; % sec
    param.sampleFq = sampleFq;

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
    
    disp('Calculate spatial information for each subsession')
    indLaps = trialNoNonStimGood; % changed on 2/19/2021

    % Collect the distance information over all the trials within a
    % subsession
    timeTr = [];
    for j = 1:length(indLaps)
        if(onlyRun == 1)
            ind = find(trialsCue.speed_MMsec{indLaps(j)} > minSpeed);
        else
            ind = (1:trialsCue.numSamples(indLaps(j)))';
        end
        ind = ind(ind <= intervalT*sampleFq);
        timeTr = [timeTr; ind/sampleFq];
    end

    % Collect the distance information over all the trials when spikes
    % occurs
    timeSpikes = ...
        concatenateSpikes(trialsCueSpikes.Time,indLaps,totClu,sampleFq,intervalT);

    for j = 1:totClu
        if(~isempty(timeSpikes{j}))
            spatialInfoSessNonStimGood.meanFR(j) = length(timeSpikes{j}) / ...
                (length(timeTr)/sampleFq);

            % Bin the spikes over space, and smooth the curve after binning
            [spatialInfoSessNonStimGood.smoothedFR(j,:), ...
                spatialInfoSessNonStimGood.smoothedBinProb(j,:), ...
                spatialInfoSessNonStimGood.normSmoothedFR(j,:),...
                spatialInfoSessNonStimGood.binnedSpikes(j,:),...
                spatialInfoSessNonStimGood.binnedTime(j,:),...
                spatialInfoSessNonStimGood.binProb(j,:)] =...
                binCoordiate(timeSpikes{j},timeTr,param);

            spatialInfoSessNonStimGood.spatialInfo(j) = ...
                getSpInfo1(spatialInfoSessNonStimGood.smoothedFR(j,:), ...
                    spatialInfoSessNonStimGood.smoothedBinProb(j,:),...
                    spatialInfoSessNonStimGood.meanFR(j));

            [spatialInfoSessNonStimGood.adaptSpatialInfo(j),...
                spatialInfoSessNonStimGood.adaptSmoothedFR(j,:)] = ...
                getSpInfo_addaptBin1(spatialInfoSessNonStimGood.binnedSpikes(j,:),...
                    spatialInfoSessNonStimGood.binnedTime(j,:), ...
                    spatialInfoSessNonStimGood.binProb(j,:),...
                    spatialInfoSessNonStimGood.meanFR(j));

            [spatialInfoSessNonStimGood.sparsity(j), spatialInfoSessNonStimGood.SNR(j)] = ...
                sparsityInfo(spatialInfoSessNonStimGood.smoothedBinProb(j,:), ...
                spatialInfoSessNonStimGood.smoothedFR(j,:),...
                spatialInfoSessNonStimGood.meanFR(j));

            if drawFig == 1
                nl=1;
                nc=2;
                subplot(nl,nc,1); cla;
                plot(1:timeTr(1:2:end),ones(1,length(1:timeTr(1:2:end))),...
                    '.', 'MarkerSize', 3,'Color', [0.6 0.6 0.6]);
                hold on;
                plot(timeSpikes{j},ones(1,length(timeSpikes{j})),...
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

function distSpikes = ...
        concatenateSpikes(distSp,indLaps,totClu,sampleFq,intervalT)
    
    distSpikes = cell(1,totClu);
    
    for i = 1:length(indLaps)
        for j = 1:totClu
            ind = distSp{j,indLaps(i)} <= intervalT*sampleFq;
%             if(sum(ind) < length(distSp{j,indLaps(i)}))
%                 a=1;
%             end
            distSpikes{j} = [distSpikes{j};distSp{j,indLaps(i)}(ind)/sampleFq];
        end
    end
end
