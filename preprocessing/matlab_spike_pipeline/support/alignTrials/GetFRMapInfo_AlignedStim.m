function GetFRMapInfo_AlignedStim(path, fileName, mazeSess, onlyRun, intervalT)
% calculate spatial information
%
% by Yingxue, 2021.05.21
    
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
    fileNameRun = [fileName '_alignRun_msess' num2str(mazeSess) '.mat']; 
    fileNameSpInfo = [fileName '_SpInfoAlignedStim_msess' num2str(mazeSess) ...
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
              ' alignToRunOnset first']);
        return;
    end
    load(fullPath);
    
    fullPath = [path fileNameSpike];
    if(exist(fullPath) == 0)
        disp(['spikes perNperT file does not exist. Please run function',...
              ' getSpikesPerNeuPerTr first']);
        return;
    end
    load(fullPath,'trialsRunSpikes');
    
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
    load(fullPath,'trialNoStim','trialNoStimCtrl','pulseMeth');
    
    GlobalConst; 
    
    param.divDist = 0.0008; % sec
    std = 0.03; % sec
    param.smooth = std*6; % sec
    param.sampleFq = sampleFq;
    param.intervalT = intervalT;
    param.minSpeed = minSpeed;

    drawFig = 0;
    if(drawFig == 1)
        figure;
    end
    
    for i = 1:length(pulseMeth)
        disp('Calculate spatial information for stimulation trials')
        indLaps = trialNoStim{i}; 
        spatialInfoSessStim{i} = calSpatioInfo(trialsRun,trialsRunSpikes,indLaps,totClu,param);
        
        disp('Calculate spatial information for stimulation control trials')
        indLaps = trialNoStimCtrl{i}; 
        spatialInfoSessStimCtrl{i} = calSpatioInfo(trialsRun,trialsRunSpikes,indLaps,totClu,param);
    end
    
    fullPath = [path fileNameSpInfo];
    save(fullPath, 'spatialInfoSessStim','spatialInfoSessStimCtrl');
    
end
    
function spatialInfoSess = calSpatioInfo(trialsRun,trialsRunSpikes,indLaps,totClu,param)

    spatialInfoSess = struct('meanFR',zeros(1,totClu),... % mean firing rate
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
                     
    % Collect the distance information over all the trials within a
    % subsession
    timeTr = [];
    for j = 1:length(indLaps)
        if(onlyRun == 1)
            ind = find(trialsRun.speed_MMsec{indLaps(j)} > param.minSpeed);
        else
            ind = (1:trialsRun.numSamples(indLaps(j)))';
        end
        ind = ind(ind <= param.intervalT*param.sampleFq);
        timeTr = [timeTr; ind/param.sampleFq];
    end

    % Collect the distance information over all the trials when spikes
    % occurs
    timeSpikes = ...
        concatenateSpikes(trialsRunSpikes.Time,indLaps,totClu,...
                            param.sampleFq,param.intervalT);

    for j = 1:totClu
        if(~isempty(timeSpikes{j}))
            spatialInfoSess.meanFR(j) = length(timeSpikes{j}) / ...
                (length(timeTr)/param.sampleFq);

            % Bin the spikes over space, and smooth the curve after binning
            [spatialInfoSess.smoothedFR(j,:), ...
                spatialInfoSess.smoothedBinProb(j,:), ...
                spatialInfoSess.normSmoothedFR(j,:),...
                spatialInfoSess.binnedSpikes(j,:),...
                spatialInfoSess.binnedTime(j,:),...
                spatialInfoSess.binProb(j,:)] =...
                binCoordiate(timeSpikes{j},timeTr,param);

            spatialInfoSess.spatialInfo(j) = ...
                getSpInfo1(spatialInfoSess.smoothedFR(j,:), ...
                    spatialInfoSess.smoothedBinProb(j,:),...
                    spatialInfoSess.meanFR(j));

            [spatialInfoSess.adaptSpatialInfo(j),...
                spatialInfoSess.adaptSmoothedFR(j,:)] = ...
                getSpInfo_addaptBin1(spatialInfoSess.binnedSpikes(j,:),...
                    spatialInfoSess.binnedTime(j,:), ...
                    spatialInfoSess.binProb(j,:),...
                    spatialInfoSess.meanFR(j));

            [spatialInfoSess.sparsity(j), spatialInfoSess.SNR(j)] = ...
                sparsityInfo(spatialInfoSess.smoothedBinProb(j,:), ...
                spatialInfoSess.smoothedFR(j,:),...
                spatialInfoSess.meanFR(j));

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
                maxC = max(spatialInfoSess.smoothedFR(j,:));
                imagesc(bins'*param.divDist, ones(1,length(bins)),...
                    spatialInfoSess.smoothedFR(j,:),[0 maxC]); 
                hold on;
                set(gca,'YDir','normal','XLim',[0 max(bins)*param.divDist]);
                title(num2str([spatialInfoSess.spatialInfo(j) ...
                            spatialInfoSess.sparsity(j)...
                            spatialInfoSess.SNR(j)]));
            end
        else
            spatialInfoSess.adaptSpatialInfo(j) = nan;
            spatialInfoSess.spatialInfo(j) = nan;
            spatialInfoSess.sparsity(j) = nan;
            spatialInfoSess.SNR(j) = nan;
        end
    end 
    
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
