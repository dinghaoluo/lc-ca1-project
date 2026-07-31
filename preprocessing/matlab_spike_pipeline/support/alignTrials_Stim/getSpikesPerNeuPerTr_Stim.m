function getSpikesPerNeuPerTr_Stim(path,fileName,onlyRun,mazeSess,isStim)
% separate the spikes for each neuron from each trial
% e.g. getSpikesPerNeuPerTr('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1,1)
    %%%%%%%%% load recording file
    
    %% added by Yingxue on 2/3/2021
    if(nargin == 4)
        isStim = 0;
    end
    %%
    
    fullPath = [path fileName '.mat']; 
    if(exist(fullPath) == 0)
        disp('The recording file does not exist');
        return;
    end
    load(fullPath,'cluList','lap');
    neuronNo = 1:length(cluList.all);
    
    fullPath = [path fileName '_alignRun_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to run file does not exist');
        return;
    end
    load(fullPath,'trialsRun');
    
    fullPath = [path fileName '_alignRew_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to reward file does not exist');
        return;
    end
    load(fullPath,'trialsRew');
    
    fullPath = [path fileName '_alignCue_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to cue file does not exist');
        return;
    end
    load(fullPath,'trialsCue');
    
    fullPath = [path fileName '_alignCueOff_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to cue file does not exist');
        return;
    end
    load(fullPath,'trialsCueOff');
    
    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
        
    GlobalConst;
    %% changed by Yingxue on 2/3/2021
    if(isStim == 1)
        trialNo = behPar.indLapsStimInSess;
    else
        trialNo = setdiff(behPar.indLapsSess,behPar.indLapsStimSess);
    end
    %%
    trialNoCueOff = intersect(trialNo,find(behPar.numSamplesCueOff ~= 0));
    trialsCueOffSpikes.Time = [];
    trialsCueOffSpikes.thPhaseInterpSpike = [];
    trialsCueOffSpikes.TimeBef = [];
    trialsCueOffSpikes.thPhaseInterpSpikeBef = [];
    trialsCueOffSpikes.MaxLenT = 0;
    trialsCueOffSpikes.MeanLenT = 0;
    trialsCueOffSpikes.Dist = 0;
    trialsCueOffSpikes.MaxLenDist = 0;
    trialsCueOffSpikes.MeanLenDist = 0;
    if(onlyRun == 0)
        trialsRunSpikes.Time = perNeuPerTr(trialsRun.res,trialsRun.totclu,neuronNo,trialNo);
        trialsRewSpikes.Time = perNeuPerTr(trialsRew.res,trialsRew.totclu,neuronNo,trialNo);
        trialsCueSpikes.Time = perNeuPerTr(trialsCue.res,trialsCue.totclu,neuronNo,trialNo);

        %%%%%% added by Yingxue on 3/11/2020
        trialsRunSpikes.Time_LasttoCurTr = perNeuPerTr(trialsRun.res_LasttoCurTr,trialsRun.totclu_LasttoCurTr,neuronNo,trialNo);
        trialsRewSpikes.Time_LasttoCurTr = perNeuPerTr(trialsRew.res_LasttoCurTr,trialsRew.totclu_LasttoCurTr,neuronNo,trialNo);
        trialsCueSpikes.Time_LasttoCurTr = perNeuPerTr(trialsCue.res_LasttoCurTr,trialsCue.totclu_LasttoCurTr,neuronNo,trialNo);
        %%%%%%

        if(~isempty(trialNoCueOff))
            trialsCueOffSpikes.Time = perNeuPerTr(trialsCueOff.res,trialsCueOff.totclu,neuronNo,trialNoCueOff);
            %%%%%% added by Yingxue on 3/11/2020
            trialsCueOffSpikes.Time_LasttoCurTr = perNeuPerTr(trialsCueOff.res_LasttoCurTr,trialsCueOff.totclu_LasttoCurTr,neuronNo,trialNoCueOff);
            %%%%%% 
        end
        
        trialsRunSpikes.thPhaseInterpSpike = perNeuPerTr(trialsRun.thPhaseInterpSpike,trialsRun.totclu,neuronNo,trialNo);
        
        % added on 4/27/2020
        trialsRunSpikes.thPhaseHilbSpike = perNeuPerTr(trialsRun.thPhaseHilbSpike,trialsRun.totclu,neuronNo,trialNo);

        trialsRunSpikes.TimeBef = perNeuPerTr(trialsRun.resBef,trialsRun.totcluBef,neuronNo,trialNo);
        trialsRunSpikes.thPhaseInterpSpikeBef = perNeuPerTr(trialsRun.thPhaseInterpSpikeBef,trialsRun.totcluBef,neuronNo,trialNo);
        
        if(~isempty(trialNoCueOff))
            trialsCueOffSpikes.thPhaseInterpSpike = perNeuPerTr(trialsCueOff.thPhaseInterpSpike,trialsCueOff.totclu,neuronNo,trialNoCueOff);

            trialsCueOffSpikes.TimeBef = perNeuPerTr(trialsCueOff.resBef,trialsCueOff.totcluBef,neuronNo,trialNoCueOff);
            trialsCueOffSpikes.thPhaseInterpSpikeBef = perNeuPerTr(trialsCueOff.thPhaseInterpSpikeBef,trialsCueOff.totcluBef,neuronNo,trialNoCueOff);
            
            trialsCueOffSpikes.MaxLenT = max(trialsCueOff.numSamples(trialNoCueOff));
            
            trialsCueOffSpikes.MeanLenT = mean(trialsCueOff.numSamples(trialNoCueOff));
            
            trialsCueOffSpikes.Dist = perNeuPerTr(trialsCueOff.xMMSpike,trialsCueOff.totclu,neuronNo,trialNoCueOff);
            
            [trialsCueOffSpikes.MaxLenDist,trialsCueOffSpikes.MeanLenDist] = maxDistance(trialsCueOff.xMM,trialNoCueOff);  
        end

        trialsRunSpikes.MaxLenT = max(trialsRun.numSamples);
        trialsRewSpikes.MaxLenT = max(trialsRew.numSamples);
        trialsCueSpikes.MaxLenT = max(trialsCue.numSamples);
        
        trialsRunSpikes.MeanLenT = mean(trialsRun.numSamples);
        trialsRewSpikes.MeanLenT = mean(trialsRew.numSamples);
        trialsCueSpikes.MeanLenT = mean(trialsCue.numSamples);
        
        trialsRunSpikes.Dist = perNeuPerTr(trialsRun.xMMSpike,trialsRun.totclu,neuronNo,trialNo);
        trialsRewSpikes.Dist = perNeuPerTr(trialsRew.xMMSpike,trialsRew.totclu,neuronNo,trialNo);
        trialsCueSpikes.Dist = perNeuPerTr(trialsCue.xMMSpike,trialsCue.totclu,neuronNo,trialNo);
        
        [trialsRunSpikes.MaxLenDist,trialsRunSpikes.MeanLenDist] = maxDistance(trialsRun.xMM,trialNo);
        [trialsRewSpikes.MaxLenDist,trialsRewSpikes.MeanLenDist] = maxDistance(trialsRew.xMM,trialNo);
        [trialsCueSpikes.MaxLenDist,trialsCueSpikes.MeanLenDist] = maxDistance(trialsCue.xMM,trialNo);  
        
    else
        trialsRunSpikes.Time = ...
            perNeuPerTrRun(trialsRun.res,trialsRun.totclu,trialsRun.speed_MMsecSpike,...
                           neuronNo,trialNo,minSpeed);
        trialsRewSpikes.Time = ...
            perNeuPerTrRun(trialsRew.res,trialsRew.totclu,trialsRew.speed_MMsecSpike,...
                           neuronNo,trialNo,minSpeed);
        trialsCueSpikes.Time = ...
            perNeuPerTrRun(trialsCue.res,trialsCue.totclu,trialsCue.speed_MMsecSpike,...
                           neuronNo,trialNo,minSpeed);

        %%%%%% added by Yingxue on 3/11/2020
        trialsRunSpikes.Time_LasttoCurTr = ...
            perNeuPerTrRun(trialsRun.res_LasttoCurTr,trialsRun.totclu_LasttoCurTr,...
                        trialsRun.speed_MMsecSpike_LasttoCurTr,neuronNo,trialNo,minSpeed);
        trialsRewSpikes.Time_LasttoCurTr = ...
            perNeuPerTrRun(trialsRew.res_LasttoCurTr,trialsRew.totclu_LasttoCurTr,...
                        trialsRun.speed_MMsecSpike_LasttoCurTr,neuronNo,trialNo,minSpeed);
        trialsCueSpikes.Time_LasttoCurTr = ...
            perNeuPerTrRun(trialsCue.res_LasttoCurTr,trialsCue.totclu_LasttoCurTr,...
                        trialsRun.speed_MMsecSpike_LasttoCurTr,neuronNo,trialNo,minSpeed);
        %%%%%% 
        
        if(~isempty(trialNoCueOff))
            trialsCueOffSpikes.Time = ...
                perNeuPerTrRun(trialsCueOff.res,trialsCueOff.totclu,trialsCueOff.speed_MMsecSpike,...
                               neuronNo,trialNoCueOff,minSpeed);
            %%%%%% added by Yingxue on 3/11/2020
            trialsCueOffSpikes.Time_LasttoCurTr = ...
                perNeuPerTrRun(trialsCueOff.res_LasttoCurTr,trialsCueOff.totclu_LasttoCurTr,...
                               trialsRun.speed_MMsecSpike_LasttoCurTr,...
                               neuronNo,trialNoCueOff,minSpeed);
            %%%%%%
        end
                       
        trialsRunSpikes.thPhaseInterpSpike = perNeuPerTrRun(trialsRun.thPhaseInterpSpike,...
                           trialsRun.totclu,trialsRun.speed_MMsecSpike,...
                           neuronNo,trialNo,minSpeed);
        
        % added on 4/27/2020
        trialsRunSpikes.thPhaseHilbSpike = perNeuPerTrRun(trialsRun.thPhaseHilbSpike,...
                           trialsRun.totclu,trialsRun.speed_MMsecSpike,...
                           neuronNo,trialNo,minSpeed);
                       
        trialsRunSpikes.TimeBef = ...
            perNeuPerTrRun(trialsRun.resBef,trialsRun.totcluBef,trialsRun.speed_MMsecSpikeBef,...
                           neuronNo,trialNo,minSpeed);
        trialsRunSpikes.thPhaseInterpSpikeBef = perNeuPerTrRun(trialsRun.thPhaseInterpSpikeBef,...
                           trialsRun.totcluBef,trialsRun.speed_MMsecSpikeBef,...
                           neuronNo,trialNo,minSpeed);
        
        if(~isempty(trialNoCueOff))
            trialsCueOffSpikes.thPhaseInterpSpike = perNeuPerTrRun(trialsCueOff.thPhaseInterpSpike,...
                               trialsCueOff.totclu,trialsCueOff.speed_MMsecSpike,...
                               neuronNo,trialNoCueOff,minSpeed);

            trialsCueOffSpikes.TimeBef = ...
                perNeuPerTrRun(trialsCueOff.resBef,trialsCueOff.totcluBef,trialsCueOff.speed_MMsecSpikeBef,...
                               neuronNo,trialNoCueOff,minSpeed);
            trialsCueOffSpikes.thPhaseInterpSpikeBef = perNeuPerTrRun(trialsCueOff.thPhaseInterpSpikeBef,...
                               trialsCueOff.totcluBef,trialsCueOff.speed_MMsecSpikeBef,...
                               neuronNo,trialNoCueOff,minSpeed);
            [trialsCueOffSpikes.MaxLenT,trialsCueOffSpikes.MeanLenT]...
                = maxTime(trialsCueOff.speed_MMsec,trialNoCueOff,minSpeed);
        
            trialsCueOffSpikes.Dist = ...
            perNeuPerTrRun(trialsCueOff.xMMSpike,trialsCueOff.totclu,trialsCueOff.speed_MMsecSpike,...
                           neuronNo,trialNoCueOff,minSpeed);
                       
            [trialsCueOffSpikes.MaxLenDist,trialsCueOffSpikes.MeanLenDist] = ...
            maxDistanceRun(trialsCueOff.xMM,trialsCueOff.speed_MMsec,trialNoCueOff,minSpeed);
        end

        [trialsRunSpikes.MaxLenT,trialsRunSpikes.MeanLenT]...
            = maxTime(trialsRun.speed_MMsec,trialNo,minSpeed);
        [trialsRewSpikes.MaxLenT,trialsRewSpikes.MeanLenT]...
            = maxTime(trialsRew.speed_MMsec,trialNo,minSpeed);
        [trialsCueSpikes.MaxLenT,trialsCueSpikes.MeanLenT]...
            = maxTime(trialsCue.speed_MMsec,trialNo,minSpeed);
                       
        trialsRunSpikes.Dist = ...
            perNeuPerTrRun(trialsRun.xMMSpike,trialsRun.totclu,trialsRun.speed_MMsecSpike,...
                           neuronNo,trialNo,minSpeed);
        trialsRewSpikes.Dist = ...
            perNeuPerTrRun(trialsRew.xMMSpike,trialsRew.totclu,trialsRew.speed_MMsecSpike,...
                           neuronNo,trialNo,minSpeed);
        trialsCueSpikes.Dist = ...
            perNeuPerTrRun(trialsCue.xMMSpike,trialsCue.totclu,trialsCue.speed_MMsecSpike,...
                           neuronNo,trialNo,minSpeed);
        
        [trialsRunSpikes.MaxLenDist,trialsRunSpikes.MeanLenDist] = ...
            maxDistanceRun(trialsRun.xMM,trialsRun.speed_MMsec,trialNo,minSpeed);
        [trialsRewSpikes.MaxLenDist,trialsRewSpikes.MeanLenDist] = ...
            maxDistanceRun(trialsRew.xMM,trialsRew.speed_MMsec,trialNo,minSpeed);
        [trialsCueSpikes.MaxLenDist,trialsCueSpikes.MeanLenDist] = ...
            maxDistanceRun(trialsCue.xMM,trialsCue.speed_MMsec,trialNo,minSpeed);
        
    end
    
    nNeurons = length(neuronNo);
    nTrials = length(trialNo);
    spikeTrainRun = convertSpikeTtoTrain(trialsRunSpikes.Time,behPar.numSamplesRun,nNeurons,nTrials);
    spikeTrainRew = convertSpikeTtoTrain(trialsRewSpikes.Time,behPar.numSamplesRew,nNeurons,nTrials);
    spikeTrainCue = convertSpikeTtoTrain(trialsCueSpikes.Time,behPar.numSamplesCue,nNeurons,nTrials);
    spikeTrainCueOff = [];
    if(~isempty(trialNoCueOff))
        spikeTrainCueOff = convertSpikeTtoTrain(trialsCueOffSpikes.Time,behPar.numSamplesCueOff,nNeurons,nTrials);
    end
    
    save([path fileName '_alignedSpikesPerNPerT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_Stim' num2str(isStim) '.mat'],...
        'spikeTrainRun','spikeTrainRew','spikeTrainCue','spikeTrainCueOff',...
        'trialsRunSpikes','trialsRewSpikes','trialsCueSpikes','trialsCueOffSpikes',...
        'trialNo','trialNoCueOff','-v7.3');
end

function spikes = perNeuPerTr(res,clu,neuronNo,trialNo)
    spikes = cell(length(neuronNo),max(trialNo));
    for i = neuronNo
        for j = trialNo
            if(~isempty(clu{j}))
                indClu = clu{j} == i; 
                spikes{i,j} = res{j}(indClu);
            end
        end
    end
end

function spikes = perNeuPerTrRun(res,clu,speed,neuronNo,trialNo,minSpeed)
    spikes = cell(length(neuronNo),max(trialNo));
    for i = neuronNo
        for j = trialNo
            if(~isempty(clu{j}))
                if(length(clu{j})~=length(speed{j}))
                    a=0;
                end
%                 disp(['neuron' num2str(i) ' trialNo' num2str(j)])
                indClu = (clu{j} == i) & (speed{j} >= minSpeed); 
                spikes{i,j} = res{j}(indClu);
            end
        end
    end
end

function [maxDist meanDist] = maxDistance(xMM,trialNo)
    maxDistPerTr = zeros(1,max(trialNo));
    for i = trialNo
        if(~isempty(xMM{i}))
            maxDistPerTr(i) = max(xMM{i});
        end
    end
    maxDist = max(maxDistPerTr);
    meanDist = mean(maxDistPerTr);
end

function [maxDist meanDist] = maxDistanceRun(xMM,speed,trialNo,minSpeed)
    maxDistPerTr = zeros(1,max(trialNo));
    for i = trialNo
        if(~isempty(xMM{i}))
            indSpeed = speed{i} >= minSpeed;
            maxDistPerTr(i) = max(xMM{i}(indSpeed));
        end
    end
    maxDist = max(maxDistPerTr);
    meanDist = mean(maxDistPerTr);
end

function [maxT meanT] = maxTime(speed,trialNo,minSpeed)
    maxTAll = zeros(1,max(trialNo));
    for i = trialNo
        if(~isempty(speed{i}))
            maxTAll(i) = find(speed{i} >= minSpeed,1,'last');
        end
    end
    maxT = max(maxTAll);
    meanT = mean(maxTAll);
end

function spikeTrain = convertSpikeTtoTrain(spikes,numSamples,neuronNo,nTrials)
    spikeTrain = cell(1,nTrials);
    for i = 1:nTrials
        if(numSamples(i) == 0)
            continue;
        end
        spikeTrainTmp = zeros(neuronNo,numSamples(i));
        
        for j = 1:neuronNo
            spikeTrainTmp(j,spikes{j,i}) = 1;
        end
        spikeTrain{i} = sparse(spikeTrainTmp);
    end
end
