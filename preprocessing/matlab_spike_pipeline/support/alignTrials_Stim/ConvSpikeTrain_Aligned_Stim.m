function ConvSpikeTrain_Aligned_Stim(path, fileName, onlyRun, mazeSess, isStim)
% convolve spike train with gaussian filter in time
% e.g. ConvSpikeTrain_Aligned('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',0)

    %% added by Yingxue on 2/3/2021
    if(nargin == 4)
        isStim = 0;
    end
    %%
    
    fullPath = [path fileName '_alignedSpikesPerNPerT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_Stim' num2str(isStim) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _alignedSpikesPerNPerT file does not exist');
        return;
    end
    load(fullPath,'trialsRunSpikes','trialsRewSpikes','trialsCueSpikes','trialsCueOffSpikes');
    
    fileNameConv = [fileName '_convSpikesAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_Stim' num2str(isStim) '.mat'];
    
    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    GlobalConst;
    
    paramC.trialLenT = 20; %sec
    paramC.timeBin = 0.0025; %sec
    std = 0.03/paramC.timeBin;
    paramC.gaussFilt = gaussFilter(12*std, std);
    lenGaussKernel = length(paramC.gaussFilt);
    normFactor = sum(paramC.gaussFilt);
    paramC.gaussFilt = paramC.gaussFilt./normFactor;
        
    %% changed by Yingxue on 2/3/2021
    if(isStim == 1)
        trialNo = behPar.indLapsStimInSess;
    else
        trialNo = setdiff(behPar.indLapsSess,behPar.indLapsStimSess);
    end
    %%
    neuronNo = size(trialsRunSpikes.Time,1);
    nMaxSample = paramC.trialLenT*sampleFq;
    nSamplePerBin = paramC.timeBin*sampleFq;
    timeStep = floor(nSamplePerBin/2):nSamplePerBin:nMaxSample;
    nBins = length(timeStep);
      
    paramC.negTrialLenT = -10; % sec, after including all the spikes from last reward to the end of current trial, and adjusting spike timing
    nMinSample = paramC.negTrialLenT*sampleFq; 
    timeStep1 = nMinSample:nSamplePerBin:nMaxSample;
    nBins1 = length(timeStep1);
    
    for i  = 1:neuronNo
        disp(['Neuron ' num2str(i)]);
        %% filter spikes aligned to start run
        filteredSpikeArrayTmp = zeros(trialNo,nBins);
        for j = 1:lengnth(trialNo)
            indTime = trialsRunSpikes.Time{i,trialNo(j)} <= nMaxSample;
            if(isempty(indTime))
                continue;
            end
            spikeTime = trialsRunSpikes.Time{i,trialNo(j)}(indTime);  
            spikeTrain = hist(spikeTime,timeStep);
            spikeArray = [spikeTrain(nBins-lenGaussKernel+1:nBins)...
                        spikeTrain spikeTrain(1:lenGaussKernel)];
                
            filteredSpikeTmp = conv(spikeArray,paramC.gaussFilt);
            filteredSpikeArrayTmp(j,:) = ...
                filteredSpikeTmp(floor(lenGaussKernel/2)+lenGaussKernel+1:...
                    (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1)); 
        end
        filteredSpikeArrayRun{i} = filteredSpikeArrayTmp./paramC.timeBin;
        
        %% filter spikes aligned to start run, including all the spikes from last reward
        % to the end of current trial. Added by Yingxue on 3/13/2020
        filteredSpikeArrayTmp = zeros(trialNo,nBins1);
        for j = 1:length(trialNo)
            indTime = trialsRunSpikes.Time_LasttoCurTr{i,trialNo(j)} <= nMaxSample & ...
                trialsRunSpikes.Time_LasttoCurTr{i,trialNo(j)} >= nMinSample;
            if(isempty(indTime))
                continue;
            end
            spikeTime = trialsRunSpikes.Time_LasttoCurTr{i,trialNo(j)}(indTime);  
            spikeTrain = hist(spikeTime,timeStep1);
            spikeArray = [spikeTrain(nBins-lenGaussKernel+1:nBins)...
                        spikeTrain spikeTrain(1:lenGaussKernel)];
                
            filteredSpikeTmp = conv(spikeArray,paramC.gaussFilt);
            filteredSpikeArrayTmp(j,:) = ...
                filteredSpikeTmp(floor(lenGaussKernel/2)+lenGaussKernel+1:...
                    (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1)); 
        end
        filteredSpikeArrayRun_LasttoCurTr{i} = filteredSpikeArrayTmp./paramC.timeBin;
        
        %% filter spikes aligned to reward
        filteredSpikeArrayTmp = zeros(trialNo,nBins);
        for j = 1:length(trialNo) 
            indTime = trialsRewSpikes.Time{i,trialNo(j)} <= nMaxSample;
            if(isempty(indTime))
                continue;
            end
            spikeTime = trialsRewSpikes.Time{i,trialNo(j)}(indTime);  
            spikeTrain = hist(spikeTime,timeStep);
            spikeArray = [spikeTrain(nBins-lenGaussKernel+1:nBins)...
                        spikeTrain spikeTrain(1:lenGaussKernel)];
                
            filteredSpikeTmp = conv(spikeArray,paramC.gaussFilt);
            filteredSpikeArrayTmp(j,:) = ...
                filteredSpikeTmp(floor(lenGaussKernel/2)+lenGaussKernel+1:...
                    (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1)); 
        end
        filteredSpikeArrayRew{i} = filteredSpikeArrayTmp./paramC.timeBin;
        
        %% filter spikes aligned to reward, including all the spikes from last reward
        % to the end of current trial. Added by Yingxue on 3/13/2020
        filteredSpikeArrayTmp = zeros(trialNo,nBins1);
        for j = 1:length(trialNo)
            indTime = trialsRewSpikes.Time_LasttoCurTr{i,trialNo(j)} <= nMaxSample & ...
                trialsRewSpikes.Time_LasttoCurTr{i,trialNo(j)} >= nMinSample;
            if(isempty(indTime))
                continue;
            end
            spikeTime = trialsRewSpikes.Time_LasttoCurTr{i,trialNo(j)}(indTime);  
            spikeTrain = hist(spikeTime,timeStep1);
            spikeArray = [spikeTrain(nBins-lenGaussKernel+1:nBins)...
                        spikeTrain spikeTrain(1:lenGaussKernel)];
                
            filteredSpikeTmp = conv(spikeArray,paramC.gaussFilt);
            filteredSpikeArrayTmp(j,:) = ...
                filteredSpikeTmp(floor(lenGaussKernel/2)+lenGaussKernel+1:...
                    (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1)); 
        end
        filteredSpikeArrayRew_LasttoCurTr{i} = filteredSpikeArrayTmp./paramC.timeBin;
        
        %% filter spikes aligned to cue
        filteredSpikeArrayTmp = zeros(trialNo,nBins);
        for j = 1:length(trialNo)
            indTime = trialsCueSpikes.Time{i,trialNo(j)} <= nMaxSample;
            if(isempty(indTime))
                continue;
            end
            spikeTime = trialsCueSpikes.Time{i,trialNo(j)}(indTime);  
            spikeTrain = hist(spikeTime,timeStep);
            spikeArray = [spikeTrain(nBins-lenGaussKernel+1:nBins)...
                        spikeTrain spikeTrain(1:lenGaussKernel)];
                
            filteredSpikeTmp = conv(spikeArray,paramC.gaussFilt);
            filteredSpikeArrayTmp(j,:) = ...
                filteredSpikeTmp(floor(lenGaussKernel/2)+lenGaussKernel+1:...
                    (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1)); 
        end
        filteredSpikeArrayCue{i} = filteredSpikeArrayTmp./paramC.timeBin;
        
        %% filter spikes aligned to cue, including all the spikes from last reward
        % to the end of current trial. Added by Yingxue on 3/13/2020
        filteredSpikeArrayTmp = zeros(trialNo,nBins1);
        for j = 1:length(trialNo)
            indTime = trialsCueSpikes.Time_LasttoCurTr{i,trialNo(j)} <= nMaxSample & ...
                trialsCueSpikes.Time_LasttoCurTr{i,trialNo(j)} >= nMinSample;
            if(isempty(indTime))
                continue;
            end
            spikeTime = trialsCueSpikes.Time_LasttoCurTr{i,trialNo(j)}(indTime);  
            spikeTrain = hist(spikeTime,timeStep1);
            spikeArray = [spikeTrain(nBins-lenGaussKernel+1:nBins)...
                        spikeTrain spikeTrain(1:lenGaussKernel)];
                
            filteredSpikeTmp = conv(spikeArray,paramC.gaussFilt);
            filteredSpikeArrayTmp(j,:) = ...
                filteredSpikeTmp(floor(lenGaussKernel/2)+lenGaussKernel+1:...
                    (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1)); 
        end
        filteredSpikeArrayCue_LasttoCurTr{i} = filteredSpikeArrayTmp./paramC.timeBin;
        
    end
    
    fullPath = [path fileNameConv];
    save(fullPath, 'filteredSpikeArrayRun','filteredSpikeArrayRew',...
        'filteredSpikeArrayCue','filteredSpikeArrayRun_LasttoCurTr',...
        'filteredSpikeArrayRew_LasttoCurTr',...
        'filteredSpikeArrayCue_LasttoCurTr','paramC','timeStep','timeStep1','-v7.3'); 
end