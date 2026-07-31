function ConvSpikeTrain_Aligned1(path, fileName, onlyRun, mazeSess)
% convolve spike train with gaussian filter in time
% used conv2 instead of conv. Slightly faster, but did not dramatically
% improve the performance
% e.g. ConvSpikeTrain_Aligned1('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',0)

    fullPathSpike = [path fileName '_alignedSpikesPerNPerT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPathSpike) == 0)
        disp('The _alignedSpikesPerNPerT file does not exist');
        return;
    end
    load(fullPathSpike,'trialsRunSpikes');
    
    fileNameConv = [fileName '_convSpikesAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    
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
        
    trialNo = size(trialsRunSpikes.Time,2);
    neuronNo = size(trialsRunSpikes.Time,1);
    nMaxSample = paramC.trialLenT*sampleFq;
    nSamplePerBin = paramC.timeBin*sampleFq;
    timeStep = floor(nSamplePerBin/2):nSamplePerBin:nMaxSample;
    nBins = length(timeStep);
      
    paramC.negTrialLenT = -10; % sec, after including all the spikes from last reward to the end of current trial, and adjusting spike timing
    nMinSample = paramC.negTrialLenT*sampleFq; 
    timeStep1 = nMinSample:nSamplePerBin:nMaxSample;
    nBins1 = length(timeStep1);
    
    tic
    
    for i  = 1:neuronNo
        disp(['Neuron ' num2str(i)]);   
        %% filter spikes aligned to start run
        spikeArray = zeros(trialNo,nBins+2*lenGaussKernel);
        for j = 1:trialNo
            indTime = trialsRunSpikes.Time{i,j} <= nMaxSample;
            if(isempty(indTime))
                continue;
            end
            spikeTime = trialsRunSpikes.Time{i,j}(indTime);  
            spikeTrain = hist(spikeTime,timeStep);
            spikeArray(j,:) = [spikeTrain(nBins-lenGaussKernel+1:nBins)...
                        spikeTrain spikeTrain(1:lenGaussKernel)];   
        end
        filteredSpikeTmp = conv2(spikeArray,paramC.gaussFilt);
        filteredSpikeTmp = ...
                filteredSpikeTmp(:,floor(lenGaussKernel/2)+lenGaussKernel+1:...
                    (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1));     
        filteredSpikeArrayRun{i} = filteredSpikeTmp./paramC.timeBin;
        
        %% filter spikes aligned to start run, including all the spikes from last reward
        % to the end of current trial. Added by Yingxue on 3/13/2020
        spikeArray = zeros(trialNo,nBins1+2*lenGaussKernel);
        for j = 1:trialNo
            indTime = trialsRunSpikes.Time_LasttoCurTr{i,j} <= nMaxSample & ...
                trialsRunSpikes.Time_LasttoCurTr{i,j} >= nMinSample;
            if(isempty(indTime))
                continue;
            end
            spikeTime = trialsRunSpikes.Time_LasttoCurTr{i,j}(indTime);  
            spikeTrain = hist(spikeTime,timeStep1);
            spikeArray(j,:) = [spikeTrain(nBins-lenGaussKernel+1:nBins)...
                        spikeTrain spikeTrain(1:lenGaussKernel)]; 
        end
        filteredSpikeTmp = conv2(spikeArray,paramC.gaussFilt);
        filteredSpikeTmp = ...
                filteredSpikeTmp(:,floor(lenGaussKernel/2)+lenGaussKernel+1:...
                    (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1)); 
        filteredSpikeArrayRun_LasttoCurTr{i} = filteredSpikeTmp./paramC.timeBin;
    end
    
    fullPath = [path fileNameConv];
    save(fullPath, 'filteredSpikeArrayRun','filteredSpikeArrayRun_LasttoCurTr',...
        'paramC','timeStep','timeStep1','-v7.3'); 
    clearvars trialsRunSpikes filteredSpikeArrayRun ...
        filteredSpikeArrayRun_LasttoCurTr 
    
    load(fullPathSpike,'trialsRewSpikes');
    for i = 1:neuronNo
        %% filter spikes aligned to reward
        spikeArray = zeros(trialNo,nBins+2*lenGaussKernel);
        for j = 1:trialNo 
            indTime = trialsRewSpikes.Time{i,j} <= nMaxSample;
            if(isempty(indTime))
                continue;
            end
            spikeTime = trialsRewSpikes.Time{i,j}(indTime);  
            spikeTrain = hist(spikeTime,timeStep);
            spikeArray(j,:) = [spikeTrain(nBins-lenGaussKernel+1:nBins)...
                        spikeTrain spikeTrain(1:lenGaussKernel)];
        end  
        filteredSpikeTmp = conv2(spikeArray,paramC.gaussFilt);
        filteredSpikeTmp = ...
            filteredSpikeTmp(:,floor(lenGaussKernel/2)+lenGaussKernel+1:...
                (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1)); 
        filteredSpikeArrayRew{i} = filteredSpikeTmp./paramC.timeBin;
        
        %% filter spikes aligned to reward, including all the spikes from last reward
        % to the end of current trial. Added by Yingxue on 3/13/2020
        spikeArray = zeros(trialNo,nBins1+2*lenGaussKernel);
        for j = 1:trialNo
            indTime = trialsRewSpikes.Time_LasttoCurTr{i,j} <= nMaxSample & ...
                trialsRewSpikes.Time_LasttoCurTr{i,j} >= nMinSample;
            if(isempty(indTime))
                continue;
            end
            spikeTime = trialsRewSpikes.Time_LasttoCurTr{i,j}(indTime);  
            spikeTrain = hist(spikeTime,timeStep1);
            spikeArray(j,:) = [spikeTrain(nBins-lenGaussKernel+1:nBins)...
                        spikeTrain spikeTrain(1:lenGaussKernel)];
        end
        filteredSpikeTmp = conv2(spikeArray,paramC.gaussFilt);
        filteredSpikeTmp = ...
                filteredSpikeTmp(:,floor(lenGaussKernel/2)+lenGaussKernel+1:...
                    (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1)); 
        filteredSpikeArrayRew_LasttoCurTr{i} = filteredSpikeTmp./paramC.timeBin;
    end
    
    save(fullPath, 'filteredSpikeArrayRew','filteredSpikeArrayRew_LasttoCurTr','-append');
    clearvars trialsRewSpikes filteredSpikeArrayRew ...
        filteredSpikeArrayRew_LasttoCurTr
    
    load(fullPathSpike,'trialsCueSpikes');
    for i = 1:neuronNo
        %% filter spikes aligned to cue
        spikeArray = zeros(trialNo,nBins+2*lenGaussKernel);
        for j = 1:trialNo 
            indTime = trialsCueSpikes.Time{i,j} <= nMaxSample;
            if(isempty(indTime))
                continue;
            end
            spikeTime = trialsCueSpikes.Time{i,j}(indTime);  
            spikeTrain = hist(spikeTime,timeStep);
            spikeArray(j,:) = [spikeTrain(nBins-lenGaussKernel+1:nBins)...
                        spikeTrain spikeTrain(1:lenGaussKernel)];
        end
        filteredSpikeTmp = conv2(spikeArray,paramC.gaussFilt);
        filteredSpikeTmp = ...
                filteredSpikeTmp(:,floor(lenGaussKernel/2)+lenGaussKernel+1:...
                    (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1)); 
        filteredSpikeArrayCue{i} = filteredSpikeTmp./paramC.timeBin;
        
        %% filter spikes aligned to cue, including all the spikes from last reward
        % to the end of current trial. Added by Yingxue on 3/13/2020
        spikeArray = zeros(trialNo,nBins1+2*lenGaussKernel);
        for j = 1:trialNo
            indTime = trialsCueSpikes.Time_LasttoCurTr{i,j} <= nMaxSample & ...
                trialsCueSpikes.Time_LasttoCurTr{i,j} >= nMinSample;
            if(isempty(indTime))
                continue;
            end
            spikeTime = trialsCueSpikes.Time_LasttoCurTr{i,j}(indTime);  
            spikeTrain = hist(spikeTime,timeStep1);
            spikeArray(j,:) = [spikeTrain(nBins-lenGaussKernel+1:nBins)...
                        spikeTrain spikeTrain(1:lenGaussKernel)];
        end
        filteredSpikeTmp = conv2(spikeArray,paramC.gaussFilt);
        filteredSpikeTmp = ...
            filteredSpikeTmp(:,floor(lenGaussKernel/2)+lenGaussKernel+1:...
                (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1)); 
        filteredSpikeArrayCue_LasttoCurTr{i} = filteredSpikeTmp./paramC.timeBin;
        
    end
    
    save(fullPath, 'filteredSpikeArrayCue','filteredSpikeArrayCue_LasttoCurTr','-append');
    clearvars trialsCueSpikes filteredSpikeArrayCue ...
        filteredSpikeArrayCue_LasttoCurTr
    
    toc    
    
end