function ConvSpikeTrainDist_Aligned2P(path, fileName,spaceBin,onlyRun,mazeSess)
% convolve spike train with gaussian filter in time
% e.g. ConvSpikeTrainDist_Aligned('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',20,1)

    fullPath = [path fileName '_timePerDistBin_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _timePerDistBin_Run file does not exist');
        return;
    end
    load(fullPath,'timePerDistBinRun','timePerDistBinRew','timePerDistBinCue');
    
    fileNameConv = [fileName '_convSpikesDistAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    
    GlobalConst2P;
    
    std = floor(spaceBin/spaceMergeBin/2);
    paramC.gaussFilt = gaussFilter2P(12*std, std);
    lenGaussKernel = length(paramC.gaussFilt);
    normFactor = sum(paramC.gaussFilt);
    paramC.gaussFilt = paramC.gaussFilt./normFactor;
    paramC.trackLen = trackLen; %mm
    
    if(spaceMergeBin ~= 0)
        paramC.spaceSteps = 0:spaceMergeBin:paramC.trackLen;
    else
        paramC.spaceSteps = 0:paramC.trackLen;
    end
    paramC.spaceBin = spaceBin;
    
    trialNo = length(timePerDistBinRun.spikesHist);
    neuronNo = size(timePerDistBinRun.spikesHist{end},1);
    nBins = length(paramC.spaceSteps);
    
    distRun = timePerDistBinRun.spikesHist;
    distRew = timePerDistBinRew.spikesHist;
    distCue = timePerDistBinCue.spikesHist;
    timeRun = timePerDistBinRun.hist;
    timeRew = timePerDistBinRew.hist;
    timeCue = timePerDistBinCue.hist;
    
    for i  = 1:neuronNo
        disp(['Neuron ' num2str(i)]);
        %% filter spikes aligned to start run
        filteredSpikeArrayTmp = zeros(trialNo,nBins);
        filteredSpikeArrayTmp1 = zeros(trialNo,nBins);
        for j = 1:trialNo
            if(isempty(timeRun{j}))
                continue;
            end
            
            timeRunTmp = timeRun{j};
            avgTimePerBin = sum(timeRunTmp)/(paramC.trackLen/spaceMergeBin)/sampleFq; % added on 1/7/2022 by Yingxue 
            
            spikeTrain = distRun{j}(i,:);  
            trLen = size(distRun{j},2);
            spikeArray = [spikeTrain(trLen-lenGaussKernel+1:trLen)...
                        spikeTrain spikeTrain(1:lenGaussKernel)]/avgTimePerBin; % added on 1/7/2022 by Yingxue
                
            filteredSpikeTmp = conv(spikeArray,paramC.gaussFilt);
            filteredSpikeTmp = filteredSpikeTmp(floor(lenGaussKernel/2)+lenGaussKernel+1:...
                    (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1));
            ind = 1:min(trLen,nBins);
            filteredSpikeArrayTmp(j,ind) = filteredSpikeTmp(ind);
                
            % normT            
            spikeTrain = spikeTrain./timeRunTmp*sampleFq; % added *sampleFq on 1/3/2022
            spikeTrain(isnan(spikeTrain)) = 0;
            spikeArray = [spikeTrain(trLen-lenGaussKernel+1:trLen)...
                        spikeTrain spikeTrain(1:lenGaussKernel)];
                
            filteredSpikeTmp = conv(spikeArray,paramC.gaussFilt);
            filteredSpikeTmp = filteredSpikeTmp(floor(lenGaussKernel/2)+lenGaussKernel+1:...
                    (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1));
            ind = 1:min(trLen,nBins);
            filteredSpikeArrayTmp1(j,ind) = filteredSpikeTmp(ind);
        end
        filteredSpikeDistArrayRun{i} = filteredSpikeArrayTmp;
        filteredSpikeDistNormTArrayRun{i} = filteredSpikeArrayTmp1;
        
        %% filter spikes aligned to reward
        filteredSpikeArrayTmp = zeros(trialNo,nBins);
        filteredSpikeArrayTmp1 = zeros(trialNo,nBins);
        for j = 1:trialNo 
            if(isempty(timeRew{j}))
                continue;
            end
            
            timeRewTmp = timeRew{j};
            avgTimePerBin = sum(timeRewTmp)/(paramC.trackLen/spaceMergeBin)/sampleFq; % added on 1/7/2022 by Yingxue 
            
            spikeTrain = distRew{j}(i,:);  
            trLen = size(distRew{j},2);
            spikeArray = [spikeTrain(trLen-lenGaussKernel+1:trLen)...
                        spikeTrain spikeTrain(1:lenGaussKernel)]/avgTimePerBin; % added on 1/7/2022 by Yingxue
                
            filteredSpikeTmp = conv(spikeArray,paramC.gaussFilt);
            filteredSpikeTmp = filteredSpikeTmp(floor(lenGaussKernel/2)+lenGaussKernel+1:...
                    (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1));
            ind = 1:min(trLen,nBins);
            filteredSpikeArrayTmp(j,ind) = filteredSpikeTmp(ind);
               
            % normT            
            spikeTrain = spikeTrain./timeRewTmp*sampleFq; % added *sampleFq on 1/3/2022
            spikeTrain(isnan(spikeTrain)) = 0;
            spikeArray = [spikeTrain(trLen-lenGaussKernel+1:trLen)...
                        spikeTrain spikeTrain(1:lenGaussKernel)];
                
            filteredSpikeTmp = conv(spikeArray,paramC.gaussFilt);
            filteredSpikeTmp = filteredSpikeTmp(floor(lenGaussKernel/2)+lenGaussKernel+1:...
                    (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1));
            ind = 1:min(trLen,nBins);
            filteredSpikeArrayTmp1(j,ind) = filteredSpikeTmp(ind);
        end
        filteredSpikeDistArrayRew{i} = filteredSpikeArrayTmp;
        filteredSpikeDistNormTArrayRew{i} = filteredSpikeArrayTmp1;
        
        
        %% filter spikes aligned to cue
        filteredSpikeArrayTmp = zeros(trialNo,nBins);
        filteredSpikeArrayTmp1 = zeros(trialNo,nBins);
        for j = 1:trialNo 
            if(isempty(timeCue{j}))
                continue;
            end
            
            timeCueTmp = timeCue{j};
            avgTimePerBin = sum(timeCueTmp)/paramC.trackLen/sampleFq; % added on 1/7/2022 by Yingxue 
                        
            spikeTrain = distCue{j}(i,:);  
            trLen = size(distCue{j},2);
            spikeArray = [spikeTrain(trLen-lenGaussKernel+1:trLen)...
                        spikeTrain spikeTrain(1:lenGaussKernel)]/avgTimePerBin; % added on 1/7/2022 by Yingxue
                
            filteredSpikeTmp = conv(spikeArray,paramC.gaussFilt);
            filteredSpikeTmp = filteredSpikeTmp(floor(lenGaussKernel/2)+lenGaussKernel+1:...
                    (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1));
            ind = 1:min(trLen,nBins);
            filteredSpikeArrayTmp(j,ind) = filteredSpikeTmp(ind);
                
            % normT
            spikeTrain = spikeTrain./timeCueTmp*sampleFq; % added *sampleFq on 1/3/2022
            spikeTrain(isnan(spikeTrain)) = 0;
            spikeArray = [spikeTrain(trLen-lenGaussKernel+1:trLen)...
                        spikeTrain spikeTrain(1:lenGaussKernel)];
                
            filteredSpikeTmp = conv(spikeArray,paramC.gaussFilt);
            filteredSpikeTmp = filteredSpikeTmp(floor(lenGaussKernel/2)+lenGaussKernel+1:...
                    (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1));
            ind = 1:min(trLen,nBins);
            filteredSpikeArrayTmp1(j,ind) = filteredSpikeTmp(ind);
        end
        filteredSpikeDistArrayCue{i} = filteredSpikeArrayTmp;
        filteredSpikeDistNormTArrayCue{i} = filteredSpikeArrayTmp1;
    end
    
    fullPath = [path fileNameConv];
    save(fullPath, 'filteredSpikeDistArrayRun','filteredSpikeDistArrayRew',...
        'filteredSpikeDistArrayCue','filteredSpikeDistNormTArrayRun',...
        'filteredSpikeDistNormTArrayRew','filteredSpikeDistNormTArrayCue',...
        'paramC','-v7.3'); 
end
