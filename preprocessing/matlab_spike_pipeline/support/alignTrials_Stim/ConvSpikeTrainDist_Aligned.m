function ConvSpikeTrainDist_Aligned(path, fileName,spaceBin,onlyRun,mazeSess)
% convolve spike train with gaussian filter in time
% e.g. ConvSpikeTrainDist_Aligned('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',20,1)

    fullPath = [path fileName '_alignedSpikesPerNPerT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _alignedSpikesPerNPerT file does not exist');
        return;
    end
    load(fullPath,'trialsRunSpikes','trialsRewSpikes','trialsCueSpikes');
    
    fullPath = [path fileName '_timePerDistBin_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _timePerDistBin_Run file does not exist');
        return;
    end
    load(fullPath,'timePerDistBinRun','timePerDistBinRew','timePerDistBinCue');
    
    fileNameConv = [fileName '_convSpikesDistAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    
%     fullPath = [path fileName '_alignRun_msess' num2str(mazeSess) '.mat']; 
%     if(exist(fullPath) == 0)
%         disp('The aligned to run file does not exist');
%         return;
%     end
%     load(fullPath,'trialsRun');
%     
%     fullPath = [path fileName '_alignRew_msess' num2str(mazeSess) '.mat']; 
%     if(exist(fullPath) == 0)
%         disp('The aligned to reward file does not exist');
%         return;
%     end
%     load(fullPath,'trialsRew');
%     
%     fullPath = [path fileName '_alignCue_msess' num2str(mazeSess) '.mat']; 
%     if(exist(fullPath) == 0)
%         disp('The aligned to cue file does not exist');
%         return;
%     end
%     load(fullPath,'trialsCue');
    
    GlobalConst;
    
    std = floor(spaceBin/spaceMergeBin/2);
    paramC.gaussFilt = gaussFilter(12*std, std);
    lenGaussKernel = length(paramC.gaussFilt);
    normFactor = sum(paramC.gaussFilt);
    paramC.gaussFilt = paramC.gaussFilt./normFactor;
    paramC.trackLen = 2000; %mm
    
    if(spaceMergeBin ~= 0)
        paramC.spaceSteps = 0:spaceMergeBin:paramC.trackLen;
    else
        paramC.spaceSteps = 0:paramC.trackLen;
    end
    paramC.spaceBin = spaceBin;
    
    trialNo = size(trialsRunSpikes.Time,2);
    neuronNo = size(trialsRunSpikes.Time,1);
    nBins = length(paramC.spaceSteps);
    
    distRun = trialsRunSpikes.Dist;
    distRew = trialsRewSpikes.Dist;
    distCue = trialsCueSpikes.Dist;
    timeRun = timePerDistBinRun.hist;
    timeRew = timePerDistBinRew.hist;
    timeCue = timePerDistBinCue.hist;
    for i  = 1:neuronNo
        disp(['Neuron ' num2str(i)]);
        %% filter spikes aligned to start run
        filteredSpikeArrayTmp = zeros(trialNo,nBins);
        filteredSpikeArrayTmp1 = zeros(trialNo,nBins);
        for j = 1:trialNo
            indDist = distRun{i,j} <= paramC.trackLen;
            if(isempty(indDist))
                continue;
            end
            spikeDist = distRun{i,j}(indDist);  
            spikeTrain = hist(spikeDist,paramC.spaceSteps);
            spikeArray = [spikeTrain(nBins-lenGaussKernel+1:nBins)...
                        spikeTrain spikeTrain(1:lenGaussKernel)];
                
            filteredSpikeTmp = conv(spikeArray,paramC.gaussFilt);
            filteredSpikeArrayTmp(j,:) = ...
                filteredSpikeTmp(floor(lenGaussKernel/2)+lenGaussKernel+1:...
                    (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1));
              
            % normT
            timeRunTmp = ones(1,nBins);
            lenTr = length(timeRun{j});
            if(lenTr <= nBins)
                timeRunTmp(1:lenTr) = timeRun{j};
            else
                timeRunTmp = timeRun{j}(1:nBins);
            end
            spikeTrain = spikeTrain./timeRunTmp;
            spikeArray = [spikeTrain(nBins-lenGaussKernel+1:nBins)...
                        spikeTrain spikeTrain(1:lenGaussKernel)];
                
            filteredSpikeTmp = conv(spikeArray,paramC.gaussFilt);
            filteredSpikeArrayTmp1(j,:) = ...
                filteredSpikeTmp(floor(lenGaussKernel/2)+lenGaussKernel+1:...
                    (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1));
        end
        filteredSpikeDistArrayRun{i} = filteredSpikeArrayTmp;
        filteredSpikeDistNormTArrayRun{i} = filteredSpikeArrayTmp1;
        
        %% filter spikes aligned to reward
        filteredSpikeArrayTmp = zeros(trialNo,nBins);
        filteredSpikeArrayTmp1 = zeros(trialNo,nBins);
        for j = 1:trialNo 
            indDist = distRew{i,j} <= paramC.trackLen;
            if(isempty(indDist))
                continue;
            end
            spikeDist = distRew{i,j}(indDist);  
            spikeTrain = hist(spikeDist,paramC.spaceSteps);
            spikeArray = [spikeTrain(nBins-lenGaussKernel+1:nBins)...
                        spikeTrain spikeTrain(1:lenGaussKernel)];
                
            filteredSpikeTmp = conv(spikeArray,paramC.gaussFilt);
            filteredSpikeArrayTmp(j,:) = ...
                filteredSpikeTmp(floor(lenGaussKernel/2)+lenGaussKernel+1:...
                    (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1)); 
               
            % normT
            timeRewTmp = ones(1,nBins);
            lenTr = length(timeRew{j});
            if(lenTr <= nBins)
                timeRewTmp(1:lenTr) = timeRew{j};
            else
                timeRewTmp = timeRew{j}(1:nBins);
            end
            spikeTrain = spikeTrain./timeRewTmp;
            spikeArray = [spikeTrain(nBins-lenGaussKernel+1:nBins)...
                        spikeTrain spikeTrain(1:lenGaussKernel)];
                
            filteredSpikeTmp = conv(spikeArray,paramC.gaussFilt);
            filteredSpikeArrayTmp1(j,:) = ...
                filteredSpikeTmp(floor(lenGaussKernel/2)+lenGaussKernel+1:...
                    (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1));
        end
        filteredSpikeDistArrayRew{i} = filteredSpikeArrayTmp;
        filteredSpikeDistNormTArrayRew{i} = filteredSpikeArrayTmp1;
        
        
        %% filter spikes aligned to cue
        filteredSpikeArrayTmp = zeros(trialNo,nBins);
        filteredSpikeArrayTmp1 = zeros(trialNo,nBins);
        for j = 1:trialNo 
            indDist = distCue{i,j} <= paramC.trackLen;
            if(isempty(indDist))
                continue;
            end
            spikeDist = distCue{i,j}(indDist);  
            spikeTrain = hist(spikeDist,paramC.spaceSteps);
            spikeArray = [spikeTrain(nBins-lenGaussKernel+1:nBins)...
                        spikeTrain spikeTrain(1:lenGaussKernel)];
                
            filteredSpikeTmp = conv(spikeArray,paramC.gaussFilt);
            filteredSpikeArrayTmp(j,:) = ...
                filteredSpikeTmp(floor(lenGaussKernel/2)+lenGaussKernel+1:...
                    (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1)); 
                
            % normT
            timeCueTmp = ones(1,nBins);
            lenTr = length(timeCue{j});
            if(lenTr <= nBins)
                timeCueTmp(1:lenTr) = timeCue{j};
            else
                timeCueTmp = timeCue{j}(1:nBins);
            end
            spikeTrain = spikeTrain./timeCueTmp;
            spikeArray = [spikeTrain(nBins-lenGaussKernel+1:nBins)...
                        spikeTrain spikeTrain(1:lenGaussKernel)];
                
            filteredSpikeTmp = conv(spikeArray,paramC.gaussFilt);
            filteredSpikeArrayTmp1(j,:) = ...
                filteredSpikeTmp(floor(lenGaussKernel/2)+lenGaussKernel+1:...
                    (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1));
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

function timePerBin = spikeTime2Dist(xMM, spaceSteps)
                        
    numBins = length(spaceSteps);
    step = spaceSteps(2) - spaceSteps(1);
    timePerBin = zeros(1,numBins);
    for i = 1:numBins
        ind = find(xMM >= spaceSteps(i)-step/2 & xMM < spaceSteps(i)+step/2);
        if(~isempty(ind))
            timePerBin(i) = length(ind);
        else
            timePerBin(i) = 1;
        end
    end 
end