function ConvSpikeTrain_AlignedReward2P(path, fileName, onlyRun, mazeSess)
% convolve spike train with gaussian filter in time
% e.g. ConvSpikeTrain_Aligned('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',0)

    fileNameRun = [fileName '_alignRew_msess' num2str(mazeSess) '.mat']; 
    fullPath = [path fileNameRun];
    if(exist(fullPath) == 0)
        disp('The alignToReward2P file does not exist');
        return;
    end
    load(fullPath,'trialsRew');
    
    fileNameConv = [fileName '_convSpikesAligned_msess' num2str(mazeSess) '_BefRew' num2str(onlyRun) '.mat'];
    
    GlobalConst2P;
    
    paramC.trialLenT = 20; %sec
    paramC.timeBin = timeStep; %sec
    std = timeBin/paramC.timeBin;
    paramC.gaussFilt = gaussFilter2P(12*std, std);
    lenGaussKernel = length(paramC.gaussFilt);
    normFactor = sum(paramC.gaussFilt);
    paramC.gaussFilt = paramC.gaussFilt./normFactor;
        
    trialNo = length(trialsRew.spikes);
    neuronNo = size(trialsRew.spikes{end},2);
    nMaxSample = paramC.trialLenT*sampleFq;
    
    for i  = 1:neuronNo
        disp(['Neuron ' num2str(i)]);
        %% filter spikes aligned to start run
        filteredSpikeArrayTmp = zeros(trialNo,nMaxSample);
        dFFArrayTmp = zeros(trialNo,nMaxSample);
        for j = 1:trialNo
            if(isempty(trialsRew.spikes{j}))
                continue;
            end
            spikeTime = [trialsRew.spikesBef{j}(:,i)' ...
                            trialsRew.spikes{j}(:,i)'];     
            dFF = [trialsRew.dFFBefGF{j}(:,i)' ...
                            trialsRew.dFFGF{j}(:,i)'];     % changed to dFF gaussian filtered instead of using the thresholded dFF on 9/6/2023
            if(onlyRun == 1)
                speed = [trialsRew.speed_MMsecBef{j}' trialsRew.speed_MMsec{j}'];
                indSpeed = speed > minSpeed;
                spikeTime = spikeTime(indSpeed);
                dFF = dFF(indSpeed);
            end
            numSamples = min(length(spikeTime),nMaxSample);
            spikeTime = spikeTime(1:numSamples);
            dFF = dFF(1:numSamples);
            
            spikeArray = [spikeTime(numSamples-lenGaussKernel+1:numSamples)...
                        spikeTime spikeTime(1:lenGaussKernel)];
                    
            filteredSpikeTmp = conv(spikeArray,paramC.gaussFilt);
            filteredSpikeArrayTmp(j,1:numSamples) = ...
                filteredSpikeTmp(floor(lenGaussKernel/2)+lenGaussKernel+1:...
                    (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1)); 
                                
            dFFArrayTmp(j,1:numSamples) = dFF;
        end
        filteredSpikeArrayReward{i} = filteredSpikeArrayTmp./paramC.timeBin;
        dFFArrayReward{i} = dFFArrayTmp;
    end
    
    fullPath = [path fileNameConv];
    save(fullPath, 'filteredSpikeArrayReward','dFFArrayReward','paramC','-v7.3'); 
end 