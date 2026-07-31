function ConvSpikeTrain_Aligned2P(path, fileName, onlyRun, mazeSess)
% convolve spike train with gaussian filter in time
% e.g. ConvSpikeTrain_Aligned('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',0)

    fileNameRun = [fileName '_alignRun_msess' num2str(mazeSess) '.mat']; 
    fullPath = [path fileNameRun];
    if(exist(fullPath) == 0)
        disp('The alignToRunOnset2P file does not exist');
        return;
    end
    load(fullPath,'trialsRun');
           
    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    GlobalConst2P;
    
    paramC.trialLenT = 20; %sec
    paramC.timeBin = timeStep; %sec
    std = timeBin/paramC.timeBin;
    paramC.gaussFilt = gaussFilter2P(12*std, std);
    lenGaussKernel = length(paramC.gaussFilt);
    normFactor = sum(paramC.gaussFilt);
    paramC.gaussFilt = paramC.gaussFilt./normFactor;
        
    trialNo = length(trialsRun.spikes);
    neuronNo = size(trialsRun.spikes{end},2);
    nMaxSample = round(paramC.trialLenT*sampleFq);
              
    disp('Convolve spike trains with Gaussian kernel (Run)')
    [filteredSpikeArrayRun,dFFArrayRun] = convSpikeTrain(trialsRun,trialNo,neuronNo,nMaxSample,onlyRun,minSpeed,...
                    paramC.gaussFilt,lenGaussKernel,paramC.timeBin);
    
    fileNameConv = [fileName '_convSpikesAlignedRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    fullPath = [path fileNameConv];
    save(fullPath, 'filteredSpikeArrayRun','dFFArrayRun','paramC','-v7.3'); 
    clear filteredSpikeArrayRun dFFArrayRun trialsRun
      
    %%
    disp('Convolve spike trains with Gaussian kernel (Rew)')
    fileNameRew = [fileName '_alignRew_msess' num2str(mazeSess) '.mat']; 
    fullPath = [path fileNameRew];
    if(exist(fullPath) == 0)
        disp('The alignToReward2P file does not exist');
        return;
    end
    load(fullPath,'trialsRew');
    
    [filteredSpikeArrayRew,dFFArrayRew] = convSpikeTrain(trialsRew,trialNo,neuronNo,nMaxSample,onlyRun,minSpeed,...
                    paramC.gaussFilt,lenGaussKernel,paramC.timeBin);       
    
    fileNameConv = [fileName '_convSpikesAlignedRew_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    fullPath = [path fileNameConv];
    save(fullPath, 'filteredSpikeArrayRew','dFFArrayRew','-v7.3'); 
    clear filteredSpikeArrayRew dFFArrayRew trialsRew
    
    %%
    disp('Convolve spike trains with Gaussian kernel (Cue)')
    fileNameCue = [fileName '_alignCue_msess' num2str(mazeSess) '.mat']; 
    fullPath = [path fileNameCue];
    if(exist(fullPath) == 0)
        disp('The alignToCue2P file does not exist');
        return;
    end
    load(fullPath,'trialsCue');
    
    [filteredSpikeArrayCue,dFFArrayCue] = convSpikeTrain(trialsCue,trialNo,neuronNo,nMaxSample,onlyRun,minSpeed,...
                    paramC.gaussFilt,lenGaussKernel,paramC.timeBin);
   
    fileNameConv = [fileName '_convSpikesAlignedCue_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    fullPath = [path fileNameConv];
    save(fullPath,'filteredSpikeArrayCue','dFFArrayCue','-v7.3');
    clear filteredSpikeArrayCue dFFArrayCue trialsCue
end

function [filteredSpikeArrayRun,dFFArrayRun] = convSpikeTrain(trialsRun,trialNo,neuronNo,nMaxSample,onlyRun,minSpeed,...
                    gaussFilt,lenGaussKernel,timeBin)

    for i  = 1:neuronNo
        disp(['Neuron ' num2str(i)]);   
        %% filter spikes aligned to start run
        filteredSpikeArrayTmp = zeros(trialNo,nMaxSample);
        dFFArrayTmp = zeros(trialNo,nMaxSample);
        for j = 1:trialNo
            if(isempty(trialsRun.spikes{j}))
                continue;
            end
            if(onlyRun == 1)
                spike = trialsRun.spikes{j}(trialsRun.speed_MMsec{j} > minSpeed,i)';
            else
                spike = trialsRun.spikes{j}(:,i)';
            end
            numSamples = min(length(spike),nMaxSample);
            spike = spike(1:numSamples);  
            spikeArray = [spike(numSamples-lenGaussKernel+1:numSamples)...
                        spike spike(1:lenGaussKernel)];
                    
            filteredSpikeTmp = conv(spikeArray,gaussFilt);
            filteredSpikeArrayTmp(j,1:numSamples) = ...
                filteredSpikeTmp(floor(lenGaussKernel/2)+lenGaussKernel+1:...
                    (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1));  
                 
            dFFArrayTmp(j,1:numSamples) = trialsRun.dFFGF{j}(1:numSamples,i)';  % changed to dFF gaussian filtered instead of using the thresholded dFF on 9/6/2023
        end
        
        filteredSpikeArrayRun{i} = filteredSpikeArrayTmp/timeBin; % added *sampleFq on 1/7/2022 ;
        dFFArrayRun{i} = dFFArrayTmp;
        clear filteredSpikeArrayTmp dFFArrayTmp spike spikeArray filteredSpikeTmp
    end

end