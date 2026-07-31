function ConvSpikeTrain_AlignedRunOnset(path, fileName, onlyRun, mazeSess)
% convolve spike train with gaussian filter in time
% e.g. ConvSpikeTrain_Aligned('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',0)

    fullPath = [path fileName '_alignedSpikesPerNPerT_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    if(exist(fullPath) == 0)
        disp('The _alignedSpikesPerNPerT file does not exist');
        return;
    end
    load(fullPath,'trialsRunSpikes');
    
    fileNameConv = [fileName '_convSpikesAligned_msess' num2str(mazeSess) '_BefRun' num2str(onlyRun) '.mat'];
    
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
    timeStepRun = -nSampBef:nSamplePerBin:nMaxSample;
    nBins = length(timeStepRun);
    for i  = 1:neuronNo
        disp(['Neuron ' num2str(i)]);
        %% filter spikes aligned to start run
        filteredSpikeArrayTmp = zeros(trialNo,nBins);
        for j = 1:trialNo
            indTime = trialsRunSpikes.Time{i,j} <= nMaxSample;
            if(isempty(indTime) & isempty(trialsRunSpikes.TimeBef{i,j}))
                continue;
            end
            spikeTime = [trialsRunSpikes.TimeBef{i,j}' ...
                            trialsRunSpikes.Time{i,j}(indTime)'];  
            spikeTrain = hist(spikeTime,timeStepRun);
            spikeArray = [spikeTrain(nBins-lenGaussKernel+1:nBins)...
                        spikeTrain spikeTrain(1:lenGaussKernel)];
                
            filteredSpikeTmp = conv(spikeArray,paramC.gaussFilt);
            filteredSpikeArrayTmp(j,:) = ...
                filteredSpikeTmp(floor(lenGaussKernel/2)+lenGaussKernel+1:...
                    (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1)); 
        end
        filteredSpikeArrayRunOnSet{i} = filteredSpikeArrayTmp./paramC.timeBin;

    end
    
    fullPath = [path fileNameConv];
    save(fullPath, 'filteredSpikeArrayRunOnSet','paramC','timeStepRun','-v7.3'); 
end 