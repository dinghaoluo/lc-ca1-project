function ConvSpikeTrain_AlignedRunBout(file_path, file_name)
    behave_lfp_path = [file_path '\' file_name(1:17) '_BehavElectrDataLFP.mat'];
    load(behave_lfp_path, 'Spike');
    
    run_bout_path = ['Z:\Dinghao\code\lc-ca1-project\data\run_bouts\'...
                     file_name(1:17) '_run_bouts.mat'];
    load(run_bout_path, 'run_bout_table');
    
    GlobalConst;
    
    paramC.trialLenT = 20; %sec
    paramC.timeBin = 0.0025; %sec
    std = 0.05/paramC.timeBin;   % matching python scripts
    paramC.gaussFilt = gaussFilter(6*std, std);
    lenGaussKernel = length(paramC.gaussFilt);
    normFactor = sum(paramC.gaussFilt);
    paramC.gaussFilt = paramC.gaussFilt./normFactor;
    
    trialNo = size(run_bout_table,1);
    neuronNo = length(unique(Spike.totclu));
    
    nMaxSample = paramC.trialLenT*sampleFq;
    nSamplePerBin = paramC.timeBin*sampleFq;
    timeStepRun = -nSampBef:nSamplePerBin:nMaxSample;
    nBins = length(timeStepRun);
    
    filteredSpikeArrayRunBoutOnSet = {};
    
    for i = 1:neuronNo
        filteredSpikeArrayTmp = zeros(trialNo,nBins);
        spike_neu_i = Spike.res(Spike.totclu == i);
        for j = 1:trialNo
            trial_j_lfp = run_bout_table.run_start_lfp(j);
            trial_j_lfp_range = trial_j_lfp - nSampBef : trial_j_lfp + nMaxSample;
            spikeTime = spike_neu_i(ismember(spike_neu_i, trial_j_lfp_range)) - trial_j_lfp;
            spikeTrain = hist(spikeTime,timeStepRun);
            spikeArray = [spikeTrain(nBins-lenGaussKernel+1:nBins)...
                        spikeTrain spikeTrain(1:lenGaussKernel)];
                
            filteredSpikeTmp = conv(spikeArray,paramC.gaussFilt);
            filteredSpikeArrayTmp(j,:) = ...
                filteredSpikeTmp(floor(lenGaussKernel/2)+lenGaussKernel+1:...
                    (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1)); 
        end
        filteredSpikeArrayRunBoutOnSet{i} = filteredSpikeArrayTmp./paramC.timeBin;
    end
    
    save_path = 'Z:\Dinghao\code\lc-ca1-project\data\run_bouts\fsa_run_bouts\';
    save_path = [save_path file_name(1:17) '_BefRunBout0.mat'];
    save(save_path, 'filteredSpikeArrayRunBoutOnSet','paramC','timeStepRun','-v7.3'); 
   
end
