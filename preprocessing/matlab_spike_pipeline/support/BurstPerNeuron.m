function burstPerNeuron = BurstPerNeuronVR(burstPerNeuPerTr,indLapList,numNeurons)
% Extract all the burst spikes of one neuron over all trials
% 
% Inputs:
% burstPerNeuPerTr:        the struct records the information of all
%                          the burst spikes per neuron per trial (referring to function BurstPerNeuronPerTrial)
% indLapList:              indices of all valid laps
%
% Outputs:
% burstPerNeuron:         the struct records the burst spikes for
%                               each neuron

numTrials = length(indLapList);
if((numTrials == 0) || (numNeurons == 0) || isempty(burstPerNeuPerTr))
    burstPerNeuron = [];
    return;
end

burstPerNeuron = struct('indLapList', indLapList,...
                      'numBurstyTrials', zeros(1,numNeurons),... % number of trials where bursts occur
                      'numBurstSp', zeros(1,numNeurons),... % total number of burst spikes
                      'numSp', zeros(1,numNeurons),... % total number of spikes
                      'burstyTrials', zeros(numTrials,numNeurons),... % bursty trials for each neuron
                      'startTimeBurst', [],... % the start time of the each burst
                      'timeBurst', [],... % time of all the spikes within each burst
                      'numSpPerBurst', [],... % number of spikes per burst
                      'startPhaseBurst', [],... % start phase of each burst
                      'phaseBurst', [],... % phase of all the spikes within each burst
                      'isiBurst',[],... % ISI of burst spikes
                      'fractBurstPerTrial',zeros(numTrials,numNeurons),... % fraction of spikes which belongs to a burst    
                      'fractBurstMean',zeros(1,numNeurons),... % mean fraction of spikes which belongs to a burst
                      'fractBurstStd',zeros(1,numNeurons),... % std of the fraction of spikes which belongs to a burst
                      'fractBurst',zeros(1,numNeurons),... % fraction of spikes which belongs to a burst over all the spikes across all the trials
                      ...
                      'StartDistBurst', [],... % the x location of where the first spike of a burst occurs for the maze run
                      'DistBurst', [],... % the x location of where the burst occurs for the maze run
                      ...
                      'meanDire', zeros(1,numNeurons),... % the mean phase direction
                      'meanResultantLen', zeros(1,numNeurons),... % the mean resultant length of the mean phase direction
                      'pRayleigh', ones(1,numNeurons),... % p-value obtained from Rayleigh test
                      'pOmnibus', ones(1,numNeurons),... % p-value obtained from Omnibus test 
                      'meanDireStart', zeros(1,numNeurons),... % the mean phase direction of the first spikes in the bursts
                      'meanResultantLenStart', zeros(1,numNeurons),... % the mean resultant length of the mean phase direction of the first spikes in the bursts
                      'pRayleighStart', ones(1,numNeurons),... % p-value obtained from Rayleigh test of the first spikes in the bursts
                      'pOmnibusStart', ones(1,numNeurons),... % p-value obtained from Omnibus test of the first spikes in the bursts
                      ...
                      'timeNonBurst',[],... % time of non burst spikes
                      'phaseNonBurst',[],... % phase of non burst spikes
                      ...
                      'meanDireNonBurst', zeros(1,numNeurons),... % the mean phase direction for non-burst spikes
                      'meanResultantLenNonBurst', zeros(1,numNeurons),... % the mean resultant length of the mean phase direction for non-burst spikes
                      'pRayleighNonBurst', ones(1,numNeurons),... % p-value obtained from Rayleigh test for non-burst spikes
                      'pOmnibusNonBurst', ones(1,numNeurons)); % p-value obtained from Omnibus test for non-burst spikes 

for i = 1:numNeurons
    burstPerNeuron.startTimeBurst{i} = [];
    burstPerNeuron.timeBurst{i} = [];
    burstPerNeuron.numSpPerBurst{i} = [];
    burstPerNeuron.startPhaseBurst{i} = [];
    burstPerNeuron.phaseBurst{i} = [];
    burstPerNeuron.isiBurst{i} = [];
    burstPerNeuron.timeNonBurst{i} = [];
    burstPerNeuron.phaseNonBurst{i} = [];
    burstPerNeuron.StartDistBurst{i} = [];
    burstPerNeuron.DistBurst{i} = [];
    
    burstPerNeuron.numSp(i) = sum(burstPerNeuPerTr.numSpPerTrial(:,i));
    for j = 1:numTrials
        burstPerNeuron.startTimeBurst{i} = [burstPerNeuron.startTimeBurst{i}...
            burstPerNeuPerTr.startTimeBurst{j,i}];
        burstPerNeuron.timeBurst{i} = [burstPerNeuron.timeBurst{i}...
            burstPerNeuPerTr.timeBurst{j,i}];
        burstPerNeuron.numSpPerBurst{i} = [burstPerNeuron.numSpPerBurst{i}...
            burstPerNeuPerTr.numSpPerBurst{j,i}];
        burstPerNeuron.startPhaseBurst{i} = [burstPerNeuron.startPhaseBurst{i}...
            burstPerNeuPerTr.startPhaseBurst{j,i}];
        burstPerNeuron.phaseBurst{i} = [burstPerNeuron.phaseBurst{i}...
            burstPerNeuPerTr.phaseBurst{j,i}];
        burstPerNeuron.isiBurst{i} = [burstPerNeuron.isiBurst{i}...
            burstPerNeuPerTr.isiBurst{j,i}];
        
        burstPerNeuron.timeNonBurst{i} = [burstPerNeuron.timeNonBurst{i}...
            burstPerNeuPerTr.timeNonBurst{j,i}];
        burstPerNeuron.phaseNonBurst{i} = [burstPerNeuron.phaseNonBurst{i}...
            burstPerNeuPerTr.phaseNonBurst{j,i}];
        
        if(~isempty(burstPerNeuPerTr.startTimeBurst{j,i}))
            burstPerNeuron.numBurstyTrials(i) = burstPerNeuron.numBurstyTrials(i) + 1;
            burstPerNeuron.burstyTrials(j,i) = 1;
        end
        
        burstPerNeuron.StartDistBurst{i} = [burstPerNeuron.StartDistBurst{i}...
            burstPerNeuPerTr.StartDistBurst{j,i}];
        burstPerNeuron.DistBurst{i} = [burstPerNeuron.DistBurst{i}...
            burstPerNeuPerTr.DistBurst{j,i}];
    end
    
    burstPerNeuron.fractBurstPerTrial(:,i) = burstPerNeuPerTr.fractBurst(:,i);
    burstPerNeuron.fractBurstMean(i) = mean(burstPerNeuPerTr.fractBurst(:,i));
    burstPerNeuron.fractBurstStd(i) = std(burstPerNeuPerTr.fractBurst(:,i));
    if(~isempty(burstPerNeuron.timeBurst{i}))
        burstPerNeuron.meanDire(i) = circ_mean(burstPerNeuron.phaseBurst{i}');
        burstPerNeuron.meanResultantLen(i) = circ_r(burstPerNeuron.phaseBurst{i}');
        burstPerNeuron.pRayleigh(i) = circ_rtest(burstPerNeuron.phaseBurst{i}');
        burstPerNeuron.pOmnibus(i) = circ_otest(burstPerNeuron.phaseBurst{i}');
        
        burstPerNeuron.meanDireStart(i) = circ_mean(burstPerNeuron.startPhaseBurst{i}');
        burstPerNeuron.meanResultantLenStart(i) = circ_r(burstPerNeuron.startPhaseBurst{i}');
        burstPerNeuron.pRayleighStart(i) = circ_rtest(burstPerNeuron.startPhaseBurst{i}');
        burstPerNeuron.pOmnibusStart(i) = circ_otest(burstPerNeuron.startPhaseBurst{i}');
        
        burstPerNeuron.numBurstSp(i) = length(burstPerNeuron.timeBurst{i});
        burstPerNeuron.fractBurst(i) = burstPerNeuron.numBurstSp(i)/...
            burstPerNeuron.numSp(i);
    end
    
    if(~isempty(burstPerNeuron.timeNonBurst{i}))
        burstPerNeuron.meanDireNonBurst(i) = circ_mean(burstPerNeuron.phaseNonBurst{i}');
        burstPerNeuron.meanResultantLenNonBurst(i) = circ_r(burstPerNeuron.phaseNonBurst{i}');
        burstPerNeuron.pRayleighNonBurst(i) = circ_rtest(burstPerNeuron.phaseNonBurst{i}');
        burstPerNeuron.pOmnibusNonBurst(i) = circ_otest(burstPerNeuron.phaseNonBurst{i}');
    end
end