function spikeThetaPhaseStruct = SpikeThetaPhaseRunOnset...
    (spikesBef, spikes, spikeThetaBef, spikeTheta, thetaCycleStruct, indLapList, ...
     numNeurons)
% calculate the theta phase of each spikes
% 1. vs sample number
% 2. vs theta cycle number
%
% Inputs: 
% spikes:               spike structure which contains a two dimensional cell
%                       structure, with spikes{i,j} representing the spike train
%                       from neuron i in trial j 
% theta:                theta phase structure including theta phase at each
%                       sampling point of each valid trial, theta{i} is the
%                       theta phase array from trial i
% thetaCycleStruct:     structure containing the information about all the
%                       theta cycles within each trials (referring to 
%                       function ThetaCycle)
% indLapList:           the array containing the index numbers of the valid
%                       trials (Caution: since theta strcuture records the 
%                       theta traces from all the valid trials in 
%                       the task, indLapList here includes the indices of
%                       all to be selected cells in the theta structure, 
%                       rather than referring back to the original trial no.)
% numNeurons:           number of neurons
% numSamples:           only consider the samples between 0 and numSamples in the
%                       theta phase calculation (when the numSamples > 
%                       minNumSamples, set the input thetaCycleStruct = [])
% timeStep:             sampling time step
% figureState:          0: figure off
%                       2: figure on
%
% Outputs: spikeThetaPhaseStruct with following field
% indLapList: 
% spTimePerTrialPerNeuron:      spike time per trial per neuron (sample <
%                               numSamples), spTimePerTrialPerNeuron{i,j} is for neuron i trial j
% spPhaseVsTPerTrialPerNeuron:  spike phases per trial per neuron (sample < numSamples)
% spCyclePerTrialPerNeuron:     spike cycle per trial per neuron (cycle < thetaCycleStruct.numCyclesMin)
% spPhaseVsCPerTrialPerNeuron:  spike phases per trial per neuron (cycle < thetaCycleStruct.numCyclesMin)
% spTimePerNeuron:              spike time per neuron (sample < numSamples)
% spPhaseVsTPerNeuron:          spike phases per neuron (sample < numSamples)
% spCyclePerNeuron:             spike cycle per neuron (cycle < thetaCycleStruct.numCyclesMin)
% spPhaseVsCPerNeuron:          spike phases per neuron (cycle < thetaCycleStruct.numCyclesMin)
% minPhaseArr:                  phase with min spike activities for each
%                               neuron (estimated based on histogram)
% minPhaseFilArr:               phase with min spike activities for each
%                               neuron (estimated based on filtered histogram)
% meanDire:                     mean phase direction
% meanResultantLen:             mean resultant length of the mean phase
%                               direction
% pRayleigh:                    p-value obtained from Rayleigh test
% pOmnibus:                     p-value obtained from Omnibus test 


%%%% first check whether spikes, theta and thetaCycleStruct contain the
%%%% same number of trials
sizeTrialsSp = size(spikes,2);
sizeTrialsTheta = size(spikeTheta,2);
spikeThetaPhaseStruct = [];

if(sizeTrialsSp ~= sizeTrialsTheta)
    disp('spikes and theta should contain the same number of trials');
    return;
end
if(~isempty(thetaCycleStruct))
    sizeTrialsThetaCycle = length(thetaCycleStruct.numCyclesPerTrial);
    if(sizeTrialsTheta < sizeTrialsThetaCycle)
        disp('thetaCycleStruct should contain less or equal number of trials than spikeTheta');
        return;
    end
end

GlobalConst;

%%%% check whether trial no. in indLapList is beyond the range of trials
%%%% in the spike structure
maxTrialNo = max(indLapList);
if(maxTrialNo > sizeTrialsSp)
    disp('indLapList contains trials which is not included in the spikes structure');
    return;
end

numTrials = length(indLapList);
if(numTrials == 0)
    return;
end

spikeThetaPhaseStruct = struct('indLapList',indLapList,...
                               ...
                               'thetaCycleAroundZero',zeros(1,numTrials),...
                               'indCycleStartAroundZero',zeros(1,numTrials),...
                               ...
                               'spTimePerNeuron',{cell(1,numNeurons)},...
                               'spPhaseVsTPerNeuron',{cell(1,numNeurons)},...
                               ...
                               'spPhaseVsCPerTrialPerNeuron',{cell(numNeurons,numTrials)},...
                               'spCyclePerTrialPerNeuron',{cell(numNeurons,numTrials)},...
                               'spCycleRunOnPerTrialPerNeuron',{cell(numNeurons,numTrials)},...
                               'spMeanPhaseVsCPerTrialPerNeuron',{cell(numNeurons,numTrials)},...
                               'spMeanCyclePerTrialPerNeuron',{cell(numNeurons,numTrials)},...
                               'spMeanCycleRunOnPerTrialPerNeuron',{cell(numNeurons,numTrials)},...
                               ...
                               'spCyclePerNeuron',{cell(1,numNeurons)},...
                               'spCycleRunOnPerNeuron',{cell(1,numNeurons)},...
                               'spPhaseVsCPerNeuron',{cell(1,numNeurons)},...
                               'spMeanPhaseVsCPerNeuron',{cell(1,numNeurons)},...
                               'spMeanCyclePerNeuron',{cell(1,numNeurons)},...
                               'spMeanCycleRunOnPerNeuron',{cell(1,numNeurons)},...
                               ...
                               'meanPhaseVsC',{cell(1,numNeurons)},...
                               'stdPhaseVsC',{cell(1,numNeurons)},...
                               ...
                               'meanPhaseVsCRunOn',{cell(1,numNeurons)},...
                               'stdPhaseVsCRunOn',{cell(1,numNeurons)});      

for i = 1:numTrials
    if(isempty(thetaCycleStruct.indCycleStart{i}))
        continue;
    end
    thetaIndLfp = thetaCycleStruct.indCycleStart{i} - nSampBef;
    ind = find(abs(thetaIndLfp) == min(abs(thetaIndLfp)),1);
    spikeThetaPhaseStruct.thetaCycleAroundZero(i) = ind;
    spikeThetaPhaseStruct.indCycleStartAroundZero(i) = ...
        thetaCycleStruct.indCycleStart{i}(ind) - nSampBef;
    for j = 1:numNeurons
        % calculate the spike timing vs theta phase for each neuron
        spTimePerTrialPerNeuron = ...
            [spikesBef{j,indLapList(i)}' spikes{j,indLapList(i)}']; 
            % spike time of each neuron each trial
        spPhaseVsTPerTrialPerNeuron = ...
            [spikeThetaBef{j,indLapList(i)}' spikeTheta{j,indLapList(i)}']; 
    
        if(~isempty(thetaCycleStruct))
        % calculate the spike cycle vs theta phase for each neuron 
        % (***** cycle is not quite right when backward == 1 *********)
            for k = 1:thetaCycleStruct.numCyclesPerTrial(i)
                indSpikesInCycle = find(spTimePerTrialPerNeuron >= ...
                    thetaCycleStruct.indCycleStart{i}(k) - nSampBef & ...
                    spTimePerTrialPerNeuron < ...
                    thetaCycleStruct.indCycleStart{i}(k) - nSampBef ...
                    +thetaCycleStruct.cycleLenPerTrial{i}(k));
                
                if(~isempty(indSpikesInCycle))
                    spikeThetaPhaseStruct.spPhaseVsCPerTrialPerNeuron{j,i} = ...
                        [spikeThetaPhaseStruct.spPhaseVsCPerTrialPerNeuron{j,i}...
                        spPhaseVsTPerTrialPerNeuron(indSpikesInCycle)]; 
                        % spike theta phases of each neuron each trial 
                    spikeThetaPhaseStruct.spCyclePerTrialPerNeuron{j,i} =...
                        [spikeThetaPhaseStruct.spCyclePerTrialPerNeuron{j,i}...
                        k*ones(1,length(indSpikesInCycle))];
                    spikeThetaPhaseStruct.spCycleRunOnPerTrialPerNeuron{j,i} =...
                        [spikeThetaPhaseStruct.spCycleRunOnPerTrialPerNeuron{j,i}...
                        (k-ind)*ones(1,length(indSpikesInCycle))];
                    spikeThetaPhaseStruct.spMeanPhaseVsCPerTrialPerNeuron{j,i} = ...
                        [spikeThetaPhaseStruct.spMeanPhaseVsCPerTrialPerNeuron{j,i} ...
                        circ_mean(spPhaseVsTPerTrialPerNeuron(indSpikesInCycle),[],2)];
                    spikeThetaPhaseStruct.spMeanCyclePerTrialPerNeuron{j,i} = ...
                        [spikeThetaPhaseStruct.spMeanCyclePerTrialPerNeuron{j,i},k];
                    spikeThetaPhaseStruct.spMeanCycleRunOnPerTrialPerNeuron{j,i} = ...
                        [spikeThetaPhaseStruct.spMeanCycleRunOnPerTrialPerNeuron{j,i},k-ind];    
                end
            end
        end

        spikeThetaPhaseStruct.spTimePerNeuron{j} = ...
            [spikeThetaPhaseStruct.spTimePerNeuron{j} ...
            spTimePerTrialPerNeuron]; 
        spikeThetaPhaseStruct.spPhaseVsTPerNeuron{j} = ...
            [spikeThetaPhaseStruct.spPhaseVsTPerNeuron{j} ...
            spPhaseVsTPerTrialPerNeuron]; 
        spikeThetaPhaseStruct.spCyclePerNeuron{j} = ...
            [spikeThetaPhaseStruct.spCyclePerNeuron{j} ...
            spikeThetaPhaseStruct.spCyclePerTrialPerNeuron{j,i}]; 
            % spike theta cycles of each neuron over all valid trials
        spikeThetaPhaseStruct.spCycleRunOnPerNeuron{j} = ...
            [spikeThetaPhaseStruct.spCycleRunOnPerNeuron{j} ...
            spikeThetaPhaseStruct.spCycleRunOnPerTrialPerNeuron{j,i}]; 
            % spike theta cycles aligned by run onset of each neuron over all valid trials
        spikeThetaPhaseStruct.spPhaseVsCPerNeuron{j} = ...
            [spikeThetaPhaseStruct.spPhaseVsCPerNeuron{j} ...
            spikeThetaPhaseStruct.spPhaseVsCPerTrialPerNeuron{j,i}]; 
            % spike theta phases of each neuron over all valid trials
        spikeThetaPhaseStruct.spMeanPhaseVsCPerNeuron{j} = ...
            [spikeThetaPhaseStruct.spMeanPhaseVsCPerNeuron{j}...
            spikeThetaPhaseStruct.spMeanPhaseVsCPerTrialPerNeuron{j,i}]; 
            % mean spike phase per cycle for each neuron over all valid trials
        spikeThetaPhaseStruct.spMeanCyclePerNeuron{j} = ...
            [spikeThetaPhaseStruct.spMeanCyclePerNeuron{j}...
            spikeThetaPhaseStruct.spMeanCyclePerTrialPerNeuron{j,i}]; 
            % cycle number corresponding to spMeanPhaseVsCPerNeuron
        spikeThetaPhaseStruct.spMeanCycleRunOnPerNeuron{j} = ...
            [spikeThetaPhaseStruct.spMeanCycleRunOnPerNeuron{j}...
            spikeThetaPhaseStruct.spMeanCycleRunOnPerTrialPerNeuron{j,i}]; 
            % cycle number (aligned to run onset) corresponding to spMeanPhaseVsCPerNeuron
    end        
end

for i = 1:numNeurons
    numCycles = max(spikeThetaPhaseStruct.spMeanCyclePerNeuron{i});
    
    for j = 1:numCycles
        ind = spikeThetaPhaseStruct.spCyclePerNeuron{i} == j;
        if(sum(ind)>0)
            spikeThetaPhaseStruct.meanPhaseVsC{i}(j) = ...
                circ_mean(spikeThetaPhaseStruct.spPhaseVsCPerNeuron{i}(ind),[],2);
            spikeThetaPhaseStruct.stdPhaseVsC{i}(j) = ...
                circ_r(spikeThetaPhaseStruct.spPhaseVsCPerNeuron{i}(ind),[],[],2);
        else
            spikeThetaPhaseStruct.meanPhaseVsC{i}(j) = nan;
            spikeThetaPhaseStruct.stdPhaseVsC{i}(j) = nan;
        end
    end
    
    numCyclesMin = min(spikeThetaPhaseStruct.spMeanCycleRunOnPerNeuron{i});
    numCyclesMax = max(spikeThetaPhaseStruct.spMeanCycleRunOnPerNeuron{i});
    
    count = 1;
    for j = numCyclesMin:numCyclesMax
        ind = spikeThetaPhaseStruct.spCycleRunOnPerNeuron{i} == j;
        if(sum(ind)>0)
            spikeThetaPhaseStruct.meanPhaseVsCRunOn{i}(count) = ...
                circ_mean(spikeThetaPhaseStruct.spPhaseVsCPerNeuron{i}(ind),[],2);
            spikeThetaPhaseStruct.stdPhaseVsCRunOn{i}(count) = ...
                circ_r(spikeThetaPhaseStruct.spPhaseVsCPerNeuron{i}(ind),[],[],2);
        else
            spikeThetaPhaseStruct.meanPhaseVsCRunOn{i}(count) = nan;
            spikeThetaPhaseStruct.stdPhaseVsCRunOn{i}(count) = nan;
        end
        count = count+1;
    end
end
