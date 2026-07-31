function spikeThetaPhaseStruct = SpikeThetaPhaseAlignedRun...
    (spikes, spikeTheta, thetaCycleStruct, indLapList, ...
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
                               'spPhaseVsCPerTrialPerNeuron',{cell(numNeurons,numTrials)},...
                               'spCyclePerTrialPerNeuron',{cell(numNeurons,numTrials)},...
                               'spCycleRunOnPerTrialPerNeuron',{cell(numNeurons,numTrials)},...
                               'spMeanPhaseVsCPerTrialPerNeuron',{cell(numNeurons,numTrials)},...
                               'spMeanCyclePerTrialPerNeuron',{cell(numNeurons,numTrials)},...
                               'spMeanCycleRunOnPerTrialPerNeuron',{cell(numNeurons,numTrials)},...
                               ...
                               'spTimePerNeuron',{cell(1,numNeurons)},...
                               'spPhaseVsTPerNeuron',{cell(1,numNeurons)},...
                               'spCyclePerNeuron',{cell(1,numNeurons)},...
                               ...
                               'spPhaseVsCPerNeuron',{cell(1,numNeurons)},...
                               'spMeanPhaseVsCPerNeuron',{cell(1,numNeurons)},...
                               'spMeanCyclePerNeuron',{cell(1,numNeurons)},...
                               ...
                               'meanPhaseVsC',{cell(1,numNeurons)},...
                               'stdPhaseVsC',{cell(1,numNeurons)},...
                               ...
                               'histPhasePerNeuron',{cell(1,numNeurons)},...      % histogram of the spike phases
                               'posPhase',[],...                % index of each histogram bar
                               'minPhaseArr',-1*ones(1,numNeurons),... % the phase which fires the least number of spikes 
                               'maxPhaseArr',-1*ones(1,numNeurons),... % the phase which fires the largest number of spikes 
                               'indMinPhasePerNeuron',zeros(1,numNeurons),... % the index of the least firing phase in the histogram
                               'indMaxPhasePerNeuron',zeros(1,numNeurons),... % the index of the max firing phase in the histogram 
                               'histPhaseFilPerNeuron',[],... % the filtered histogram
                               'indMinPhaseFilPerNeuron',zeros(1,numNeurons),... % index of the min
                               'indMaxPhaseFilPerNeuron',zeros(1,numNeurons),... % index of the max
                               'minPhaseFilArr',-1*ones(1,numNeurons),...   % the min phase of the filtered histogram
                               'maxPhaseFilArr',-1*ones(1,numNeurons),...   % the max phase of the filtered histogram
                               'thetaMod',-1*ones(1,numNeurons),... % theta modulation calculated based on filtered histogram
                               ...
                               'meanDire', zeros(1,numNeurons),... % the mean phase direction
                               'meanResultantLen', zeros(1,numNeurons),... % the mean resultant length of the mean phase direction
                               'pRayleigh', zeros(1,numNeurons),... % p-value obtained from Rayleigh test
                               'pOmnibus', zeros(1,numNeurons)); % p-value obtained from Omnibus test       

for i = 1:numTrials
    if(isempty(thetaCycleStruct.indCycleStart{i}))
        continue;
    end
    for j = 1:numNeurons
        % calculate the spike timing vs theta phase for each neuron
        spTimePerTrialPerNeuron = ...
            [spikes{j,indLapList(i)}']; 
            % spike time of each neuron each trial
        spPhaseVsTPerTrialPerNeuron = ...
            [spikeTheta{j,indLapList(i)}']; 
    
        if(~isempty(thetaCycleStruct))
        % calculate the spike cycle vs theta phase for each neuron 
        % (***** cycle is not quite right when backward == 1 *********)
            for k = 1:thetaCycleStruct.numCyclesPerTrial(i)
                indSpikesInCycle = find(spTimePerTrialPerNeuron >= ...
                    thetaCycleStruct.indCycleStart{i}(k) & ...
                    spTimePerTrialPerNeuron < ...
                    thetaCycleStruct.indCycleStart{i}(k) ...
                    +thetaCycleStruct.cycleLenPerTrial{i}(k));
                
                if(~isempty(indSpikesInCycle))
                    spikeThetaPhaseStruct.spPhaseVsCPerTrialPerNeuron{j,i} = ...
                        [spikeThetaPhaseStruct.spPhaseVsCPerTrialPerNeuron{j,i}...
                        spPhaseVsTPerTrialPerNeuron(indSpikesInCycle)]; 
                        % spike theta phases of each neuron each trial 
                    spikeThetaPhaseStruct.spCyclePerTrialPerNeuron{j,i} =...
                        [spikeThetaPhaseStruct.spCyclePerTrialPerNeuron{j,i}...
                        k*ones(1,length(indSpikesInCycle))];
                    spikeThetaPhaseStruct.spMeanPhaseVsCPerTrialPerNeuron{j,i} = ...
                        [spikeThetaPhaseStruct.spMeanPhaseVsCPerTrialPerNeuron{j,i} ...
                        circ_mean(spPhaseVsTPerTrialPerNeuron(indSpikesInCycle),[],2)];
                    spikeThetaPhaseStruct.spMeanCyclePerTrialPerNeuron{j,i} = ...
                        [spikeThetaPhaseStruct.spMeanCyclePerTrialPerNeuron{j,i},k]; 
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
end


%%%%%%%%%%% for each neuron, estimate the phase with the lowest spike
%%%%%%%%%%% activities (least preferred phase)\
%% changed by Yingxue on 4/28/2020
[b,a] = butter(2,0.1); % low pass filtered phase histogram
%%
filDelay = 0; % the emperical delay of the filter
for i = 1:numNeurons
    stepPhase = 5; %(degree)
    stepPhaseDefault = 5;
    if(mod(360,stepPhase) ~= 0) 
        % if 360 is not divisible by the stepPhase, then use the default stepPhase
        stepPhase = stepPhaseDefault;
    end
    binNum = 360/stepPhase;
    if(~isempty(spikeThetaPhaseStruct.spTimePerNeuron{i}))
        spikeThetaPhaseStruct.meanDire(i) = ...
            circ_mean(spikeThetaPhaseStruct.spPhaseVsTPerNeuron{i}');
        spikeThetaPhaseStruct.meanResultantLen(i) = ...
            circ_r(spikeThetaPhaseStruct.spPhaseVsTPerNeuron{i}');
        spikeThetaPhaseStruct.pRayleigh(i) = ...
            circ_rtest(spikeThetaPhaseStruct.spPhaseVsTPerNeuron{i}');
        spikeThetaPhaseStruct.pOmnibus(i) = ...
            circ_otest(spikeThetaPhaseStruct.spPhaseVsTPerNeuron{i}');
        
        % test spike phase against gaussian distribution
    %             [h,pTheta(i)] = chi2gof(thetaSpikes{i});

        % construct spike trains which repeat the neuronal activity over
        % 3 phase cyclesq
        [thetaTimeTriple{i},thetaSpikesTriple{i}] = ...
            getMultiCycles(spikeThetaPhaseStruct.spTimePerNeuron{i}', ...
            spikeThetaPhaseStruct.spPhaseVsTPerNeuron{i}', 3);
            % three cycles are used because in the second method of min 
            % estimation, the left most part of the filtered
            % histogram is influenced by the initial condition of the 
            % filtering process, thus is different from the filtered histograms of 
            % the second and third cycle. Therefore, the left most part of 
            % the filtered histogram is not used in the analysis 

        %%%% method1: plot the phase histogram and find the phase with
        %%%% minimum and maximum activity
        % min
        [spikeThetaPhaseStruct.histPhasePerNeuron{i},spikeThetaPhaseStruct.posPhase] =...
            hist(thetaSpikesTriple{i},[stepPhase/2:stepPhase:1080-stepPhase/2]);
        indMinPhase = ...
            find(spikeThetaPhaseStruct.histPhasePerNeuron{i}(binNum+1:2*binNum)...
            == min(spikeThetaPhaseStruct.histPhasePerNeuron{i}(binNum+1:2*binNum))); 
            % find the min spike count within a cycle
        if(length(indMinPhase)> 1) 
            % if there are several min elements, select the first whose 
            % closest neighbor also reaches the min value
            diffIndex = diff(indMinPhase);
            minDiffIndex = find(diffIndex == min(diffIndex));
            indMinPhase = indMinPhase(minDiffIndex(1));
        end
        spikeThetaPhaseStruct.minPhaseArr(i) = ...
            spikeThetaPhaseStruct.posPhase(indMinPhase); 
            % the starting phase of each neuron
        spikeThetaPhaseStruct.indMinPhasePerNeuron(i) = indMinPhase; 
            % the index of the min in the histogram
       
        % max
        indMaxPhase = ...
            find(spikeThetaPhaseStruct.histPhasePerNeuron{i}(binNum+1:2*binNum)...
            == max(spikeThetaPhaseStruct.histPhasePerNeuron{i}(binNum+1:2*binNum)));
            % find the min spike count within a cycle
        if(length(indMaxPhase)> 1) 
            % if there are several min elements, select the first whose 
            % closest neighbor also reaches the min value
            diffIndex = diff(indMaxPhase);
            minDiffIndex = find(diffIndex == min(diffIndex));
            indMaxPhase = indMaxPhase(minDiffIndex(1));
        end
        spikeThetaPhaseStruct.maxPhaseArr(i) = ...
            spikeThetaPhaseStruct.posPhase(indMaxPhase); 
            % the starting phase of each neuron
        spikeThetaPhaseStruct.indMaxPhasePerNeuron(i) = indMaxPhase; 
            % the index of the min in the histogram

        %%%% method2: use filtered histogram and compute the phase of
        %%%% the global minimum and maximum
        %%% changed by Yingxue from using "filter" to "filtfilt" on
        %%% 4/7/2020
        [spikeThetaPhaseStruct.histPhaseFilPerNeuron{i}] = ...
            filtfilt(b,a,spikeThetaPhaseStruct.histPhasePerNeuron{i}); 
            % filter with initial condition to guarantee that the first
            % output value = the first input value
        % max  
        peakVal = max(spikeThetaPhaseStruct.histPhaseFilPerNeuron{i}(binNum+1:2*binNum));
        peakInd = ...
            find(spikeThetaPhaseStruct.histPhaseFilPerNeuron{i}(binNum+1:2*binNum)...
            == peakVal); 
        % find the max of the filtered histogram
%         if(peakInd(1) < 2*binNum) 
%             % use the second and third cycles of the filtered histogram
%             peakInd1 = peakInd(1);
%             peakInd2 = peakInd(1) + binNum;
%         else
%             peakInd1 = peakInd(1) - binNum;
%             peakInd2 = peakInd(1);
%         end
                
        indMaxPhaseFil = peakInd;
        spikeThetaPhaseStruct.maxPhaseFilArr(i) = ...
            (indMaxPhaseFil(1) - 1/2 - filDelay)*stepPhase;
                    % find the phase value: since the filter has a
                    % delay by itself, the index in the
                    % original histogram should be the index of max
                    % value - the filter delay. Since the phase
                    % at the middle of the phase bin is considered, 
                    % another 1/2 is substracted from the index
        spikeThetaPhaseStruct.indMaxPhaseFilPerNeuron(i) = ...
            indMaxPhaseFil(1) - filDelay;
        if(spikeThetaPhaseStruct.maxPhaseFilArr(i) > 360) 
            % the starting phase should be smaller than 360 degree
            spikeThetaPhaseStruct.maxPhaseFilArr(i) = ...
                mod(spikeThetaPhaseStruct.maxPhaseFilArr(i),360);
        end
           
        % min
        minVal = min(spikeThetaPhaseStruct.histPhaseFilPerNeuron{i}(binNum+1:2*binNum));
        indMinPhaseFil = ...
            find(spikeThetaPhaseStruct.histPhaseFilPerNeuron{i}(binNum+1:2*binNum)...
            == minVal);
                    % find the min within a cycle
        spikeThetaPhaseStruct.minPhaseFilArr(i) = ...
            (indMinPhaseFil(1) - 1/2 - filDelay)*stepPhase;
                    % find the phase value: since the filter has a
                    % delay by itself, the index in the
                    % original histogram should be the index of min
                    % value - the filter delay. Since the phase
                    % at the middle of the phase bin is considered, 
                    % another 1/2 is substracted from the index
        spikeThetaPhaseStruct.indMinPhaseFilPerNeuron(i) = ...
            indMinPhaseFil(1) - filDelay;
        if(spikeThetaPhaseStruct.minPhaseFilArr(i) > 360) 
            % the starting phase should be smaller than 360 degree
            spikeThetaPhaseStruct.minPhaseFilArr(i) = ...
                mod(spikeThetaPhaseStruct.minPhaseFilArr(i),360);
        end
        
        %% added by Yingxue on 4/7/2020
        spikeThetaPhaseStruct.thetaMod(i) = (peakVal - minVal)./(peakVal + minVal);
        %%
        
%         spikeThetaPhaseStruct.meanCount(i) = mean(spikeThetaPhaseStruct.histPhasePerNeuron{i});
%         indLg = spikeThetaPhaseStruct.histPhasePerNeuron{i}>=spikeThetaPhaseStruct.meanCount(i);
%         indSm = spikeThetaPhaseStruct.histPhasePerNeuron{i}<spikeThetaPhaseStruct.meanCount(i);
%         histPhasePerNeuronMod = spikeThetaPhaseStruct.histPhasePerNeuron{i};
%         histPhasePerNeuronMod(indLg) = histPhasePerNeuronMod(indLg) - spikeThetaPhaseStruct.meanCount(i);
%         histPhasePerNeuronMod(indSm) = 0;
%         spikeThetaPhaseStruct.meanDireVar(i) = circ_mean(histPhasePerNeuronMod');
%         spikeThetaPhaseStruct.meanResultantLenVar(i) = circ_r(histPhasePerNeuronMod');
%         spikeThetaPhaseStruct.pRayleighVar(i) = circ_rtest(histPhasePerNeuronMod');
%         spikeThetaPhaseStruct.pOmnibusVar(i) = circ_otest(histPhasePerNeuronMod');
    else
        thetaTimeTriple{i} = [];
        thetaSpikesTriple{i} = [];
        spikeThetaPhaseStruct.histPhasePerNeuron{i} = [];
        spikeThetaPhaseStruct.histPhaseFilPerNeuron{i} = [];
        disp(['Neuron ' num2str(i) ' does not have any spikes within this time interval.']);
    end
end
