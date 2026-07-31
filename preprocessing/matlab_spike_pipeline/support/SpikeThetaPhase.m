function spikeThetaPhaseStruct = SpikeThetaPhaseVR...
    (spikes, spikesDist, spikeTheta, thetaCycleStruct, indLapList, ...
     numNeurons, timeStep, figureState)
% calculate the theta phase of each spikes
% 1. vs sample number
% 2. vs theta cycle number
%
% Inputs: 
% spikes:               spike structure which contains a two dimensional cell
%                       structure, with spikes{i,j} representing the spike train
%                       from neuron i in trial j 
% spikesDist:           distance at which a spike fires
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
        disp('thetaCycleStruct should contain smaller or equal number of trials than spikeTheta');
        return;
    end
end

%%%% check whether trial no. in indLapList is beyond the range of trials
%%%% in the spike structure
maxTrialNo = max(indLapList);
if(maxTrialNo > sizeTrialsSp)
    disp('indLapList contains trials which is not included in the spikes structure');
    return;
end

numTrials = length(indLapList);
if(numTrials == 0)
    spikeThetaPhaseStruct = [];
    return;
end

spikeThetaPhaseStruct = struct('indLapList',indLapList,...
                               'spTimePerTrialPerNeuron',{cell(numNeurons,numTrials)},...
                               'spDistPerTrialPerNeuron',{cell(numNeurons,numTrials)},...
                               'spPhaseVsTPerTrialPerNeuron',{cell(numNeurons,numTrials)},...
                               ...
                               'spAllCyclePerTrialPerNeuron',{cell(numNeurons,numTrials)},...
                               'spAllPhaseVsCPerTrialPerNeuron',{cell(numNeurons,numTrials)},...
                               'spAllTimeVsCPerTrialPerNeuron',{cell(numNeurons,numTrials)},...
                               ...
                               'spTimePerNeuron',{cell(1,numNeurons)},...
                               'spDistPerNeuron',{cell(1,numNeurons)},...
                               'spPhaseVsTPerNeuron',{cell(1,numNeurons)},...
                               ...
                               'spCyclePerNeuron',{cell(1,numNeurons)},...
                               'spPhaseVsCPerNeuron',{cell(1,numNeurons)},...
                               'spTimeVsCPerNeuron',{cell(1,numNeurons)},...
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
    for j = 1:numNeurons
        if(~isempty(spikes{indLapList(i)}{j}))
            % calculate the spike timing vs theta phase for each neuron
            spikeThetaPhaseStruct.spTimePerTrialPerNeuron{j,i} = ...
                spikes{indLapList(i)}{j}; 
                % spike time of each neuron each trial
            spikeThetaPhaseStruct.spDistPerTrialPerNeuron{j,i} = ...
                spikesDist{indLapList(i)}{j}; 
                % spike distance of each neuron each trial
            spikeThetaPhaseStruct.spPhaseVsTPerTrialPerNeuron{j,i} = ...
                spikeTheta{indLapList(i)}{j}; 
                % spike phases of each neuron each trial

            spikeThetaPhaseStruct.spTimePerNeuron{j} = ...
                [spikeThetaPhaseStruct.spTimePerNeuron{j}; ...
                spikeThetaPhaseStruct.spTimePerTrialPerNeuron{j,i}]; 
                % spike time of each neuron over all valid trials
            spikeThetaPhaseStruct.spPhaseVsTPerNeuron{j} = ...
                [spikeThetaPhaseStruct.spPhaseVsTPerNeuron{j}; ...
                spikeThetaPhaseStruct.spPhaseVsTPerTrialPerNeuron{j,i}]; 
                % spike phases of each neuron over all valid trials
            spikeThetaPhaseStruct.spDistPerNeuron{j} = ...
                [spikeThetaPhaseStruct.spDistPerNeuron{j}; ...
                spikeThetaPhaseStruct.spDistPerTrialPerNeuron{j,i}]; 
                % spike distance of each neuron over all valid trials
            
            if(~isempty(thetaCycleStruct))
            % calculate the spike cycle vs theta phase for each neuron 
            % (***** cycle is not quite right when backward == 1 *********)
                for k = 1:thetaCycleStruct.numCyclesPerTrial(i)
                    indSpikesInCycle = find(spikes{indLapList(i)}{j} >= ...
                        thetaCycleStruct.indCycleStart{i}(k) & ...
                        spikes{indLapList(i)}{j} < ...
                        thetaCycleStruct.indCycleStart{i}(k)...
                        +thetaCycleStruct.cycleLenPerTrial{i}(k));
                    if(~isempty(indSpikesInCycle))
                        spikeThetaPhaseStruct.spAllCyclePerTrialPerNeuron{j,i} = ...
                            [spikeThetaPhaseStruct.spAllCyclePerTrialPerNeuron{j,i}; ...
                            k*ones(length(indSpikesInCycle),1)]; 
                            % spike theta cycles of each neuron each trial 
                        spikeThetaPhaseStruct.spAllPhaseVsCPerTrialPerNeuron{j,i} = ...
                            [spikeThetaPhaseStruct.spAllPhaseVsCPerTrialPerNeuron{j,i};...
                            spikeTheta{indLapList(i)}{j}(indSpikesInCycle)]; 
                            % spike theta phases of each neuron each trial 
                        spikeThetaPhaseStruct.spAllTimeVsCPerTrialPerNeuron{j,i} =...
                            [spikeThetaPhaseStruct.spAllTimeVsCPerTrialPerNeuron{j,i};...
                            spikes{indLapList(i)}{j}(indSpikesInCycle)];
                    
                        if(length(spikeThetaPhaseStruct.spAllCyclePerTrialPerNeuron{j,i})...
                           ~= length(spikeThetaPhaseStruct.spAllPhaseVsCPerTrialPerNeuron{j,i}))
                            disp('length error');
                        end
                    end
                end
            end

            spikeThetaPhaseStruct.spCyclePerNeuron{j} = ...
                [spikeThetaPhaseStruct.spCyclePerNeuron{j}; ...
                spikeThetaPhaseStruct.spAllCyclePerTrialPerNeuron{j,i}]; 
                % spike theta cycles of each neuron over all valid trials
            spikeThetaPhaseStruct.spPhaseVsCPerNeuron{j} = ...
                [spikeThetaPhaseStruct.spPhaseVsCPerNeuron{j}; ...
                spikeThetaPhaseStruct.spAllPhaseVsCPerTrialPerNeuron{j,i}]; 
                % spike theta phases of each neuron over all valid trials
            spikeThetaPhaseStruct.spTimeVsCPerNeuron{j} = ...
                [spikeThetaPhaseStruct.spTimeVsCPerNeuron{j};...
                spikeThetaPhaseStruct.spAllTimeVsCPerTrialPerNeuron{j,i}]; 
                % spike time of each neuron over all valid trials
        end
        
%       if(dispDebug == 1)
%           disp(['trial = ' num2str(i) ', neuron = ' num2str(j) ', cycle = ' num2str(k)]);
%       end
    end
end

        
%%%%%%%%%%% for each neuron, estimate the phase with the lowest spike
%%%%%%%%%%% activities (least preferred phase)\
%% changed by Yingxue on 4/7/2020
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
            circ_mean(spikeThetaPhaseStruct.spPhaseVsTPerNeuron{i});
        spikeThetaPhaseStruct.meanResultantLen(i) = ...
            circ_r(spikeThetaPhaseStruct.spPhaseVsTPerNeuron{i});
        spikeThetaPhaseStruct.pRayleigh(i) = ...
            circ_rtest(spikeThetaPhaseStruct.spPhaseVsTPerNeuron{i});
        spikeThetaPhaseStruct.pOmnibus(i) = ...
            circ_otest(spikeThetaPhaseStruct.spPhaseVsTPerNeuron{i});
        
        % test spike phase against gaussian distribution
    %             [h,pTheta(i)] = chi2gof(thetaSpikes{i});

        % construct spike trains which repeat the neuronal activity over
        % 3 phase cyclesq
        [thetaTimeTriple{i},thetaSpikesTriple{i}] = ...
            getMultiCycles(spikeThetaPhaseStruct.spTimePerNeuron{i}, ...
            spikeThetaPhaseStruct.spPhaseVsTPerNeuron{i}, 3);
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

if(figureState == 2)
    minTimeInterval = numSamples * timeStep;
    for i = 1:numNeurons
        if(mod(i-1,16) == 0)
            [figNew,pos] = CreateFig();
            set(0,'Units','pixels') 
            figTitle = 'Phase histogram';
            set(figure(figNew),'OuterPosition',...
                [pos(1) pos(2) pos(3)*2 pos(4)*2.2],'Name',figTitle)
        end

        subplot(4,4,mod(i-1,16)+1)
        fidsub1 = gca;

        % plot the histogram of theta phase
        if(~isempty(thetaSpikesTriple{i}))
            hist(thetaSpikesTriple{i},[stepPhase/2:stepPhase:1080-stepPhase/2]);
            h = findobj(gca,'Type','patch');
            set(h,'FaceColor','k','EdgeColor','k','FaceAlpha', 0.7);
        end
        hold on;
        if(spikeThetaPhaseStruct.indMinPhasePerNeuron(i)~= 0)
            h = plot([spikeThetaPhaseStruct.posPhase...
                            (spikeThetaPhaseStruct.indMinPhasePerNeuron(i)),...
                        spikeThetaPhaseStruct.posPhase...
                            (spikeThetaPhaseStruct.indMinPhasePerNeuron(i)+binNum)],...
                        [spikeThetaPhaseStruct.histPhasePerNeuron{i}...
                            (spikeThetaPhaseStruct.indMinPhasePerNeuron(i)),...
                        spikeThetaPhaseStruct.histPhasePerNeuron{i}...
                            (spikeThetaPhaseStruct.indMinPhasePerNeuron(i)+binNum)],...
                        'ro');
                % label the global minimum of the histogram
            set(h,'LineWidth',2.0);
        end
        
        if(~isempty(spikeThetaPhaseStruct.histPhasePerNeuron{i}))
            h = plot(stepPhase*(1/2-filDelay):stepPhase:1080-stepPhase*(1/2+filDelay),...
                spikeThetaPhaseStruct.histPhaseFilPerNeuron{i},'b-');
                % plot the filtered histogram
            set(h,'LineWidth',2.0);
        end
        if(spikeThetaPhaseStruct.indMinPhaseFilPerNeuron(i) ~= 0)
            h = plot([spikeThetaPhaseStruct.minPhaseFilArr(i) ...
                spikeThetaPhaseStruct.minPhaseFilArr(i)+360],...
                        [spikeThetaPhaseStruct.histPhaseFilPerNeuron{i}...
                            (spikeThetaPhaseStruct.indMinPhaseFilPerNeuron(i)) ...
                        spikeThetaPhaseStruct.histPhaseFilPerNeuron{i}
                            (spikeThetaPhaseStruct.indMinPhaseFilPerNeuron(i))],...
                        'm^');
                    % label the global minimum of the filtered histogram
            set(h,'LineWidth',2.0);
        end
        set(fidsub1,'FontSize',14.0,'Box','on','XLim',[0 1080]);
        if(mod(i,16) == 13)
            xlabel('Phase (Deg)')
            ylabel('Spike count')
        end
        title(['Neuron ' num2str(i)]);
    end
    pause;
    close all;

    for i = 1:numNeurons
        if(mod(i-1,16) == 0)
            [figNew,pos] = CreateFig();
            set(0,'Units','pixels') 
            figTitle = 'Phase vs time';
            set(figure(figNew),'OuterPosition',...
                [pos(1) pos(2) pos(3)*2 pos(4)*2.2],'Name',figTitle)
        end

        subplot(4,4,mod(i-1,16)+1)
        fidsub2 = gca;
        % plot the spike theta phase changes over time
        if(~isempty(thetaSpikesTriple{i}))
            h = plot(thetaTimeTriple{i}*timeStep,thetaSpikesTriple{i},'.');
            set(h,'LineWidth',2.0,'Color',[0 0 0]);
        end
        hold on;
        if(spikeThetaPhaseStruct.indMinPhasePerNeuron(i)~= 0)
            h = plot([min(thetaTimeTriple{i})*timeStep max(thetaTimeTriple{i})*timeStep],...
                [spikeThetaPhaseStruct.minPhaseArr(i) spikeThetaPhaseStruct.minPhaseArr(i)],...
                'r:');
            set(h,'LineWidth',1.5);
            h = plot([min(thetaTimeTriple{i})*timeStep max(thetaTimeTriple{i})*timeStep],...
                [spikeThetaPhaseStruct.minPhaseArr(i)+360 spikeThetaPhaseStruct.minPhaseArr(i)+360],...
                'r:'); 
            set(h,'LineWidth',1.5);
        
            % plot the boundary extracted using method 1 (histogram method)
            h = plot([min(thetaTimeTriple{i})*timeStep max(thetaTimeTriple{i})*timeStep],...
                [spikeThetaPhaseStruct.minPhaseFilArr(i) spikeThetaPhaseStruct.minPhaseFilArr(i)],...
                'b:');
            set(h,'LineWidth',1.5);
            h = plot([min(thetaTimeTriple{i})*timeStep max(thetaTimeTriple{i})*timeStep],...
                [spikeThetaPhaseStruct.minPhaseFilArr(i)+360 spikeThetaPhaseStruct.minPhaseFilArr(i)+360],...
                'b:'); 
            set(h,'LineWidth',1.5);
        end
            % plot the boundary extracted using method 2 (filtered
            % histogram method)
        set(fidsub2,'FontSize',14.0,'Box','on',...
            'XLim',[0 minTimeInterval],'YLim', [0 1080]); 
        title(['Neuron' num2str(i)])
        if(mod(i,16) == 13)
            xlabel('Time (sec)')
            ylabel('Spike phase (Deg)')
        end
    end
    pause;
    close all;
    
    for i = 1:numNeurons
        if(mod(i-1,16) == 0)
            [figNew,pos] = CreateFig();
            set(0,'Units','pixels') 
            figTitle = 'Phase vs cycle';
            set(figure(figNew),'OuterPosition',...
                [pos(1) pos(2) pos(3)*2 pos(4)*2.2],'Name',figTitle)
        end

        subplot(4,4,mod(i-1,16)+1)
        fidsub3 = gca;
        if(~isempty(spikeThetaPhaseStruct.spCyclePerNeuron{i}))
            [thetaCyclesTriple,thetaSpikesTriple] = ...
                getMultiCycles(spikeThetaPhaseStruct.spCyclePerNeuron{i},...
                    spikeThetaPhaseStruct.spPhaseVsCPerNeuron{i}, 3);
            % plot the spike theta phase changes over time
            h = plot(thetaCyclesTriple,thetaSpikesTriple,'.');
            set(h,'LineWidth',2.0,'Color',[0 0 0]);
            set(fidsub3,'FontSize',14.0,'Box','on','XLim',...
                [0 max(thetaCyclesTriple)],'YLim', [0 1080]); 
        end
        title(['Neuron' num2str(i)])
        if(mod(i,16) == 13)
            xlabel('Cum theta (rad)')
            ylabel('Spike phase (Deg)')
        end
    end
    pause;
    close all;
end