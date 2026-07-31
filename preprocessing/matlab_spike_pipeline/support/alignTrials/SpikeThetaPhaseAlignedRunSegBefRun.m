function spikeThetaPhaseSegStruct = SpikeThetaPhaseAlignedRunSegBefRun...
    (spikeThetaPhaseStruct, timeInt,sampleFq)
 
    numNeurons = length(spikeThetaPhaseStruct.spTimePerNeuron);
    
    spikeThetaPhaseSegStruct = struct(...
                                       'indLapList',spikeThetaPhaseStruct.indLapList,...
                                       'timeInterval',timeInt,...
                                       ...
                                       'spTimePerNeuron', {cell(1,numNeurons)},... 
                                       'spPhaseVsTPerNeuron',{cell(1,numNeurons)},...
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
                                       'pOmnibus', zeros(1,numNeurons));
      
     
     %% changed by Yingxue on 4/28/2020
    [b,a] = butter(2,0.1); % low pass filtered phase histogram
    %%
    filDelay = 0; % the emperical delay of the filter
    for i = 1:numNeurons
        ind = spikeThetaPhaseStruct.spTimePerNeuron{i} >= timeInt(1)*sampleFq &...
             spikeThetaPhaseStruct.spTimePerNeuron{i} < timeInt(2)*sampleFq;
        spikeThetaPhaseSegStruct.spTimePerNeuron{i} = spikeThetaPhaseStruct.spTimePerNeuron{i}(ind);
        spikeThetaPhaseSegStruct.spPhaseVsTPerNeuron{i} = spikeThetaPhaseStruct.spPhaseVsTPerNeuron{i}(ind);
        
        stepPhase = 5; %(degree)
        stepPhaseDefault = 5;
        if(mod(360,stepPhase) ~= 0) 
            % if 360 is not divisible by the stepPhase, then use the default stepPhase
            stepPhase = stepPhaseDefault;
        end
        binNum = 360/stepPhase;
        if(~isempty(spikeThetaPhaseSegStruct.spTimePerNeuron{i}))
            spikeThetaPhaseSegStruct.meanDire(i) = ...
                circ_mean(spikeThetaPhaseSegStruct.spPhaseVsTPerNeuron{i}');
            spikeThetaPhaseSegStruct.meanResultantLen(i) = ...
                circ_r(spikeThetaPhaseSegStruct.spPhaseVsTPerNeuron{i}');
            spikeThetaPhaseSegStruct.pRayleigh(i) = ...
                circ_rtest(spikeThetaPhaseSegStruct.spPhaseVsTPerNeuron{i}');
            spikeThetaPhaseSegStruct.pOmnibus(i) = ...
                circ_otest(spikeThetaPhaseSegStruct.spPhaseVsTPerNeuron{i}');

            % construct spike trains which repeat the neuronal activity over
            % 3 phase cyclesq
            [thetaTimeTriple{i},thetaSpikesTriple{i}] = ...
                getMultiCycles(spikeThetaPhaseSegStruct.spTimePerNeuron{i}', ...
                spikeThetaPhaseSegStruct.spPhaseVsTPerNeuron{i}', 3);
                % three cycles are used because in the second method of min 
                % estimation, the left most part of the filtered
                % histogram is influenced by the initial condition of the 
                % filtering process, thus is different from the filtered histograms of 
                % the second and third cycle. Therefore, the left most part of 
                % the filtered histogram is not used in the analysis 

            %%%% method1: plot the phase histogram and find the phase with
            %%%% minimum and maximum activity
            % min
            [spikeThetaPhaseSegStruct.histPhasePerNeuron{i},spikeThetaPhaseSegStruct.posPhase] =...
                hist(thetaSpikesTriple{i},[stepPhase/2:stepPhase:1080-stepPhase/2]);
            indMinPhase = ...
                find(spikeThetaPhaseSegStruct.histPhasePerNeuron{i}(binNum+1:2*binNum)...
                == min(spikeThetaPhaseSegStruct.histPhasePerNeuron{i}(binNum+1:2*binNum))); 
                % find the min spike count within a cycle
            if(length(indMinPhase)> 1) 
                % if there are several min elements, select the first whose 
                % closest neighbor also reaches the min value
                diffIndex = diff(indMinPhase);
                minDiffIndex = find(diffIndex == min(diffIndex));
                indMinPhase = indMinPhase(minDiffIndex(1));
            end
            spikeThetaPhaseSegStruct.minPhaseArr(i) = ...
                spikeThetaPhaseSegStruct.posPhase(indMinPhase); 
                % the starting phase of each neuron
            spikeThetaPhaseSegStruct.indMinPhasePerNeuron(i) = indMinPhase; 
                % the index of the min in the histogram

            % max
            indMaxPhase = ...
                find(spikeThetaPhaseSegStruct.histPhasePerNeuron{i}(binNum+1:2*binNum)...
                == max(spikeThetaPhaseSegStruct.histPhasePerNeuron{i}(binNum+1:2*binNum)));
                % find the min spike count within a cycle
            if(length(indMaxPhase)> 1) 
                % if there are several min elements, select the first whose 
                % closest neighbor also reaches the min value
                diffIndex = diff(indMaxPhase);
                minDiffIndex = find(diffIndex == min(diffIndex));
                indMaxPhase = indMaxPhase(minDiffIndex(1));
            end
            spikeThetaPhaseSegStruct.maxPhaseArr(i) = ...
                spikeThetaPhaseSegStruct.posPhase(indMaxPhase); 
                % the starting phase of each neuron
            spikeThetaPhaseSegStruct.indMaxPhasePerNeuron(i) = indMaxPhase; 
                % the index of the min in the histogram

            %%%% method2: use filtered histogram and compute the phase of
            %%%% the global minimum and maximum
            %%% changed by Yingxue from using "filter" to "filtfilt" on
            %%% 4/7/2020
            [spikeThetaPhaseSegStruct.histPhaseFilPerNeuron{i}] = ...
                filtfilt(b,a,spikeThetaPhaseSegStruct.histPhasePerNeuron{i}); 
                % filter with initial condition to guarantee that the first
                % output value = the first input value
            % max  
            peakVal = max(spikeThetaPhaseSegStruct.histPhaseFilPerNeuron{i}(binNum+1:2*binNum));
            peakInd = ...
                find(spikeThetaPhaseSegStruct.histPhaseFilPerNeuron{i}(binNum+1:2*binNum)...
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
            spikeThetaPhaseSegStruct.maxPhaseFilArr(i) = ...
                (indMaxPhaseFil(1) - 1/2 - filDelay)*stepPhase;
                        % find the phase value: since the filter has a
                        % delay by itself, the index in the
                        % original histogram should be the index of max
                        % value - the filter delay. Since the phase
                        % at the middle of the phase bin is considered, 
                        % another 1/2 is substracted from the index
            spikeThetaPhaseSegStruct.indMaxPhaseFilPerNeuron(i) = ...
                indMaxPhaseFil(1) - filDelay;
            if(spikeThetaPhaseSegStruct.maxPhaseFilArr(i) > 360) 
                % the starting phase should be smaller than 360 degree
                spikeThetaPhaseSegStruct.maxPhaseFilArr(i) = ...
                    mod(spikeThetaPhaseSegStruct.maxPhaseFilArr(i),360);
            end

            % min
            minVal = min(spikeThetaPhaseSegStruct.histPhaseFilPerNeuron{i}(binNum+1:2*binNum));
            indMinPhaseFil = ...
                find(spikeThetaPhaseSegStruct.histPhaseFilPerNeuron{i}(binNum+1:2*binNum)...
                == minVal);
                        % find the min within a cycle
            spikeThetaPhaseSegStruct.minPhaseFilArr(i) = ...
                (indMinPhaseFil(1) - 1/2 - filDelay)*stepPhase;
                        % find the phase value: since the filter has a
                        % delay by itself, the index in the
                        % original histogram should be the index of min
                        % value - the filter delay. Since the phase
                        % at the middle of the phase bin is considered, 
                        % another 1/2 is substracted from the index
            spikeThetaPhaseSegStruct.indMinPhaseFilPerNeuron(i) = ...
                indMinPhaseFil(1) - filDelay;
            if(spikeThetaPhaseSegStruct.minPhaseFilArr(i) > 360) 
                % the starting phase should be smaller than 360 degree
                spikeThetaPhaseSegStruct.minPhaseFilArr(i) = ...
                    mod(spikeThetaPhaseSegStruct.minPhaseFilArr(i),360);
            end

            %% added by Yingxue on 4/7/2020
            spikeThetaPhaseSegStruct.thetaMod(i) = (peakVal - minVal)./(peakVal + minVal);
            
        else
            thetaTimeTriple{i} = [];
            thetaSpikesTriple{i} = [];
            spikeThetaPhaseSegStruct.histPhasePerNeuron{i} = [];
            spikeThetaPhaseSegStruct.histPhaseFilPerNeuron{i} = [];
            disp(['Neuron ' num2str(i) ' does not have any spikes within this time interval.']);
        end 
    end
