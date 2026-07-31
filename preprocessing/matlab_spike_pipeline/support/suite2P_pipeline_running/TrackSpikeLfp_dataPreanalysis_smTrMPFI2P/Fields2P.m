function fieldStruct = Fields2P(pFRStruct,indNeurons,paramFW,paramGen,figureState, numTrials)
% extract episode fields from the neuron activity
% 
% Inputs:
% pFRStruct:                peak firing rate structure (referring to function
%                           PeakFR)
% spikeThetaPhaseStruct:    spike theta phase structure (referring to function SpikeThetaPhase)
% meanDist:                 the distance the animal has travelled, this is
%                           averaged over all the trials, this entry should
%                           be empty if the max trial length is used to
%                           extract fields
% indNeurons:               indices of neurons
% paramFW:                  field width parameters (referring to FieldWidth3/FieldWidthLR)
% paramGen:                 general parameters, including numSamples,
%                           spaceBin, mahaldist and mahaldistThre
% figureState:              0: figure off
%                           2: figure on
% backward:             align the trials from the end
% lenTrials:            the actual length of each trial listed in indLapList
%
% Output: fieldStruct with the following fields
% indNeuron:                indices of neurons where a field occurs
% FW:                       field width (unit: samples)
% indStartField:            start index of a field
% FWDistance:               field width (unit: distance)
% numSpPerField:            number of spikes within a field
% indPeakField:             index of the peak firing rate within a field
% p2MFRField:               peak instantaneous firing rate within a field
%                           over mean instantaneous firing rate
% p2MinFRField:             peak instantaneous firing rate within a field
%                           over min instantaneous firing rate
% fDoubleField:             1: double field                        
%                           0: single field
% distSpikesPerField:       time of each spike within each field
% phaseSpikesPerField:      phase of each spike within each field
% neuronDoubleField:        neurons contain double fields
% neuronSingleField:        neurons contain single field
% neuronFieldAll:           neurons contain at least one field
% neuronNoField:            neurons without any field

if(isempty(pFRStruct))
    fieldStruct = [];
    return;
end

numNeurons = length(indNeurons);
if(numNeurons == 0)
   fieldStruct = [];
   return;
end

fieldStruct = struct(...
                'indNeuron',[],... % neuron index
                'numTrials', numTrials,... % number of trials
                'FW',[],... % field width (samples)
                'indStartField',[],... % the start index of the field
                'numSpPerField',[],... % number of spikes per field
                'indPeakField',[],... % index of peak firing rate within the field
                'peakLocInField',[],... % the location of the peak with respect to the whole field
                'peakInstFiringRate',[],... % peak instantaneous firing rate
                'meanInstFRArr',[],... % mean inst firing rate array
                'p2MFRField',[],... % peak/mean inst FR ratio
                'fDoubleField',[],... % single/double field (0/1) 
                'distSpikesPerField',[],... % time of each spike within each field
                'phaseSpikesPerField',[],... % phase of each spike within each field
                'numSpikesPerTrialPerField',[],... % number of spikes per trial in each field
                'percSpikesPerField',[],... % percentage of spikes per field
                'numTrialsPerField',[],... % number of active trials in each field
                'neuronDoubleField',[],... % neurons contain double fields
                'neuronSingleField',[],... % neurons contain single field
                'neuronFieldAll',[],... % neurons contain at least one field
                'neuronNoField',[],...
                'percActiveTrials',paramGen.percActiveTrials); % perc. of active trials


%%%%%%%%% estimate field width (only used unnormalized peak firing
%%%%%%%%% rate)
numFieldTotal = 1; % total number of fields  
meanInstFR = pFRStruct.meanInstFR;
avgFRProfile = pFRStruct.avgFRProfile;
% spDistPerNeuron = spikeThetaPhaseStruct.spDistPerNeuron;
% spDistPerTrialPerNeuron = spikeThetaPhaseStruct.spDistPerTrialPerNeuron;

for i = 1:numNeurons
    if(indNeurons(i) == 125 || indNeurons(i) == 48 || indNeurons(i) == 152)
       a = 1; 
    end
    
    % set parameters
    if(meanInstFR(indNeurons(i)) < paramFW.highInstFR)
        paramFW.minNumSpikesPerField = ...
            paramFW.minNumSpikesPerFieldPerTr*fieldStruct.numTrials;
        paramFW.minSpikeDenPerFieldPerHzTmp = ...
            paramFW.minSpikeDenPerFieldPerHzPerTr;
        paramFW.percSpikesInFieldsTmp = paramFW.percSpikesInFields;
        paramFW.percSpikesInFieldsPerUnitTmp = paramFW.percSpikesInFieldsPerUnit;
        paramFW.maxPercSpikesInFieldsTmp = paramFW.maxPercSpikesInFields;
    else
        paramFW.minNumSpikesPerField = ...
            paramFW.minNumSpikesPerFieldPerTrHighFR*fieldStruct.numTrials;
        paramFW.minSpikeDenPerFieldPerHzTmp =...
            paramFW.minSpikeDenPerFieldPerHzPerTrHighFR;
        paramFW.percSpikesInFieldsTmp = paramFW.percSpikesInFieldsHigh;
        paramFW.percSpikesInFieldsPerUnitTmp = paramFW.percSpikesInFieldsPerUnitHigh;
        paramFW.maxPercSpikesInFieldsTmp = paramFW.maxPercSpikesInFieldsHigh;
    end
    paramFW.minNumSpikesPerField = max(paramFW.minNumSpikesPerField, ...
        paramFW.minNumSpikesPerFieldPerHz*pFRStruct.peakFR(indNeurons(i)));
                % calculate the min number of spikes per field based on the
                % peak firing rate

    if(pFRStruct.nonZeroNumTrials(indNeurons(i)) < ...
            floor(paramGen.percActiveTrials*fieldStruct.numTrials))
            % guarantee enough number of trials has elicited spikes
        disp(['Neuron ' num2str(indNeurons(i)) ': only ' ...
            num2str(pFRStruct.nonZeroNumTrials(indNeurons(i))) ...
            ' trials have non-zero number of spikes. Threshold: ' ...
            num2str(floor(paramGen.percActiveTrials*fieldStruct.numTrials))]);
        continue;
    end

    % extract the potential field areas based on the mean firing rate of
    % neurons (neurons with high firing rate tend to have smaller peak to 
    % mean firing rate ratio)
    if(meanInstFR(indNeurons(i)) < paramFW.highInstFR) 
        indexHighFR = find(avgFRProfile(indNeurons(i),:) >= ...
            meanInstFR(indNeurons(i))*paramFW.pToMFRRatioFieldMinLowFR); 
    else
        indexHighFR = find(avgFRProfile(indNeurons(i),:) >= ...
            meanInstFR(indNeurons(i))*paramFW.pToMFRRatioFieldMinHighFR); 
    end
                        
    indexLowFRLT = find(avgFRProfile(indNeurons(i),:) <= ...
        meanInstFR(indNeurons(i))*paramFW.lowThreFieldMeanInstFR); 
        % guarantee the firing rate falls back to 0
    indexLowFRHT = find(avgFRProfile(indNeurons(i),:) <= ...
        meanInstFR(indNeurons(i))*paramFW.highThreFieldMeanInstFR);
        % guarantee the firing rate falls back to 0 (double thresholds)

    if(~isempty(indexHighFR) && ~isempty(indexLowFRLT))
        [indFieldBoundTmp, numFieldsTmp] = ...
            SearchContRegionsDoubleThre(...
            indexHighFR,indexLowFRLT,indexLowFRHT,avgFRProfile(indNeurons(i),:),...
            paramGen.numSamples,paramFW.MaxFringeLen/paramGen.spaceMergeBin); 
            % find the potential fields
        numFieldTotalCurNeuronStart = numFieldTotal;
        
        % check whether the field is invading the other potential fields
        if(numFieldsTmp > 1)
            validInd = zeros(1,numFieldsTmp);
            indFieldBoundStartNextTmp = ...
                [indFieldBoundTmp(:,1);length(avgFRProfile(indNeurons(i),:))];
            indFieldBoundStopPrevTmp = [1;indFieldBoundTmp(:,2)];
            for k = 1:numFieldsTmp
                frArrField = ...
                    avgFRProfile(indNeurons(i),indFieldBoundTmp(k,1):indFieldBoundTmp(k,2)); 
                    % the firing profile within each potential field
                maxIndField = ...
                    find(frArrField == max(frArrField),1) + indFieldBoundTmp(k,1) - 1; 
                    % index of the max firing rate within the potential field
                if(maxIndField - paramFW.FieldMaxLen/paramGen.spaceMergeBin/2 > ...
                        indFieldBoundStopPrevTmp(k))  
                        % the index of the starting point of a field should
                        % not be smaller than max index - half of the max 
                        % field length
                    indFieldBoundStopPrevTmp(k) = floor(maxIndField - ...
                        paramFW.FieldMaxLen/paramGen.spaceMergeBin/2);
                end
                if(maxIndField + paramFW.FieldMaxLen/paramGen.spaceMergeBin/2 <...
                        indFieldBoundStartNextTmp(k+1))
                    indFieldBoundStartNextTmp(k+1) = ceil(maxIndField + ...
                        paramFW.FieldMaxLen/paramGen.spaceMergeBin/2); 
                        % the index of the stopping point of a field should
                        % not be larger than max index + half of the max 
                        % field length
                end
                
                if(indFieldBoundStopPrevTmp(k) > indFieldBoundTmp(k,1) &&...
                        indFieldBoundStartNextTmp(k+1) < indFieldBoundTmp(k,2))
                    validInd(k) = 0;
                else
                    if(~isempty(...
                        find(...
                            avgFRProfile(indNeurons(i),...
                            indFieldBoundStopPrevTmp(k):indFieldBoundTmp(k,1))...
                                < paramFW.lowerBoundFRFieldNeuron...
                                    *meanInstFR(indNeurons(i)),1))...
                            ...% the firing rate at the beginning of the 
                            ...% field should be around 0
                        || ~isempty(...
                            find(...
                                avgFRProfile(indNeurons(i),...
                                    indFieldBoundTmp(k,2):indFieldBoundStartNextTmp(k+1))...
                                    < paramFW.lowerBoundFRFieldNeuron...
                                    *meanInstFR(indNeurons(i)),1)) ...  
                            ...% the firing rate towards the end of the 
                            ...% field should return to 0
                        || (isempty(...
                            find(...
                                avgFRProfile(indNeurons(i),...
                                    indFieldBoundTmp(k,2):indFieldBoundStartNextTmp(k+1))...
                                    < paramFW.lowerBoundFRFieldNeuron...
                                    *meanInstFR(indNeurons(i)),1))...
                                && indFieldBoundTmp(k,2) == size(avgFRProfile,2))) 
                            % the end of the field is aligned with the end 
                            % of the recording trace
                        validInd(k) = 1;
                    end
                end
            end
            
            indFieldBoundTmp = indFieldBoundTmp(validInd == 1,:);
            numFieldsTmp = size(indFieldBoundTmp,1);
        elseif(numFieldsTmp == 0)
            disp(['Neuron ' num2str(indNeurons(i)) ': potential field has a long tail']);
        end
               
        % find the peak within each potential field
        neuronDoubleFieldCount = 0;
        neuronSingFieldCount = 0;
        neuronAllFieldCount = 0;
        for j = 1:numFieldsTmp  
            peakFRTmp = ...
                max(avgFRProfile(indNeurons(i),...
                    indFieldBoundTmp(j,1):indFieldBoundTmp(j,2))); 
                    % the peak firing rate within the field
            indPeakField = ...
                find(avgFRProfile(indNeurons(i),...
                    indFieldBoundTmp(j,1):indFieldBoundTmp(j,2)) == peakFRTmp); 
                    % indices of peaks within the field
            indPeak = indFieldBoundTmp(j,1) + indPeakField(1) - 1; 
                    % the actual index of the peak response
            indSpikeInField = ...
                find(...
                    spDistPerNeuron{indNeurons(i)} >= ...
                        paramGen.spaceSteps(indFieldBoundTmp(j,1)) ...
                    & spDistPerNeuron{indNeurons(i)} <= ...
                        paramGen.spaceSteps(indFieldBoundTmp(j,2)));
            distSpikesInField = ...
                spDistPerNeuron{indNeurons(i)}(indSpikeInField);
            [distSpikesInField,indTmp] = sort(distSpikesInField); 
                % sort the spikes in ascending order of distance
%             phaseSpikeInField = ...
%                 spikeThetaPhaseStruct.spPhaseVsTPerNeuron{indNeurons(i)}(indSpikeInField);
%             phaseSpikeInField = phaseSpikeInField(indTmp); 
                % sort the phases of spikes in ascending order of spike
                % distance
                
            if(indFieldBoundTmp(j,2)-indFieldBoundTmp(j,1) < ...
                    paramFW.FieldMaxLen/paramGen.spaceMergeBin)
                if(indFieldBoundTmp(j,2)-indFieldBoundTmp(j,1) >= ...
                        paramFW.overlapFieldWidthMin/paramGen.spaceMergeBin) 
                        % for wide field, analyze whether it contains two
                        % overlapping fields
                                                                                        
                    [indFieldBoundNew, indPeakPerFieldNew, ...
                        numFieldNew, numSpikesPerFieldNew] ...
                            = AnalyzeField(...
                                avgFRProfile(indNeurons(i),...
                                        indFieldBoundTmp(j,1):indFieldBoundTmp(j,2)),...
                                indPeakField(1),distSpikesInField,...
                                paramFW,paramGen.spaceMergeBin);
                                    % analyzing if there are overlapping fields

                    % record the field information
                    if(numFieldNew ~= 0)
                        % count the number of spikes from each trial and
                        % decide whether a field is good based on how many
                        % trials have spikes within the field
                        numFieldGood = 0;
                        goodField = [];
                        indATmp = indFieldBoundNew(:,1)' + indFieldBoundTmp(j,1) - 1;
                        indBTmp = indFieldBoundNew(:,2)' + indFieldBoundTmp(j,1) - 1;
                        for k = 1:numFieldNew
                            numActTrial = 0;
                            numSpikePerTrialPerField = zeros(1,fieldStruct.numTrials);
                            for m = 1:fieldStruct.numTrials
                                spikePerNeuronPerTrial = ...
                                    spDistPerTrialPerNeuron{indNeurons(i),m};
                                numSpikePerTrialPerField(m) = ...
                                    length(find(spikePerNeuronPerTrial >= ...
                                            paramGen.spaceSteps(indATmp(k))...
                                        & spikePerNeuronPerTrial <= ...
                                            paramGen.spaceSteps(indBTmp(k))));
                                if(numSpikePerTrialPerField(m) > 0)
                                    numActTrial = numActTrial + 1;
                                end
                            end
                            numSpikesPerField = sum(numSpikePerTrialPerField);
%                             if(numSpikesPerField ~= fieldStruct.numSpPerField)
%                                 disp('Total number of spikes in the field does not match.')
%                             end

                            % check whether there are other significant peaks
                            % immediately before the detected field
                            peakOutsideField = 0;
                            if(indFieldBoundTmp(j,1) > 1)
                                indBefFieldStart = max(indFieldBoundTmp(j,1)-paramFW.ReboundCheckRegion,1);
                                if(sum(avgFRProfile(indNeurons(i),indBefFieldStart:indFieldBoundTmp(j,1)) >= ...
                                        paramFW.ReboundHeight*avgFRProfile(indNeurons(i),indPeak(1))) > 0)
                                    peakOutsideField = 1;
                                end
                            end
                            % check whether there are other significant peaks
                            % immediately after the detected field
                            trLen = size(avgFRProfile,2);
                            if(indFieldBoundTmp(j,2) < trLen)
                                indAfterFieldEnd = min(indFieldBoundTmp(j,2)+paramFW.ReboundCheckRegion,trLen);
                                if(sum(avgFRProfile(indNeurons(i),indFieldBoundTmp(j,2):indAfterFieldEnd) >= ...
                                        paramFW.ReboundHeight*avgFRProfile(indNeurons(i),indPeak(1))) > 0)
                                    peakOutsideField = 1;
                                end
                            end
                            if(peakOutsideField == 1)
                                disp(['Neuron ' num2str(indNeurons(i)) ': there are other peaks close to the field']);
                            end
                        
                            if(numActTrial >= ...
                                    floor(paramGen.percActiveTrials*fieldStruct.numTrials) && ...
                                    peakOutsideField == 0)
                                numFieldGood = numFieldGood + 1;
                                fieldStruct.numSpikesPerTrialPerField{numFieldTotal+numFieldGood-1} = ...
                                    numSpikePerTrialPerField;
                                fieldStruct.numTrialsPerField(numFieldTotal+numFieldGood-1) = ...
                                    numActTrial;
                                percSpikeIn = ...
                                    numSpikesPerField...
                                    /length(spDistPerNeuron{indNeurons(i)});
                                fieldStruct.percSpikesPerField(numFieldTotal+numFieldGood-1) = ...
                                    percSpikeIn;
                                goodField = [goodField, k];
                            else
                                disp(['Neuron ' num2str(indNeurons(i)) ': only ' num2str(numActTrial) ...
                                    ' trials have non-zero number of spikes in the field. Threshold: ' ...
                                    num2str(floor(paramGen.percActiveTrials*fieldStruct.numTrials))]);
                            end
                        end
                        
                        if(numFieldGood ~= 0)
                            indPeakTmp = indPeakPerFieldNew(goodField) + indFieldBoundTmp(j,1) - 1;
                            indATmp = indFieldBoundNew(goodField,1)' + indFieldBoundTmp(j,1) - 1;
                            indBTmp = indFieldBoundNew(goodField,2)' + indFieldBoundTmp(j,1) - 1;
                                             
                            fieldStruct.indNeuron = ...
                                [fieldStruct.indNeuron indNeurons(i)*ones(1,numFieldGood)];
                            fieldStruct.FW = ...
                                [fieldStruct.FW (indFieldBoundNew(:,2) - indFieldBoundNew(:,1))'];
                            fieldStruct.indStartField = ...
                                [fieldStruct.indStartField indATmp];    

                            fieldStruct.numSpPerField = ...
                                [fieldStruct.numSpPerField numSpikesPerFieldNew];
                            fieldStruct.indPeakField = ...
                                [fieldStruct.indPeakField indPeakTmp];
                            fieldStruct.peakLocInField = ...
                                [fieldStruct.peakLocInField ...
                                (indPeakPerFieldNew - indFieldBoundNew(:,1)')...
                                ./(indFieldBoundNew(:,2) - indFieldBoundNew(:,1))'];
                            fieldStruct.peakInstFiringRate = ...
                                [fieldStruct.peakInstFiringRate ...
                                avgFRProfile(indNeurons(i),indPeakTmp)];
                            fieldStruct.meanInstFRArr = ...
                                [fieldStruct.meanInstFRArr ...
                                meanInstFR(indNeurons(i))*ones(1,numFieldGood)];
                            if(meanInstFR(indNeurons(i)) ~= 0)
                                fieldStruct.p2MFRField = ...
                                    [fieldStruct.p2MFRField ...
                                    avgFRProfile(indNeurons(i),indPeakTmp)/meanInstFR(indNeurons(i))];
                            else
                                fieldStruct.p2MFRField = ...
                                    [fieldStruct.p2MFRField zeros(1,numFieldGood)];
                            end
                            
                            if(numFieldGood == 1)
                                fieldStruct.fDoubleField = ...
                                    [fieldStruct.fDoubleField 0]; % single field
                            else
                                fieldStruct.fDoubleField = ...
                                    [fieldStruct.fDoubleField ones(1,numFieldGood)]; 
                                    % double field
                            end

                            for k = 1:numFieldGood
                                if(k == 1) 
                                    fieldStruct.distSpikesPerField{numFieldTotal} = ...
                                        distSpikesInField(1:numSpikesPerFieldNew(1));
                                    fieldStruct.phaseSpikesPerField{numFieldTotal} = ...
                                        phaseSpikeInField(1:numSpikesPerFieldNew(1));
                                else
                                    fieldStruct.distSpikesPerField{numFieldTotal+1} = ...
                                        distSpikesInField(numSpikesPerFieldNew(1)+1:end);
                                    fieldStruct.phaseSpikesPerField{numFieldTotal+1} = ...
                                        phaseSpikeInField(numSpikesPerFieldNew(1)+1:end);
                                end
                            end
                        end
                                               
                        numFieldTotal = numFieldTotal + numFieldGood;
                        numFieldNew = numFieldGood;
                    end 

                else % for normal field, it contains only one field
                    numSpikesInField = length(indSpikeInField);
                    if(numSpikesInField >= paramFW.minNumSpikesPerField)                    
                        % count the number of spikes from each trial
                        numActTrial = 0;
                            
                        numSpikePerTrialPerField = zeros(1,fieldStruct.numTrials);
                        for m = 1:fieldStruct.numTrials
                            spikePerNeuronPerTrial = spDistPerTrialPerNeuron{indNeurons(i),m};
                            numSpikePerTrialPerField(m) = ...
                                length(find(spikePerNeuronPerTrial >= ...
                                            paramGen.spaceSteps(indFieldBoundTmp(j,1))...
                                & spikePerNeuronPerTrial <= ...
                                            paramGen.spaceSteps(indFieldBoundTmp(j,2))));
                            if(numSpikePerTrialPerField(m) > 0)
                                numActTrial = numActTrial + 1;
                            end
                        end
                        numSpikesPerField = sum(numSpikePerTrialPerField);
%                         if(numSpikesPerField ~= fieldStruct.numSpPerField)
%                             disp('Total number of spikes in the field does not match.')
%                         end

                        % check whether there are other significant peaks
                        % immediately before the detected field
                        peakOutsideField = 0;
                        if(indFieldBoundTmp(j,1) > 1)
                            indBefFieldStart = max(indFieldBoundTmp(j,1)-paramFW.ReboundCheckRegion,1);
                            if(sum(avgFRProfile(indNeurons(i),indBefFieldStart:indFieldBoundTmp(j,1)) >= ...
                                    paramFW.ReboundHeight*avgFRProfile(indNeurons(i),indPeak(1))) > 0)
                                peakOutsideField = 1;
                            end
                        end
                        % check whether there are other significant peaks
                        % immediately after the detected field
                        trLen = size(avgFRProfile,2);
                        if(indFieldBoundTmp(j,2) < trLen)
                            indAfterFieldEnd = min(indFieldBoundTmp(j,2)+paramFW.ReboundCheckRegion,trLen);
                            if(sum(avgFRProfile(indNeurons(i),indFieldBoundTmp(j,2):indAfterFieldEnd) >= ...
                                    paramFW.ReboundHeight*avgFRProfile(indNeurons(i),indPeak(1))) > 0)
                                peakOutsideField = 1;
                            end
                        end
                        if(peakOutsideField == 1)
                            disp(['Neuron ' num2str(indNeurons(i)) ': there are other peaks close to the field']);
                        end
                        
                        if(numSpikesInField/(distSpikesInField(end)-distSpikesInField(1))...
                            >= paramFW.minSpikeDenPerFieldPerHzTmp*numActTrial && ...
                            peakOutsideField == 0)
                            % check that every trial has enough spikes
                            if(numActTrial >= floor(paramGen.percActiveTrials*fieldStruct.numTrials))
                                percSpikeIn = ...
                                    numSpikesPerField/length(spDistPerNeuron{indNeurons(i)});
                                fieldStruct.percSpikesPerField(numFieldTotal) = percSpikeIn;

                                fieldStruct.numSpikesPerTrialPerField{numFieldTotal} = ...
                                    numSpikePerTrialPerField;
                                fieldStruct.numTrialsPerField(numFieldTotal) = numActTrial;

                                fieldStruct.indNeuron = [fieldStruct.indNeuron indNeurons(i)];
                                fieldStruct.FW = ...
                                    [fieldStruct.FW indFieldBoundTmp(j,2) - indFieldBoundTmp(j,1)];
                                fieldStruct.indStartField = ...
                                    [fieldStruct.indStartField indFieldBoundTmp(j,1)];

                                fieldStruct.numSpPerField = [...
                                    fieldStruct.numSpPerField length(indSpikeInField)];
                                fieldStruct.indPeakField = ...
                                    [fieldStruct.indPeakField indPeak(1)];
                                fieldStruct.peakLocInField = ...
                                    [fieldStruct.peakLocInField ...
                                    (fieldStruct.indPeakField(end) ...
                                    - fieldStruct.indStartField(end))/fieldStruct.FW(end)];
                                fieldStruct.peakInstFiringRate = ...
                                    [fieldStruct.peakInstFiringRate ...
                                    avgFRProfile(indNeurons(i),indPeak(1))];
                                fieldStruct.meanInstFRArr = ...
                                    [fieldStruct.meanInstFRArr meanInstFR(indNeurons(i))];
                                if(meanInstFR(indNeurons(i))~=0)
                                    fieldStruct.p2MFRField = ...
                                        [fieldStruct.p2MFRField ...
                                        avgFRProfile(indNeurons(i),indPeak(1))/meanInstFR(indNeurons(i))];
                                else
                                    fieldStruct.p2MFRField = [fieldStruct.p2MFRField 0];
                                end
                                fieldStruct.fDoubleField = [fieldStruct.fDoubleField 0];

                                fieldStruct.distSpikesPerField{numFieldTotal} = distSpikesInField;
                                fieldStruct.phaseSpikesPerField{numFieldTotal} = phaseSpikeInField;

                                numFieldTotal = numFieldTotal + 1;
                                numFieldNew = 1;
                            else
                                numFieldNew = 0;
                                disp(['Neuron ' num2str(indNeurons(i)) ': only ' ...
                                    num2str(numActTrial) ...
                                    ' trials have non-zero number of spikes in the field. Threshold: ' ...
                                    num2str(floor(paramGen.percActiveTrials*fieldStruct.numTrials))]);
                            end
                        else
                            numFieldNew = 0;
                            disp(['Neuron ' num2str(indNeurons(i)) ' spikes per Hz is ' ...
                                num2str(numSpikesInField/(distSpikesInField(end)-distSpikesInField(1))) ' < ' ...
                                num2str(paramFW.minSpikeDenPerFieldPerHzTmp*numActTrial)]);
                        end
                    else
                        disp(['Neuron ' num2str(indNeurons(i)) ': potential field number of spikes: ' ...
                            num2str(numSpikesInField)]);
                        numFieldNew = 0;
                    end
                end
            else
                disp(['Neuron ' num2str(indNeurons(i)) ': potential field width is: ' ...
                    num2str(indFieldBoundTmp(j,2)-indFieldBoundTmp(j,1))]);
                numFieldNew = 0;
            end

            % classify the neuron according to the property
            % of its field
            if(isempty(find(fieldStruct.neuronDoubleField == indNeurons(i), 1)) &&...
                    numFieldNew == 2) 
                % if there are two subfields, classify the neuron into 
                % neuronDoubleField
                fieldStruct.neuronDoubleField = ...
                    [fieldStruct.neuronDoubleField indNeurons(i)];
                neuronDoubleFieldCount = neuronDoubleFieldCount+1;
            end
            if(isempty(find(fieldStruct.neuronSingleField == indNeurons(i), 1)) && ...
                    numFieldNew == 1) 
                % if there is only one field, classify the neuron into 
                % neuronSingleField
                fieldStruct.neuronSingleField = ...
                    [fieldStruct.neuronSingleField indNeurons(i)];
                neuronSingFieldCount = neuronSingFieldCount+1;
            end   
            if(isempty(find(fieldStruct.neuronFieldAll == indNeurons(i), 1)) && ...
                    numFieldNew ~= 0) 
                % classify the neuron into neuronFieldAll
                fieldStruct.neuronFieldAll = ...
                    [fieldStruct.neuronFieldAll indNeurons(i)];
                neuronAllFieldCount = neuronAllFieldCount + 1;
            end

            numFieldNew = 0;
        end
        
        numFieldTotalCurNeuronEnd = numFieldTotal;
        numFieldCurNeu = numFieldTotalCurNeuronEnd-numFieldTotalCurNeuronStart;
        
        % check the total percentage of spikes within a field
        totInFieldSpikes = 0;
        totFieldWidth = 0;
        if(numFieldCurNeu ~= 0)
            totInFieldSpikes = sum(fieldStruct.percSpikesPerField...
                (numFieldTotalCurNeuronStart:numFieldTotalCurNeuronEnd-1));
            totFieldWidth = sum(fieldStruct.FW...
                (numFieldTotalCurNeuronStart:numFieldTotalCurNeuronEnd-1));
            if(totInFieldSpikes < paramFW.percSpikesInFieldsTmp || ...
                    (totInFieldSpikes < paramFW.percSpikesInFieldsPerUnitTmp*totFieldWidth && ...
                    totInFieldSpikes < paramFW.maxPercSpikesInFieldsTmp))
                disp(['Neuron ' num2str(indNeurons(i)) ...
                    ': percentage of in field spikes is ' ...
                    num2str(totInFieldSpikes)]);
                corrInd = 1:numFieldTotalCurNeuronStart-1;
                fieldStruct.numSpikesPerTrialPerField = ...
                    fieldStruct.numSpikesPerTrialPerField(corrInd);
                fieldStruct.numTrialsPerField = ...
                    fieldStruct.numTrialsPerField(corrInd);
                fieldStruct.percSpikesPerField = ...
                    fieldStruct.percSpikesPerField(corrInd);
                fieldStruct.indNeuron = ...
                    fieldStruct.indNeuron(corrInd);
                fieldStruct.FW = ...
                    fieldStruct.FW(corrInd);
                fieldStruct.indStartField = ...
                    fieldStruct.indStartField(corrInd);    
                fieldStruct.numSpPerField = ...
                    fieldStruct.numSpPerField(corrInd);
                fieldStruct.indPeakField = ...
                    fieldStruct.indPeakField(corrInd);
                fieldStruct.peakLocInField = ...
                    fieldStruct.peakLocInField(corrInd);
                fieldStruct.peakInstFiringRate = ...
                    fieldStruct.peakInstFiringRate(corrInd);
                fieldStruct.meanInstFRArr = ...
                    fieldStruct.meanInstFRArr(corrInd);
                fieldStruct.p2MFRField = ...
                    fieldStruct.p2MFRField(corrInd);
                fieldStruct.fDoubleField = ...
                    fieldStruct.fDoubleField(corrInd);
                fieldStruct.distSpikesPerField = ...
                    fieldStruct.distSpikesPerField(corrInd);
                fieldStruct.phaseSpikesPerField = ...
                    fieldStruct.phaseSpikesPerField(corrInd);
                if(neuronDoubleFieldCount > 0)
                    fieldStruct.neuronDoubleField = ...
                        fieldStruct.neuronDoubleField(1:end-neuronDoubleFieldCount);
                end
                if(neuronSingFieldCount > 0)
                    fieldStruct.neuronSingleField = ...
                        fieldStruct.neuronSingleField(1:end-neuronSingFieldCount);
                end  
                if(neuronAllFieldCount > 0)
                    fieldStruct.neuronFieldAll = ...
                        fieldStruct.neuronFieldAll(1:end-neuronAllFieldCount);
                end  
                numFieldCurNeu = 0;
                numFieldTotalCurNeuronEnd = numFieldTotalCurNeuronStart;
                numFieldTotal = numFieldTotalCurNeuronStart;
            end
        end
        
        if(numFieldCurNeu ~= 0)
            numSpPerField = fieldStruct.numSpPerField(end-numFieldCurNeu+1:end);
            numTrialsL2Sp = zeros(1,numFieldCurNeu);
            numTrialsL1Sp = zeros(1,numFieldCurNeu);
            for k = 1:numFieldCurNeu
                numTrialsL2Sp(k) = ...
                    sum(fieldStruct.numSpikesPerTrialPerField{numFieldTotal-k} > 1);
                    % number of trials which have at least 2 spikes
                numTrialsL1Sp(k) = ...
                    sum(fieldStruct.numSpikesPerTrialPerField{numFieldTotal-k} > 0);
                    % number of trials which have at least 1 spikes
            end
        else
            numSpPerField = 0;
            numTrialsL2Sp = 0;
            numTrialsL1Sp = 0;
        end
        percInFieldThr = min(paramFW.percSpikesInFieldsPerUnitTmp*totFieldWidth,paramFW.maxPercSpikesInFieldsTmp);
        percInFieldThr = max(percInFieldThr,paramFW.percSpikesInFieldsTmp);
        disp(['Number of field of Neuron ' num2str(indNeurons(i)) ': ' ...
                num2str(numFieldCurNeu) ' minNumSpikesPerField = ' ...
                num2str(paramFW.minNumSpikesPerField) ' numSpPerField = ' ...
                num2str(numSpPerField) ' numTrials>2Spikes [' ...
                num2str(numTrialsL2Sp) '] numTrials>1Spikes [' ...
                num2str(numTrialsL1Sp) '] percSpikesInFields [' ...
                num2str(totInFieldSpikes) '] percInFieldThr [' ...
                num2str(percInFieldThr) ']']);
        if(numFieldTotalCurNeuronEnd <= numFieldTotalCurNeuronStart) 
            % if there is no valid field for this neuron, classify it into 
            % neuronNoField
            fieldStruct.neuronNoField = [fieldStruct.neuronNoField indNeurons(i)];
        end
    elseif(isempty(indexHighFR))
        disp(['Maximum inst firing rate of neuron ' num2str(indNeurons(i))...
            ' is ' num2str(meanInstFR(indNeurons(i))) ', peak to mean is ' ...
             num2str(pFRStruct.peakFR(indNeurons(i))...
                /meanInstFR(indNeurons(i))) ', likely did not' ...
            ' reach the inst FR threshold']);      
    end
% 
%     % plot the firing profile
%     if(figureState == 2)
%         fieldInfo = getFieldInfoIndNeuron(indNeurons(i),fieldStruct);
%         if(~isempty(fieldInfo))
%             neuronClass = 'With field';
%         else
%             neuronClass = 'No field';
%         end
%         newplot;
%         clf;
%         plotThetaPhaseIndNeuron(gca,spDistPerNeuron{indNeurons(i)},...
%               spikeThetaPhaseStruct.spPhaseVsTPerNeuron{indNeurons(i)},...
%               minTimeInterval,paramGen.spaceMergeBin);
%         ax = axes('Position',get(gca,'Position'),...
%            'XAxisLocation','top',...
%            'YAxisLocation','right',...
%            'XTickLabel',[],...
%            'Xlim',[0 paramGen.numSamples*paramGen.spaceMergeBin],...
%            'Color','none',...
%            'XColor','k','YColor','k',...
%            'FontSize',14); 
%         figTitle = ['Neuron ' num2str(indNeurons(i)) ' ' neuronClass];
%         plotFRProfileIndNeuron(ax,avgFRProfile(indNeurons(i),:),...
%               meanInstFR(indNeurons(i)),fieldInfo,figTitle,...
%               paramGen.numSamples,paramGen.spaceMergeBin);
% 
%         pause(1);
%     end
end
  
function [indFieldBound, numFields,fringeLen] = ...
    SearchContRegionsDoubleThre(indexHighFR,indexLowFRLT,indexLowFRHT,...
                                FRArr,numSamples,MaxFringeLen)
% find potential fields using double low FR threshold
% indexHighFR: the indices where bumps occur
% indexLowFRHT: the indices where the firing rate is lower than or equal to
% the mean inst firing rate * high threshold
% indexLowFRLT: the indices where the firing rate is lower than or equal to
% the mean inst firing rate * low threshold
% numSamples: the total number of samples in the neuronal response

    indStart = indexHighFR(1); % start index of the current bump
    indFieldBound = []; % the boundary (start and ending indices) of each candidate field
    numFields = 0; % number of fields
    for j = 1:length(indexHighFR)-1
        if(indexHighFR(j+1) - indexHighFR(j) ~= 1 || j == length(indexHighFR)-1)
            indEnd = indexHighFR(j); % end index of the current bump
            indexAHT = indexLowFRHT(indexLowFRHT < indStart);
            indexBHT = indexLowFRHT(indexLowFRHT > indEnd);

            % find the first point before the current bump which
            % touches the low threshold line 
            if(isempty(indexAHT)) 
                % assume this is the case that the bump start from the 
                % beginning of the response
                indexAHT = 1;
            else
                % if it is not the beginning of the response, the take the
                % point closest to the bump
                indexAHT = indexAHT(end);
            end

            % find the first point after the current bump which touches
            % the low threshold line 
            if(isempty(indexBHT)) 
                % assume this is the case that the bump ends at the end of
                % the response
                indexBHT = numSamples;
            else
                % if it is not the end of the response, the take the point
                % closest to the bump
                indexBHT = indexBHT(1);
            end
            
            indexALT = indexLowFRLT(indexLowFRLT < indStart);
            indexBLT = indexLowFRLT(indexLowFRLT > indEnd);

            % find the first point before the current bump which
            % touches the high threshold line 
            if(isempty(indexALT)) 
                % assume this is the case that the bump starts from the 
                % beginning of the response
                indexALT = 1;
            else % if it is not the beginning of the response, the take the 
                 % point closest to the bump
                indexALT = indexALT(end);
            end

            % find the first point after the current bump which touches
            % the high threshold line 
            if(isempty(indexBLT)) 
                % assume this is the case that the bump ends at the end of  
                % the response
                indexBLT = numSamples;
            else % if it is not the end of the response, the take the point
                 % closest to the bump
                indexBLT = indexBLT(1);
            end
            
            % if between the indices calculated with the lower threshold 
            % and the higher threshold, the inst FR rebounds above the high
            % threshold, then set the boundary to be the min inst FR between 
            % indexAHT and the first rebound
            fringeLen = 0; 
            if(indexALT < indexAHT) % before the field bump
                fringeLen = indexAHT - indexALT;
                tmpIndHTtoLT = intersect(indexLowFRHT,indexALT:indexAHT);
                diffIndexLowFR = diff(tmpIndHTtoLT);
                indDiscont = find(diffIndexLowFR ~= 1); % find rebounds
                if(~isempty(indDiscont))
                    indDiscont = tmpIndHTtoLT(indDiscont(end) + 1); 
                    % find the rebound closest to the field bump (the last of all)
                    minInstFR = min(FRArr(indDiscont:indexAHT)); 
                    % find the min inst FR between the end of rebound and indexALT
                    indexALT1 = ...
                        find(FRArr(indDiscont:indexAHT) == minInstFR,1)...
                        + indDiscont(end) - 1; 
                        % set the boundary to the min inst FR location
                    
                end 
            end
            
            if(indexBLT > indexBHT)  % after the field bump
                fringeLen = fringeLen + indexBLT - indexBHT;
                tmpIndHTtoLT = intersect(indexLowFRHT,indexBHT:indexBLT);
                diffIndexLowFR = diff(tmpIndHTtoLT);  % find rebounds
                indDiscont = find(diffIndexLowFR ~= 1);
                if(~isempty(indDiscont))
                    indDiscont = tmpIndHTtoLT(indDiscont(1));  ...
                    % find the rebound closest to the field bump (the first of all)
                    minInstFR = min(FRArr(indexBHT:indDiscont(1))); 
                    % find the min inst FR between indexBLT and the beginning of rebound
                    indexBLT1 = ...
                        find(FRArr(indexBHT:indDiscont(1)) == minInstFR,1,'last')...
                        + indexBHT - 1; 
                        % set the boundary to the min inst FR location
                    
                end
            end
            
            % if the fringe of the field is larger than MaxFringeLen, then
            % it is not a valid field
            if(fringeLen > MaxFringeLen)
                validField = -1;
                disp(['fringe length = ' num2str(fringeLen)]);
            else
                validField = 0;
                % check whether the current peak can be combined with the
                % previous one
                if(numFields ~= 0)
                    if(indFieldBound(numFields,2) > indexALT)
                        validField = -1;
                        indFieldBound(numFields,2) = indexBLT; ...
                            % ending index of the potential field
                    end
                end
            end
            
            % if the field is valid, record its start and ending indices
            if(validField == 0)
                numFields = numFields+1;
                indFieldBound(numFields,1) = indexALT; ...
                    % starting index of the potential field
                indFieldBound(numFields,2) = indexBLT; ...
                    % ending index of the potential field
            end

            indStart = indexHighFR(j+1);
        end
    end
        
function [indFieldBound, indPeakPerField, numField, numSpikesPerField] = ...
        AnalyzeField(frField,indPeakField,spikesInd,paramFW,spaceBin)
% analyze whether the potential field contains overlapping fields, and
% assure that each field contain enough number of spikes
% frField:              the instantaneous firing rate array of the potential field
% indPeakField:         the index of the peak within the field
% spikesInd:            the indices of all the spikes within a field
% paramFW:
%           minNumSpikesPerField: the min number of spikes a field should have
%           minSpikeDenPerField:  the min number of spikes per second within a field
%           threshPeakRatio:      when the higher peak is located in the second half of
%                                 the field, this threshold determines the ratio between the first and
%                                 second peak, larger than which the valley between peaks is used as the
%                                 boundary of two fields, otherwise, the peak of the first field is used as
%                                 the boundary
% minPeakDist:          the min distance between two peaks
% spaceBin:             the spatial bin

    indFieldBound = [];
    indPeakPerField = [];
    numSpikesPerField = [];
    
    numField = 1;
    fieldLen = length(frField); % the length of the field 
    spikesInd = spikesInd - spikesInd(1) + 1; 
        % align the indices to index the first spike as number 1
    halfFieldInd = floor(fieldLen/2);
    numSpikesInField = length(spikesInd); 
        % number of spikes within the potential field
                
    if(numSpikesInField >= paramFW.minNumSpikesPerField && ...
            numSpikesInField/(spikesInd(end)-spikesInd(1))/spaceBin >= ...
                paramFW.minSpikeDenPerFieldPerHzTmp) 
                % if the field is indeed a true field, check whether the 
                % subfields are also true fields
        if(indPeakField <= halfFieldInd) 
                    % if the peak is located at the first half of the field,
                    % then look for the second peak in the second half of the 
                    % field
            indSecondPeak = ...
                find(frField(halfFieldInd+1:end) == max(frField(halfFieldInd+1:end))) ...
                        + halfFieldInd;
                                % find the indices of the second peak 
            if(length(indSecondPeak) > 1) 
                % if there are multiple peaks with the same height, choose
                % the first peak besides the first element of the second half field 
                indTmp = find(indSecondPeak ~= halfFieldInd + 1);
                indSecondPeak = indSecondPeak(indTmp(1));
            end
            indMinInField = ...
                find(...
                frField(indPeakField:indSecondPeak) == ...
                    min(frField(indPeakField:indSecondPeak)))...
                    + indPeakField - 1;
                                % find the min between the first and second peaks
            if(length(indMinInField) ~= 1) 
                % if there are multiple minimums, choose the first one which
                % is not the second peak
                indTmp = find(indMinInField ~= indSecondPeak);
                indMinInField = indMinInField(indTmp(1));
            end
            if(indMinInField == indSecondPeak) 
                % if the second peak is also the minimum, then look for the 
                % true second peak in the second half of the field
                indMinInFieldTmp = [];
                indSecondPeakTmp = [];
                numPeak = 0;
                for i = indSecondPeak+1:fieldLen 
                    % search one by one, until the maximum is not the first ..
                    % element of the shrinking field 
                    maxTmp = max(frField(i:end));
                    if(maxTmp > frField(i))
                        indSecondPeak = find(frField(i:end) == maxTmp) + i - 1; 
                            % find the index of the second peak
                        indSecondPeak = indSecondPeak(1);
                        indMinInField = ...
                            find(frField(indPeakField:indSecondPeak) == ...
                                min(frField(indPeakField:indSecondPeak))) ...
                                + indPeakField - 1;
                        indMinInField = indMinInField(1); 
                            % redefine the index of the min
                        break;
                    end
                end
                
                if(i == fieldLen+1)
                    indMinInField = fieldLen;
                end
            end

            % separate the fields according to indMinInField when the first
            % peak is the higher peak and the distance between the two peak
            % is larger than minPeakDist
            if(indMinInField ~= fieldLen && ...
                    indSecondPeak-indPeakField > paramFW.minPeakDist)
                indFieldBound = [1,indMinInField;indMinInField+1,fieldLen];
                indPeakPerField = [indPeakField, indSecondPeak];
                numField = 2;
            else
                indFieldBound = [1 fieldLen];
                indPeakPerField = indPeakField;
                numField = 1;
            end
            
        else % if the peak is located at the second half of the field, then
             % look for the second peak in the first half of the field
            indSecondPeak = ...
                find(frField(1:halfFieldInd) == max(frField(1:halfFieldInd)));
                                % find the indices of the second peak 
            if(length(indSecondPeak) > 1) 
                    % if there are multiple peaks with the same height,
                    % choose the first peak besides the last element of the
                    % first half field 
                indTmp = find(indSecondPeak ~= halfFieldInd);
                indSecondPeak = indSecondPeak(indTmp(1));
            end
            indMinInField = find(frField(indSecondPeak:indPeakField) == ...
                                min(frField(indSecondPeak:indPeakField))) ...
                                + indSecondPeak - 1;
                                % find the min between the first and second peaks
            if(length(indMinInField) ~= 1) 
                % if there are multiple minimums, choose the first one which
                % is not the second peak
                indTmp = find(indMinInField ~= indSecondPeak);
                indMinInField = indMinInField(indTmp(1));
            end
            if(indMinInField == indSecondPeak) 
                % if the second peak is also the minimum, then look for the
                % true second peak in the second half of the field
                indMinInFieldTmp = [];
                indSecondPeakTmp = [];
                numPeak = 0;
                for i = indSecondPeak-1:-1:1 
                    % search one by one, until the maximum is not the last
                    % element of the shrinking field 
                    maxTmp = max(frField(1:i));
                    if(maxTmp > frField(i))
                        indSecondPeak = find(frField(1:i) == maxTmp); 
                        % find the index of the second peak
                        indSecondPeak = indSecondPeak(end);
                        indMinInField = ...
                            find(frField(indSecondPeak:indPeakField) == ...
                                min(frField(indSecondPeak:indPeakField))) ...
                                + indSecondPeak - 1;
                        indMinInField = indMinInField(1); 
                            % redefine the index of the min
                        break;
                    end
                end
                
                if(i == 0)
                    indMinInField = 1;
                end
            end

            % separate the fields according to indSecondPeak when the second
            % peak is the higher peak and the distance between the two peak
            % is larger than minPeakDist
            if(indMinInField ~= 1 && ...
                    indPeakField - indSecondPeak > paramFW.minPeakDist)
                if(frField(indSecondPeak)/frField(indPeakField) < ...
                        paramFW.threshPeakRatio)
                    indFieldBound = [1,indSecondPeak;indSecondPeak+1,fieldLen];
                else
                    indFieldBound = [1,indMinInField;indMinInField+1,fieldLen];
                end
                indPeakPerField = [indSecondPeak, indPeakField];
                numField = 2;
            else
                indFieldBound = [1 fieldLen];
                indPeakPerField = indPeakField;
                numField = 1;
            end
        end
 
        numSpikesPerField = numSpikesInField; 
        % in the default case, the field contains only one field (or the 
        % two subfields might combined into one field)
        if(numField > 1) % if there are overlapping fields
            fieldQualified = 0;
            numSpikesPerFieldTmp = [];
            for k = 1:numField
                indSpikesInField = ...
                    find(spikesInd >= indFieldBound(k,1) & ...
                         spikesInd <= indFieldBound(k,2));
                numSpikesPerFieldTmp(1,k) = length(indSpikesInField); 
                % record the number of spikes within a subfield
                if(numSpikesPerFieldTmp(1,k) < paramFW.minNumSpikesPerField  &&...
                        numSpikesInField/(spikesInd(end)-spikesInd(1))/spaceBin >=...
                        paramFW.minSpikeDenPerFieldPerHzTmp) 
                    % if the number of spikes in the subfield is lower than
                    % the threshold, combine the two subfields
                    fieldQualified = -1;
                    indFieldBound = [1 fieldLen];
                    numField = 1;
                    break;
                end
            end
            if(fieldQualified ~= -1) 
                % the subfields are true subfields, the record the number 
                % of spikes within each subfield
                numSpikesPerField = numSpikesPerFieldTmp;
            else
                indFieldBound = [1 fieldLen];
                indPeakPerField = indPeakField;
                numField = 1;
            end
        end
    else % if the number of spikes is lower than the threshold, the field 
         % is not a true field
        indFieldBound = [];
        numField = 0;
        numSpikesPerField = 0;
    end
