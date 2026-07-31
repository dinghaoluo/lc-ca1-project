function fieldStruct = FieldSpInfoCorr2P(pFRStruct,spatialInfo,meanCorr,paramF)
% Firing field detection based on correlation distance and spatial information
    
    if(isempty(pFRStruct))
        fieldStruct = [];
        return;
    end

    numNeurons = length(pFRStruct.peakFR);
    if(numNeurons == 0)
        fieldStruct = [];
        return;
    end
    
    GlobalConst2P;
    
    std = floor(spaceBin/spaceMergeBin/2)*2;
    gaussFilt = gaussFilter2P(12*std, std);
    lenGaussKernel = length(gaussFilt);
    normFactor = sum(gaussFilt);
    gaussFilt = gaussFilt./normFactor;
    
    fieldStruct = struct(...
                'indNeuron',[],... % neuron index
                'numTrials', length(pFRStruct.indLapList),... % number of trials
                'FW',[],... % field width (samples)
                'indStartField',[],... % the start index of the field
                'indEndField',[],... % the end index of the field
                'indPeakField',[],... % index of peak firing rate within the field
                'peakInstFiringRate',[],... % peak instantaneous firing rate
                'indPeakFieldNorm',[],... % index of peak firing rate within the field (normalized FR)
                'peakInstFiringRateNorm',[],... % peak instantaneous firing rate (normalized FR)
                'meanInstFRArr',[],... % mean inst firing rate array
                'spatialInfo',[],... % spatial information
                'meanCorrNZ',[],... % mean correlation non-zero trials
                'numActiveFieldTr',[],... % number of trials that active within the field
                'numActiveTrials',[]); % number of active trials
            
    fieldStruct.indNeuron = find(pFRStruct.meanInstFR > paramF.minInstFR & ...
        ((spatialInfo >= paramF.minSpInfo & meanCorr.meanGoodNZ >= paramF.minCorr) | ...
        (spatialInfo >= paramF.minHighSpInfo & meanCorr.meanGoodNZ >= paramF.minHighSpInfoCorr) | ...
        (spatialInfo >= paramF.minCorrHighSpInfo & meanCorr.meanGoodNZ >= paramF.minCorrHigh)) & ...
        fieldStruct.numTrials >= paramF.minNumTr & ...
        meanCorr.nGoodNonZeroTr >= paramF.percNumActiveTr*fieldStruct.numTrials);
    fieldStruct.indPeakField = pFRStruct.peakFRInd(fieldStruct.indNeuron);
    fieldStruct.peakInstFiringRate = pFRStruct.peakFR(fieldStruct.indNeuron);
    fieldStruct.meanInstFRArr = pFRStruct.meanInstFR(fieldStruct.indNeuron);
    fieldStruct.spatialInfo = spatialInfo(fieldStruct.indNeuron);
    fieldStruct.meanCorrNZ = meanCorr.meanGoodNZ(fieldStruct.indNeuron);
    fieldStruct.numActiveTrials = meanCorr.nGoodNonZeroTr(fieldStruct.indNeuron);
    
    indBadField = zeros(1,length(fieldStruct.indNeuron));
    filteredSpike = zeros(length(fieldStruct.indNeuron),size(pFRStruct.avgFRProfileNorm,2));
    for i = 1:length(fieldStruct.indNeuron)
        [peakFRTmp,indPeakTmp] = max(pFRStruct.avgFRProfileNorm(fieldStruct.indNeuron(i),:)); 
        fieldStruct.indPeakFieldNorm(i) = indPeakTmp;
        fieldStruct.peakInstFiringRateNorm(i) = peakFRTmp;
        
        filteredSpikeTmp = pFRStruct.avgFRProfileNorm(fieldStruct.indNeuron(i),:);
        numSamples = size(pFRStruct.avgFRProfileNorm,2);
        filteredSpikeTmp = [filteredSpikeTmp(numSamples-lenGaussKernel+1:numSamples)...
                        filteredSpikeTmp filteredSpikeTmp(1:lenGaussKernel)];
        filteredSpikeTmp = conv(filteredSpikeTmp,...
                gaussFilt);
        filteredSpikeTmp = ...
                filteredSpikeTmp(floor(lenGaussKernel/2)+lenGaussKernel+1:...
                    (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1)); 
        filteredSpike(i,:) = filteredSpikeTmp;
        [peakFRTmp,indPeakTmp] = max(filteredSpikeTmp);
        fieldStruct.indPeakFieldNorm(i) = indPeakTmp;
        fieldStruct.peakInstFiringRateNorm(i) = peakFRTmp;
        
        indStart = ...
            find(filteredSpikeTmp(1:indPeakTmp(1))...
            <= paramF.lowThreFieldMeanInstFR*peakFRTmp,1,'last');
        if(isempty(indStart))
            indStart = 1;
        end
        fieldStruct.indStartField(i) = indStart;
        indEnd = ...
            find(filteredSpikeTmp(indPeakTmp+1:end)...
            <= paramF.lowThreFieldMeanInstFR*peakFRTmp(1),1,'first');
        
        if(isempty(indEnd))
            fieldStruct.indEndField(i) = size(pFRStruct.avgFRProfileNorm,2);
            fieldStruct.FW(i) = size(pFRStruct.avgFRProfileNorm,2)-indStart+1;
        else
            fieldStruct.indEndField(i) = indPeakTmp+indEnd; % removed +1 on 5/23/2023 
            fieldStruct.FW(i) = indPeakTmp+indEnd+1-indStart+1;
        end
        
        indNeuronTmp = [375];
        if(sum(indNeuronTmp == fieldStruct.indNeuron(i)) > 0)
            a = 1;
        end
        
        %% number of trials that having spikes within the field
        nonZeroTrField = sum(pFRStruct.filteredSpArrAll{fieldStruct.indNeuron(i)}(:,indStart:fieldStruct.indEndField(i)),2);
        fieldStruct.numActiveFieldTr(i) = sum(nonZeroTrField > 0);
        if(fieldStruct.numActiveFieldTr(i) < paramF.percNumActiveTr*fieldStruct.numTrials)
            indBadField(i) = 1;
        end

            %% check rebound before the start of field
        if(indStart > 1) 
            indBefFieldStart = max(indStart-paramF.reboundCheckRegion,1);
            if((mean(filteredSpikeTmp(indBefFieldStart:indStart)) >= ...
                    paramF.maxReboundMean*fieldStruct.peakInstFiringRateNorm(i) ||...
                    max(filteredSpikeTmp(indBefFieldStart:indStart)) >= ...
                    paramF.reboundHeight*fieldStruct.peakInstFiringRateNorm(i)))
                indBadField(i) = 1;
            end
        end
        
        %% check rebound after the end of field
        if(fieldStruct.indEndField(i) < paramF.intervalTSpInfo)
            indAfterFieldEnd = min(fieldStruct.indEndField(i)+paramF.reboundCheckRegion,paramF.intervalTSpInfo);
            if((mean(filteredSpikeTmp(fieldStruct.indEndField(i):indAfterFieldEnd)) >= ...
                    paramF.maxReboundMean*fieldStruct.peakInstFiringRateNorm(i) ||...
                    max(filteredSpikeTmp(fieldStruct.indEndField(i):indAfterFieldEnd)) >= ...
                    paramF.reboundHeight*fieldStruct.peakInstFiringRateNorm(i)))
                indBadField(i) = 1;
            end
        end      
                
        figure(1)
        plot(pFRStruct.avgFRProfileNorm(fieldStruct.indNeuron(i),:))
        hold on
        plot(fieldStruct.indStartField(i),...
            pFRStruct.avgFRProfileNorm(fieldStruct.indNeuron(i),fieldStruct.indStartField(i)),'ro');
        plot(fieldStruct.indEndField(i),...
            pFRStruct.avgFRProfileNorm(fieldStruct.indNeuron(i),fieldStruct.indEndField(i)),'ro');
        hold off
    end
    
    % remove neurons with a wide field
    ind = find(fieldStruct.FW < min(paramF.maxFieldWidth,paramF.maxFieldWidth1) & indBadField == 0);
    if(length(ind) < length(fieldStruct.indNeuron))
        fieldStruct.indNeuron = fieldStruct.indNeuron(ind);
        fieldStruct.FW = fieldStruct.FW(ind);
        fieldStruct.indStartField = fieldStruct.indStartField(ind);
        fieldStruct.indEndField = fieldStruct.indEndField(ind);
        fieldStruct.indPeakField = fieldStruct.indPeakField(ind);
        fieldStruct.peakInstFiringRate = fieldStruct.peakInstFiringRate(ind);
        fieldStruct.meanInstFRArr = fieldStruct.meanInstFRArr(ind);
        fieldStruct.spatialInfo = fieldStruct.spatialInfo(ind);
        fieldStruct.meanCorrNZ = fieldStruct.meanCorrNZ(ind);
        fieldStruct.numActiveTrials = fieldStruct.numActiveTrials(ind);
        fieldStruct.indPeakFieldNorm = fieldStruct.indPeakFieldNorm(ind);
        fieldStruct.peakInstFiringRateNorm = fieldStruct.peakInstFiringRateNorm(ind);
        fieldStruct.numActiveFieldTr = fieldStruct.numActiveFieldTr(ind);
        fieldStruct.filteredSpike = filteredSpike(ind,:);
    end
end
