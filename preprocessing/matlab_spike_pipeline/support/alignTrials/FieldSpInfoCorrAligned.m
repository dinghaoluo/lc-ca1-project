function fieldStruct = FieldSpInfoCorrAligned(pFRStruct,spatialInfo,meanCorr,paramF)
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
    
    fieldStruct = struct(...
                'indNeuron',[],... % neuron index
                'numTrials', length(pFRStruct.indLapList),... % number of trials
                'FW',[],... % field width (samples)
                'indStartField',[],... % the start index of the field
                'indEndField',[],... % the end index of the field
                'indPeakField',[],... % index of peak firing rate within the field
                'peakInstFiringRate',[],... % peak instantaneous firing rate
                'meanInstFRArr',[],... % mean inst firing rate array
                'spatialInfo',[],... % spatial information
                'meanCorrNZ',[],... % mean correlation non-zero trials
                'numActiveTrials',[]); % perc. of active trials
            
    fieldStruct.indNeuron = find(pFRStruct.meanInstFR > paramF.minInstFR & ...
        (spatialInfo >= paramF.minSpInfo & meanCorr.meanGoodNZ >= paramF.minCorr | ...
        spatialInfo >= paramF.minHighSpInfo & meanCorr.meanGoodNZ >= paramF.minCorrHighSpInfo) & ...
        fieldStruct.numTrials >= paramF.minNumTr & ...
        meanCorr.nGoodNonZeroTr >= paramF.percNumActiveTr*fieldStruct.numTrials);
    fieldStruct.indPeakField = pFRStruct.peakFRInd(fieldStruct.indNeuron);
    fieldStruct.peakInstFiringRate = pFRStruct.peakFR(fieldStruct.indNeuron);
    fieldStruct.meanInstFRArr = pFRStruct.meanInstFR(fieldStruct.indNeuron);
    fieldStruct.spatialInfo = spatialInfo(fieldStruct.indNeuron);
    fieldStruct.meanCorrNZ = meanCorr.meanGoodNZ(fieldStruct.indNeuron);
    fieldStruct.numActiveTrials = meanCorr.nGoodNonZeroTr(fieldStruct.indNeuron);
    
    indBadField = zeros(1,length(fieldStruct.indNeuron));
    for i = 1:length(fieldStruct.indNeuron)
        indStart = ...
            find(pFRStruct.avgFRProfile(fieldStruct.indNeuron(i),1:fieldStruct.indPeakField(i))...
            <= paramF.lowThreFieldMeanInstFR*fieldStruct.peakInstFiringRate(i),1,'last');
        if(isempty(indStart))
            indStart = 1;
        end
        fieldStruct.indStartField(i) = indStart;
        indEnd = ...
            find(pFRStruct.avgFRProfile(fieldStruct.indNeuron(i),fieldStruct.indPeakField(i)+1:end)...
            <= paramF.lowThreFieldMeanInstFR*fieldStruct.peakInstFiringRate(i),1,'first');
        
        if(isempty(indEnd))
            fieldStruct.indEndField(i) = size(pFRStruct.avgFRProfile,2);
            fieldStruct.FW(i) = size(pFRStruct.avgFRProfile,2)-indStart+1;
        else
            fieldStruct.indEndField(i) = fieldStruct.indPeakField(i)+indEnd+1;
            fieldStruct.FW(i) = fieldStruct.indPeakField(i)+indEnd+1-indStart+1;
        end
        
        indNeuronTmp = [66 75 81 94 104];
        if(sum(indNeuronTmp == fieldStruct.indNeuron(i)) > 0)
            a = 1;
        end
        
        %% check rebound before the start of field
        if(indStart > 1) 
            indBefFieldStart = max(indStart-paramF.reboundCheckRegion,1);
            if((mean(pFRStruct.avgFRProfile(fieldStruct.indNeuron(i),indBefFieldStart:indStart)) >= ...
                    paramF.maxReboundMean*fieldStruct.peakInstFiringRate(i) ||...
                    max(pFRStruct.avgFRProfile(fieldStruct.indNeuron(i),indBefFieldStart:indStart)) >= ...
                    paramF.reboundHeight*fieldStruct.peakInstFiringRate(i)) && ...
                    (indStart-indBefFieldStart+1) >= paramF.reboundCheckRegion*0.25)
                indBadField(i) = 1;
            end
        end
        
        %% check rebound after the end of field
        if(fieldStruct.indEndField(i) < paramF.intervalTSpInfo)
            indAfterFieldEnd = min(fieldStruct.indEndField(i)+paramF.reboundCheckRegion,paramF.intervalTSpInfo);
            if((mean(pFRStruct.avgFRProfile(fieldStruct.indNeuron(i),fieldStruct.indEndField(i):indAfterFieldEnd)) >= ...
                    paramF.maxReboundMean*fieldStruct.peakInstFiringRate(i) ||...
                    max(pFRStruct.avgFRProfile(fieldStruct.indNeuron(i),fieldStruct.indEndField(i):indAfterFieldEnd)) >= ...
                    paramF.reboundHeight*fieldStruct.peakInstFiringRate(i)) && ...
                    (indAfterFieldEnd-fieldStruct.indEndField(i)+1) >= paramF.reboundCheckRegion*0.25)
                indBadField(i) = 1;
            end
        end       
%         figure(1)
%         plot(pFRStruct.avgFRProfile(fieldStruct.indNeuron(i),:))
%         hold on
%         plot(fieldStruct.indStartField(i),...
%             pFRStruct.avgFRProfile(fieldStruct.indNeuron(i),fieldStruct.indStartField(i)),'ro');
%         plot(fieldStruct.indEndField(i),...
%             pFRStruct.avgFRProfile(fieldStruct.indNeuron(i),fieldStruct.indEndField(i)),'ro');
%         hold off
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
    end
end
