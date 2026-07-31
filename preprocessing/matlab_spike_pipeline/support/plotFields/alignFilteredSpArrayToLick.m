function filteredSpikeArrayByLick = alignFilteredSpArrayToLick(indLapOrderedByLick,firstLickPos,filteredSpikeArray,alignDist)
% order the filtered spike array by the time of the first lick

    numNeurons = size(filteredSpikeArray{1},1);
    
    numSess = length(indLapOrderedByLick);
    trackLen = size(filteredSpikeArray{1},2);
    for neu = 1:numNeurons
        filteredSpikeArrayByLick{neu} = ...
        zeros(length(filteredSpikeArray),trackLen);
        for i = 1:numSess
            for j = 1:length(indLapOrderedByLick{i})
                indLap = indLapOrderedByLick{i}(j);
                tempSpikeArr = zeros(1,size(filteredSpikeArray{1},2));
                lickPos = floor(firstLickPos(indLap));
                if(lickPos < alignDist)
                    tempSpikeArr(alignDist-lickPos+1:end) = ...
                        filteredSpikeArray{indLap}(neu,1:trackLen-alignDist+lickPos);    
                elseif(firstLickPos(indLap) > alignDist)
                    tempSpikeArr(1:trackLen-lickPos+alignDist) = ...
                        filteredSpikeArray{indLap}(neu,lickPos-alignDist+1:end);  
                else
                    tempSpikeArr = filteredSpikeArray{indLap}(neu,:);
                end
                filteredSpikeArrayByLick{neu}(indLap,:) = tempSpikeArr;
            end
        end
    end
end
