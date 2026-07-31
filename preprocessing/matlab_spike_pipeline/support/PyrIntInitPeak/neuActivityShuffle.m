function ratioAftBef = neuActivityShuffle(filteredSpikeArray,indAft,indBef,numShuffle)
%     numShuffle = 1000;
%     p = [99.9,99];
    
    rowArray = size(filteredSpikeArray,1);
    colArray = size(filteredSpikeArray,2);
    shufMeanArray = zeros(numShuffle,colArray); 
    parfor i = 1:numShuffle
        shufSpikeArray = zeros(rowArray,colArray);
        randShift = randi(floor(colArray/2),1,rowArray)-floor(colArray/2);
        for j = 1:rowArray
            shiftTmp = circshift(filteredSpikeArray(j,:)',randShift(j));
            shufSpikeArray(j,:) = shiftTmp';
        end
        shufMeanArray(i,:) = mean(shufSpikeArray);
    end
%     avgShuf = mean(shufMeanArray);
%     sigShuf = prctile(shufMeanArray,p);
    
    meanAft = mean(shufMeanArray(:,indAft),2);
    meanBef = mean(shufMeanArray(:,indBef),2);
    ratioAftBef = meanAft./meanBef;
end