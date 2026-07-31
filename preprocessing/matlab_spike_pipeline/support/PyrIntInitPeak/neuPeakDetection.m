function [isPeakNeu,ratioAftBefShuf] = neuPeakDetection(filteredSpikeArray,indAft,indBef,p,numShuffle)
    isPeakNeu = zeros(1,length(p));
    ratioAftBefShuf = neuActivityShuffle(filteredSpikeArray,indAft,indBef,numShuffle);
    
    avgProfile = mean(filteredSpikeArray);
    ratioAftBef = mean(avgProfile(indAft))/mean(avgProfile(indBef));
    
    isPeakNeu = zeros(1,length(p));
    for i = 1:length(p)
        sigShuf = prctile(ratioAftBefShuf,[p(i) 100-p(i)]);

        if(ratioAftBef > sigShuf(1))
            isPeakNeu(i) = 1;
        elseif(ratioAftBef < sigShuf(2))
            isPeakNeu(i) = -1;
        end
    end
   
end
