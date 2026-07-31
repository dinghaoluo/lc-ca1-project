function [smFR, smBinP, normSmoothedFR, binnedSpikes, binnedTime, binP] =...
            binCoordiate(distSp,distTr,param)
    %% bin the coordinates
    
    trBin = round(distTr/param.divDist)+1;
    binnedTime = accumarray(trBin,1)/param.sampleFq;
    bins = unique(trBin);
    spBin = round(distSp/param.divDist)+1;
    binnedSpikes = accumarray(spBin,1);
    maxSpikeBin = max(spBin);
    M = max(bins);
    if(maxSpikeBin < M)
        binnedSpikes = [binnedSpikes; zeros(M-maxSpikeBin,1)];
    end
    binnedFR = binnedSpikes ./ binnedTime;
    binP = binnedTime ./ (sum(binnedTime));

    binnedFR(isnan(binnedFR))=0; 
    binnedFR(isinf(binnedFR))=0;
    smFR = SmoothPix(binnedFR,param.smooth/param.divDist); % changed from param.smooth to param.smooth/param.divDist on 11/7/2021
    smBinP = SmoothPix(binP,param.smooth/param.divDist); % changed from param.smooth to param.smooth/param.divDist on 11/7/2021 by Yingxue
    normSmoothedFR = smFR ./ nanmax(nanmax(smFR));
    
end
