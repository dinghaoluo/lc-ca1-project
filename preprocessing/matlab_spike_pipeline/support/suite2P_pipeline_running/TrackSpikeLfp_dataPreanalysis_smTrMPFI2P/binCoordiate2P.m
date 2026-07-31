function [smFR, smBinP, normSmoothedFR, binnedSpikes, binnedTime, binP] =...
            binCoordiate2P(distSp,distTr,param)
    %% bin the coordinates
    
    trBin = round(distTr/param.divDist)+1;
    binnedTime = accumarray(trBin,1)/param.sampleFq;
    bins = unique(trBin);
    binnedSpikes = accumarray(trBin,distSp);
%     maxSpikeBin = max(spBin);
%     M = max(bins);
%     if(maxSpikeBin < M)
%         binnedSpikes = [binnedSpikes; zeros(M-maxSpikeBin,1)];
%     end
    binnedFR = binnedSpikes ./ binnedTime;
    binP = binnedTime ./ (sum(binnedTime));

    binnedFR(isnan(binnedFR))=0; 
    binnedFR(isinf(binnedFR))=0;
    smFR = SmoothPix2P(binnedFR,param.smooth/param.divDist);
    smBinP = SmoothPix2P(binP,param.smooth/param.divDist);
    normSmoothedFR = smFR ./ nanmax(nanmax(smFR));
    
end
