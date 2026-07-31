function [smSpInfo, addSmFR] = getSpInfo_addaptBin1(binnedSpikes, binnedTime, binP, meanFR)
% function [smSpInfo, addSmFR] = getSpInfo_addaptBin1(binnedSpikes, binnedTime, meanFR)
% adaptive smoothing based on Skaggs et al., 1992

addSmFR = zeros(size(binnedSpikes, 1),size(binnedSpikes, 2));
smR = zeros(size(binnedSpikes, 1),size(binnedSpikes, 2));
alfa = 10e2;
for nX = 1 : size(binnedSpikes, 1)
    for nY = 1 : size(binnedSpikes, 2)
        Nsp = binnedSpikes(nX, nY);
        Nocc = binnedTime(nX, nY);
        r=1;
        minX = nX;
        maxX = nX;
        minY = nY;
        maxY = nY;
        while Nsp <= alfa / (Nocc^2 * r^2)
            % expand r around the central bin until.....
            minX = nX-r;
            maxX = nX+r;
            minY = nY-r;
            maxY = nY+r;
            if minX <= 0, minX=1; end
            if maxX > size(binnedSpikes, 1), maxX=size(binnedSpikes, 1); end
            if minY <= 0, minY=1; end
            if maxY > size(binnedSpikes, 2), maxY=size(binnedSpikes, 2); end
            
            Nsp = sum(sum(binnedSpikes(minX:maxX, minY:maxY)));
            Nocc = sum(sum(binnedTime(minX:maxX, minY:maxY)));
            
            if r > size(binnedSpikes, 1) / 2
                addSmFR(nX,nY) = NaN;
                smR(nX,nY) = NaN;
                break;
            end
            
            r = r + 1;
            
        end
        smR(nX,nY) = r;
%         addSmFR(nX,nY) =  20 * sum(sum(binnedSpikes(minX:maxX, minY:maxY))) / sum(sum(binnedTime(minX:maxX, minY:maxY)));
        addSmFR(nX,nY) =  1 * sum(sum(binnedSpikes(minX:maxX, minY:maxY))) / sum(sum(binnedTime(minX:maxX, minY:maxY)));
    end
end

smSpInfoMap = (log2(addSmFR./meanFR)) .* addSmFR./meanFR .* binP;
smSpInfoMap(smSpInfoMap<=0) = NaN;
%smSpInfoMap(isnan(smSpInfoMap)) = 0;
smSpInfo = nansum(smSpInfoMap(:));
if isinf(smSpInfo)
    smSpInfo = nan;
end
