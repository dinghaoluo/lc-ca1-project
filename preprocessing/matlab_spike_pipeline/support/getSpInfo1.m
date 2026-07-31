function [smSpInfo, varSmFR] = getSpInfo1(smFR, smBinP, meanFR)
    % calculate spatial information
    
    smSpInfoMap = (log2(smFR./meanFR)) .* smFR./meanFR .* smBinP;
    smSpInfoMap(isnan(smSpInfoMap)) = 0;
    smSpInfo = sum(smSpInfoMap(:));
    if isinf(smSpInfo)
        smSpInfo = nan;
    end
%         smSpInfo(divXY, smooth) = sum(sum(smSpInfoMap));
    varSmFR = var(smFR(:)-meanFR);
end
