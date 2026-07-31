function [sparsity, SNR] = sparsityInfo(smBinP, smFR, meanFR)
    % calculate the sparsity
    
    sparsity = sum(smBinP .* (smFR.^2)) ./ meanFR^2;

    % selectivity: bin with max mean FR / overall mean FR = SNR
    SNR = max(smFR(:)) / meanFR;

end
