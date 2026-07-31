function [trough,peak,troughT,peakT,thetaModInd] = thetaMod3(ACG,T)
    % calculate theta modulation index
    % theta modulation index is defined as  the difference between the 
    % theta modulation trough (defined as mean of autocorrelogram bins, 
    % 50-70 msec) and the theta modulation peak (mean of autocorrelogram 
    % bins, 100-140 msec) over their sum.
    % Yingxue
    
    absT = abs(T);
    
    troughInd = find(absT >= 30 & absT < 100); 
    troughACG = ACG(troughInd);
    [trough,troughIdx] = min(troughACG);
    troughT = T(troughInd(troughIdx));
    
    peakInd = find(absT >= 100 & absT <= 200); 
    peakACG = ACG(peakInd);
    [peak,peakIdx] = max(peakACG);
    peakT = T(peakInd(peakIdx));
    
    thetaModInd = (peak - trough)/(peak + trough);
    
end
