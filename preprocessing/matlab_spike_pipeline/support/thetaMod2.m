function [trough,peak,thetaModInd] = thetaMod2(ACG,T)
    % calculate theta modulation index
    % theta modulation index is defined as  the difference between the 
    % theta modulation trough (defined as mean of autocorrelogram bins, 
    % 50-70 msec) and the theta modulation peak (mean of autocorrelogram 
    % bins, 100-140 msec) over their sum.
    % Cacucci, J.Neurosci., 2004
    
    absT = abs(T);
    
    troughInd = absT >= 50 & absT < 80; % changed from 70 to 80 on 4/7/2020
    troughACG = ACG(troughInd);
    trough = mean(troughACG);
    
    peakInd = absT >= 100 & absT <= 160; % changed from 140 to 160 on 4/7/2020
    peakACG = ACG(peakInd);
    peak = mean(peakACG);
    
    thetaModInd = (peak - trough)/(peak + trough);
    
end
