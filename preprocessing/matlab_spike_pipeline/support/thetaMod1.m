function [P,f,thetaModul] = thetaMod1(ACG,T)
    % calculate theta modulation
    % theta modulation is defined as the power between 6 to 10 Hz divided
    % by the power between 1 and 50 Hz
    % Deshmukh, J. Neurophysiology, 2009
    
    shortAutoCorr = ACG(T>-10000 & T<10000);
    nfft = 2^nextpow2(length(shortAutoCorr));
    [P,f] = pmtm(shortAutoCorr, 4, nfft, 1000);
    thetaModul = sum(P(f>6&f<10)) / sum(P(f>=1&f<=50));   
    
end
