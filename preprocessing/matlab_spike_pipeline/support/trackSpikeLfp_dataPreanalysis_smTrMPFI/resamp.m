function arrTimeNew = resamp(arrTime,sampleFreq)
% resample all the time stamps in an array according to the sampling
% frequency
% arrTime:      an array of time stamps to be resamples (in millisecond)
% sampleFreq:   sampling frequency

    lenArr = length(arrTime);
    arrTimeNew = zeros(lenArr,1);
    
    for i = 1:lenArr
        arrTimeNew(i) = round(arrTime(i)/1000*sampleFreq);
    end
