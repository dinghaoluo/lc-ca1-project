function arrTimeNew = resamp(arrTime,sampleFreq)
% resample all the time stamps in an array according to the sampling
% frequency
% arrTime:      an array of time stamps to be resamples (in millisecond)
% sampleFreq:   sampling frequency

    lenArr = length(arrTime);
    arrTimeNew = zeros(lenArr,1);
    
    for i = 1:lenArr
        a = arrTime(i);
        b = round(a/1000*sampleFreq);
        arrTimeNew(i) = b;
    end
