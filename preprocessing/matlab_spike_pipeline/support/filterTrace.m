function filtData = filterTrace(data, samplingRate, freqRange)

[x y z] = size(data);
[b,a]=butter(2,freqRange/(samplingRate/2));

for n = 1:y
    for nn=1:z
        filtData(:,n,nn)=filtfilt(b,a,data(:,n,nn));
    end
end

return;