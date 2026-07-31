function [filtDataRMSstruct] = getThetaRms(data, samplFq, freqRange, ...
                                            indLap, indSpeed, indTrace)     

if(nargin<6)
    indTrace = 1;
end

filtDataRMSstruct = [];
numTrials = length(indLap);
if(numTrials == 0)
    return;
end

if(~isempty(indTrace))
    numEegTrace = length(indTrace);
else
    numEegTrace = 1;
    indTrace = 1;
end

filtDataRMSstruct = struct('indLap',indLap,...
    'indTrace', indTrace,...
    'RMSperTrial',zeros(numEegTrace,numTrials),...
    'meanRMS',zeros(numEegTrace),...
    'RMSperTrialFL',zeros(numEegTrace,numTrials),...
    'meanRMSFL',zeros(numEegTrace));

for j = 1:numEegTrace
    for i = 1:numTrials  % calculate the RMS for each trial and then average over trials
        
        eeg = data{indLap(i)}(indSpeed{i},indTrace(j));
        filtData = filterTrace(eeg, samplFq, freqRange);
        filtDataRMSstruct.RMSperTrial(j,i) = rms(filtData);
        eegFullLen = data{indLap(i)}(:,indTrace(j));
        filtDataFullLen = filterTrace(eegFullLen, samplFq, freqRange);
        filtDataRMSstruct.RMSperTrialFL(j,i) = rms(filtDataFullLen);
        
    end
    %     figure;
    %     plot(eegFullLen,'k')
    %     hold on;plot(filtDataFullLen,'g')
    %     title(num2str(rms(filtDataFullLen)));
    
    filtDataRMSstruct.meanRMS(j) = mean(filtDataRMSstruct.RMSperTrial(j,:));
    filtDataRMSstruct.meanRMSFL(j) = mean(filtDataRMSstruct.RMSperTrialFL(j,:));
end

return