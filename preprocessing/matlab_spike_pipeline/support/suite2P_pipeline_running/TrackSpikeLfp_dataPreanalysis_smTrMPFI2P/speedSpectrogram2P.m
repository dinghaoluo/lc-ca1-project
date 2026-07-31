function speedSpectro = speedSpectrogram(speed,indLaps,numTrials,param)
% Calculate the spectrogram of speed per trial

    speedSpectro = struct(...
        'lapList',indLaps',... % trials involved
        'spectrogram',[],... % spectrogram of speed in each trial
        'timeSpectro',[],... % time instants at which the spectrogram is computed
        'freqSpectro',[],... % frequencies of the spectrogram
        'maxAmpSpectro3_5',[],... % max power over 3-5 Hz changing over time
        'maxAmpSpectro6_10',[],... % max power over 6-10 Hz changing over time
        'maxFreqSpectro3_5',[],... % frequency of the max power over 3-5 hz chaning over time
        'maxFreqSpectro6_10',[],... % frequency of the max power over 6-10 hz chaning over time
        'meanAmpSpectro3_5',[],... % mean power over 3-5 Hz changing over time
        'meanAmpSpectro6_10',[],... % mean power over 6-10 Hz changing over time
        'mmeanAmpSpectro3_5',[],... % mean power over 3-5 Hz for each trial
        'mmeanAmpSpectro6_10',[],... % mean power over 6-10 Hz for each trial
        'speedMovingAvg',[],... % speed moving avg
        'meanPower3_5Hz',[],... % mean power for 3-5 Hz
        'stdPower3_5Hz',[],... % std power for 3-5 Hz
        'meanPower6_10Hz',[],... % mean power for 6-10 Hz
        'stdPower6_10Hz',[],... % std power for 6-10 Hz
        'lengthHighPower3_5',[],... % total number of samples with power larger than threshold
        'lengthHighPower6_10',[],... % total number of samples with power larger than threshold
        'percHighPower3_5',[],... % percent of time in the trial with high power
        'percHighPower6_10',[]); % percent of time in the trial with high power
    
    tr = 0;
    meanPowerAll3_5 = [];
    meanPowerAll6_10 = [];
    for i = 1:numTrials
        indTr = find(indLaps == i);
        if(isempty(indTr))
            continue;
        end
        tr = tr+1;
        speedTmp = speed{i};
        if(length(speedTmp) < param.sampleFq*100)
            [spectroTmp,freq,speedSpectro.timeSpectro{tr}] = ...
                spectrogram(speedTmp,param.spectroWin,param.spectroOverlap,...
                            2^nextpow2(length(speedTmp)),param.sampleFq);
            speedSpectro.spectrogram{tr} = spectroTmp;

            indPower3_5 = find(freq >= 3 & freq <= 5);
            spectro3_5 = spectroTmp(indPower3_5,:);
            ampSpectro3_5 = abs(spectro3_5).^2;
            [maxAmpSpecto3_5,maxIndSpectro3_5] = max(ampSpectro3_5,[],1);
            maxFreqSpectro3_5 = freq(indPower3_5(maxIndSpectro3_5));
            meanAmpSpectro3_5 = mean(ampSpectro3_5,1);

            indPower6_10 = find(freq >=6 & freq <= 10);
            spectro6_10 = spectroTmp(indPower6_10,:);
            ampSpectro6_10 = abs(spectro6_10).^2;
            [maxAmpSpecto6_10,maxIndSpectro6_10] = max(ampSpectro6_10,[],1);
            maxFreqSpectro6_10 = freq(indPower6_10(maxIndSpectro6_10));
            meanAmpSpectro6_10 = mean(ampSpectro6_10,1);

            speedSpectro.maxAmpSpectro3_5{tr} = maxAmpSpecto3_5;
            speedSpectro.maxFreqSpectro3_5{tr} = maxFreqSpectro3_5;
            speedSpectro.meanAmpSpectro3_5{tr} = meanAmpSpectro3_5;
            speedSpectro.mmeanAmpSpectro3_5(tr) = mean(meanAmpSpectro3_5);
            meanPowerAll3_5 = [meanPowerAll3_5 meanAmpSpectro3_5];

            speedSpectro.maxAmpSpectro6_10{tr} = maxAmpSpecto6_10;
            speedSpectro.maxFreqSpectro6_10{tr} = maxFreqSpectro6_10;
            speedSpectro.meanAmpSpectro6_10{tr} = meanAmpSpectro6_10;
            speedSpectro.mmeanAmpSpectro6_10(tr) = mean(meanAmpSpectro6_10);
            meanPowerAll6_10 = [meanPowerAll6_10 meanAmpSpectro6_10];

            movingAvg = ...
                zeros(1,length(speedSpectro.timeSpectro{tr}));
            for j = 1:length(speedSpectro.timeSpectro{tr})
                indStart = round(speedSpectro.timeSpectro{tr}(j)*param.sampleFq...
                            -param.spectroWin/2);
                if(indStart < 1)
                    indStart = 1;
                end
                indEnd = round(speedSpectro.timeSpectro{tr}(j)*param.sampleFq...
                            +param.spectroWin/2);
                if(indEnd > length(speedTmp))
                    indEnd = length(speedTmp);
                end
                movingAvg(j) = mean(speedTmp(indStart:indEnd));
            end
            speedSpectro.speedMovingAvg{tr} = movingAvg;
        else
            speedSpectro.spectrogram{tr} = spectroTmp;
            
        end
    end
    speedSpectro.freqSpectro = freq;
    speedSpectro.indPower3_5 = indPower3_5;
    speedSpectro.indPower6_10 = indPower6_10;
    speedSpectro.meanPower3_5Hz = mean(meanPowerAll3_5);
    speedSpectro.stdPower3_5Hz = std(meanPowerAll3_5);
    speedSpectro.meanPower6_10Hz = mean(meanPowerAll6_10);
    speedSpectro.stdPower6_10Hz = std(meanPowerAll6_10);
    
    speedSpectro.threPower3_5Hz = speedSpectro.meanPower3_5Hz + ...
                0*speedSpectro.stdPower3_5Hz;
    speedSpectro.threPower6_10Hz = speedSpectro.meanPower6_10Hz + ...
                0*speedSpectro.stdPower6_10Hz;
    for i = 1:length(indLaps)
        speedSpectro.lengthHighPower3_5(i) = ...
            sum(speedSpectro.meanAmpSpectro3_5{i} > speedSpectro.threPower3_5Hz);
        speedSpectro.percHighPower3_5(i) = ...
            speedSpectro.lengthHighPower3_5(i)/length(speedSpectro.timeSpectro{tr});
        speedSpectro.lengthHighPower6_10(i) = ...
            sum(speedSpectro.meanAmpSpectro6_10{i} > speedSpectro.threPower6_10Hz);
        speedSpectro.percHighPower6_10(i) = ...
            speedSpectro.lengthHighPower6_10(i)/length(speedSpectro.timeSpectro{tr});
    end
end
