function spikeThetaFreqAmpBefRunOnset(path, fileName, onlyRun, mazeSess)
% plot spike phase and theta frequency over time
    
    fullPath = [path fileName '_ThetaAlignRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _ThetaAlignRun_msess file does not exist');
        return;
    end
    load(fullPath);
        
    fullPath = [path fileName '_ThetaAlignRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat']; 
    save(fullPath,'meanThetaNonStimGood','meanThetaNonStimBad','meanThetaStim','t');
   
    meanTheta.totThetaAmpGoodTr(n,:)
end

function meanTheta =calTheta(trialsRun,trialNoG,trialLenT,nSampBef,sampleFq)
    vLen = length(-nSampBef+1:trialLenT*sampleFq);
    meanTheta.totThetaFreqGoodTr = zeros(length(trialNoG),vLen);
    meanTheta.totThetaAmpGoodTr = zeros(length(trialNoG),vLen);
    meanTheta.totSpeedGoodTr = zeros(length(trialNoG),vLen);
    
    n = 1;
    for j = trialNoG   
        if(isempty(trialsRun.ThetaFreq{j}))
            continue;
        end
        vSamp = [-nSampBef+1:trialsRun.numSamples(j)];
        vFreq = [trialsRun.ThetaFreqBef{j}' trialsRun.ThetaFreq{j}'];
        if(length(vSamp) > vLen)
            meanTheta.totThetaFreqGoodTr(n,:) = vFreq(1:vLen);
        else
            meanTheta.totThetaFreqGoodTr(n,1:length(vSamp)) = vFreq;
        end
        n = n+1;
    end
    
    n = 1;
    for j = trialNoG  
        if(isempty(trialsRun.ThetaFreq{j}))
            continue;
        end
        vSamp = [-nSampBef+1:trialsRun.numSamples(j)]/sampleFq;
        vAmp = [trialsRun.ThetaAmpBef{j}' trialsRun.ThetaAmp{j}'];
        if(length(vSamp) > vLen)
            meanTheta.totThetaAmpGoodTr(n,:) = vAmp(1:vLen);
        else
            meanTheta.totThetaAmpGoodTr(n,1:length(vSamp)) = vAmp;
        end      
        n = n+1;
    end
    
    n = 1;
    for j = trialNoG   
        if(isempty(trialsRun.ThetaFreq{j}))
            continue;
        end
        vSamp = [-nSampBef+1:trialsRun.numSamples(j)]/sampleFq;
        speed = [trialsRun.speed_MMsecBef{j}' trialsRun.speed_MMsec{j}'];
        indSpeed = speed < -1;
        speed(indSpeed) = 0;
        if(length(vSamp) > vLen)
            meanTheta.totSpeedGoodTr(n,:) = speed(1:vLen);
        else
            meanTheta.totSpeedGoodTr(n,1:length(vSamp)) = speed;
        end
        n = n+1;
    end
    
    [b,a] = butter(2,0.01); % low pass filtered phase histogram
    meanTheta.meanThetaFreqGoodTr = mean(meanTheta.totThetaFreqGoodTr);
    if(length(meanTheta.meanThetaFreqGoodTr) > 1)
        meanTheta.meanThetaFreqFiltGoodTr = filtfilt(b,a,meanTheta.meanThetaFreqGoodTr); 
    else
        meanTheta.meanThetaFreqFiltGoodTr = meanTheta.meanThetaFreqGoodTr;
    end
    meanTheta.stdThetaFreqGoodTr = std(meanTheta.totThetaFreqGoodTr)/sqrt(length(trialNoG));
    meanTheta.meanThetaAmpGoodTr = mean(meanTheta.totThetaAmpGoodTr);
    if(length(meanTheta.meanThetaFreqGoodTr) > 1)
        meanTheta.meanThetaAmpFiltGoodTr = filtfilt(b,a,meanTheta.meanThetaAmpGoodTr); 
    else
        meanTheta.meanThetaAmpFiltGoodTr = meanTheta.meanThetaAmpGoodTr;
    end
    meanTheta.stdThetaAmpGoodTr = std(meanTheta.totThetaAmpGoodTr)/sqrt(length(trialNoG));
    meanTheta.meanSpeedGoodTr = mean(meanTheta.totSpeedGoodTr);
    meanTheta.stdSpeedGoodTr = std(meanTheta.totSpeedGoodTr)/sqrt(length(trialNoG));
end
