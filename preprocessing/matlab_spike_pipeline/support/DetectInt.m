function DetectInt(path,fileName, onlyRun)
% detect interneurons based on their autocorrelogram
%
% e.g.: detectInt('./','A111-20150301-01_DataStructure_mazeSection1_TrialType1')
% 
% by Yingxue, 08/24/2017

    %%%%%%%% check arguments
    if nargin<3
        disp('At least two arguments are needed for this function.');
        return;
    elseif nargin > 3
        disp('Too many arguments');        
        return;
    end

    %%%%%%%%% load recording file
    indexFileName = strfind(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    fileNameInfo = [fileName '_Info.mat'];
    fileNameThetaMod = [fileName '_ThetaMod_Run' num2str(onlyRun) '.mat']; 
    fileName = [fileName '.mat'];
    
    fullPath = [path fileName];
    if(exist(fullPath) == 0)
        disp('File does not exist.');
        return;
    end
    load(fullPath,'cluList');    
        
    indexFileName = strfind(fileName,'_');
    fileNameCCGT = [fileName(1:indexFileName(1)-1)...
                        '_BehavElectrDataLFP_CCG.mat'];
    fullPath = [path fileNameCCGT];
    if(exist(fullPath) == 0)
        disp('_BehavElectrDataLFP_CCG.mat file does not exist.');
        return;
    end
    load(fullPath,'ccgT','ccgVal'); 
    
    fileNameCCGT = [fileName(1:indexFileName(1)-1)...
                        '_BehavElectrDataLFP.mat'];
    fullPath = [path fileNameCCGT];
    if(exist(fullPath) == 0)
        disp('_BehavElectrDataLFP.mat file does not exist.');
        return;
    end
    load(fullPath,'Clu'); 
    
    fileNameDepth = [fileName(1:end-4) '_Depth.mat'];
    fullPath = [path fileNameDepth];
    if(exist(fullPath) == 0)
        disp('_Depth.mat file does not exist.');
        return;
    end
    load(fullPath); 
    
    GlobalConst;

    paramInt.peakTimeThr = 9; % threshold for the peak time of autocorrelogram
    paramInt.peak2MeanThr = 5; % threshold for the peak to mean ratio of autocorrelogram
    paramInt.peak2MeanThr1 = 2; % threshold for the peak to mean ratio of autocorrelogram
    paramInt.maxFR1 = 1.5; % firing rate threshold for low firing rate intern neurons
    paramInt.minFRInt = 0.3;  % minimum firing rate of interneurons
    paramInt.maxSpkWidthR = 0.65; % max spike width on the right side

    numNeurons = length(cluList.all);    
    indMean = find(abs(ccgT) > 100);    
    indPeak = find(abs(ccgT) < 50);
    ind40ms = find(abs(ccgT) >= 40 & abs(ccgT) <= 50);
    ind10ms = find(abs(ccgT) <= 10);
    ind0ms = find(ccgT == 0);
    peakTimes = ccgT(indPeak);
    
    % generate Gaussian kernel
    std1 = 1.5;
    paramT.gaussFilt = gaussFilter(12*std1, std1);
    lenGaussKernel = length(paramT.gaussFilt);
    normFactor = sum(paramT.gaussFilt);
    paramT.gaussFilt = paramT.gaussFilt./normFactor;
    
    autoCorr = struct('indMean', indMean',...  % indices used to calculate mean autocorrelogram
                  'indPeak', indPeak',...  % indices used to calculate peak autocorrelogram
                  'ind40ms', ind40ms',...  % indices at 20 ms of the autocorelogram
                  'refract',zeros(1,numNeurons),... % refractory period
                  'mean',zeros(1,numNeurons),... % mean of autocorrelogram
                  'peakAmp',zeros(1,numNeurons),... % peak amplitude of autocorrelogram
                  'peakTime',zeros(1,numNeurons),... % time of the peak of autocorrelogram
                  'peakToMean',zeros(1,numNeurons),... % peak to mean ratio of autocorrelogram
                  'peakTo40ms',zeros(1,numNeurons),... % peak to 20 ms ratio of autocorrelogram
                  'isInterneuron',zeros(1,numNeurons),... % is the neuron an interneuron?
                  'isPyrneuron',zeros(1,numNeurons),... % is the neuron a pyramidal neuron?
                  'relDepthNeuHDef',zeros(1,numNeurons),... % the location of the neuron in the pyramidal layer
                  ...
                  'realMean',zeros(1,numNeurons),... % mean of autocorrelogram
                  'realPeakAmp',zeros(1,numNeurons),... % peak amplitude of autocorrelogram
                  'realPeakTime',zeros(1,numNeurons),... % time of the peak of autocorrelogram
                  'realPeakToMean',zeros(1,numNeurons),... % peak to mean ratio of autocorrelogram
                  'realPeakTo40ms',zeros(1,numNeurons),... % peak to 20 ms ratio of autocorrelogram
                  'burstInd',zeros(1,numNeurons)); % burst index
    
    for i = 1:numNeurons
        autocor = conv(ccgVal(:,i,i),paramT.gaussFilt);
        if(mod(lenGaussKernel,2) == 0)
            autocor = ...
                autocor(floor(lenGaussKernel/2)+1:...
                    (end-floor(lenGaussKernel/2))+1); 
        else
            autocor = ...
                autocor(floor(lenGaussKernel/2)+1:...
                    (end-floor(lenGaussKernel/2))); 
        end
        autoCorr.mean(i) = mean(autocor(indMean));
        [autoCorr.peakAmp(i), indPeakCorr] = max(autocor(indPeak));
        autoCorr.amp40ms(i) = mean(autocor(ind40ms));
        autoCorr.peakTime(i) = abs(peakTimes(indPeakCorr(1)));
        
        autoCorr.realMean(i) = mean(ccgVal(indMean,i,i));
        [autoCorr.realPeakAmp(i), indPeakCorr] = max(ccgVal(indPeak,i,i));
        autoCorr.realAmp40ms(i) = mean(ccgVal(ind40ms,i,i));
        autoCorr.realPeakTime(i) = abs(peakTimes(indPeakCorr(1)));
        
        indPeakTmp1 = find(ccgT >= 0 & ccgT <= autoCorr.realPeakTime(i));
        timeTmp1 = ccgT(indPeakTmp1);
        diffCorr = diff(ccgVal(indPeakTmp1,i,i));
        stdDiffCorr = std(diffCorr);
        indRef1 = find(diffCorr >= stdDiffCorr,1,'first');
        indPeakTmp2 = find(ccgT <= 0 & ccgT >= -autoCorr.realPeakTime(i));
        timeTmp2 = ccgT(indPeakTmp2);
        diffCorr = diff(ccgVal(indPeakTmp2(end:-1:1),i,i));
        stdDiffCorr = std(diffCorr);
        indRef2 = find(diffCorr >= stdDiffCorr,1,'first');
        autoCorr.refract(i) = mean(1+abs([timeTmp1(indRef1) timeTmp2(length(timeTmp2)-indRef2+1)]));
        
        if(autoCorr.mean(i) ~= 0)
            autoCorr.peakToMean(i) = autoCorr.peakAmp(i)/autoCorr.mean(i);
            autoCorr.realPeakToMean(i) = autoCorr.realPeakAmp(i)/autoCorr.realMean(i);
        else
            autoCorr.peakToMean(i) = autoCorr.peakAmp(i)/0.01;
            autoCorr.realPeakToMean(i) = autoCorr.realPeakAmp(i)/0.01;
        end
        if(autoCorr.amp40ms(i) ~= 0)
            autoCorr.peakTo40ms(i) = autoCorr.peakAmp(i)/autoCorr.amp40ms(i);
            autoCorr.realPeakTo40ms(i) = autoCorr.realPeakAmp(i)/autoCorr.realAmp40ms(i);
        else
            autoCorr.peakTo40ms(i) = autoCorr.peakAmp(i)/0.01;
            autoCorr.realPeakTo40ms(i) = autoCorr.realPeakAmp(i)/0.01;
        end
        
        peakAmpIn10ms = max(ccgVal(ind10ms,i,i));
        if(peakAmpIn10ms >= autoCorr.realAmp40ms(i))
            autoCorr.burstInd(i) = (peakAmpIn10ms - autoCorr.realAmp40ms(i))/peakAmpIn10ms;
        else
            autoCorr.burstInd(i) = (peakAmpIn10ms - autoCorr.realAmp40ms(i))/autoCorr.realAmp40ms(i);
        end
           
%         if((autoCorr.peakTime(i) > peakTimeThr && ...
%                 autoCorr.peakToMean(i) < peak2MeanThr && ...
%                 cluList.firingRate(i) > minFR) || ...
%                 cluList.firingRate(i) > maxFR)
        if(((autoCorr.peakTime(i) >= paramInt.peakTimeThr && ...
                cluList.firingRate(i) > minFR) || ...
                (autoCorr.peakTo40ms(i) < paramInt.peak2MeanThr && ...
                cluList.firingRate(i) > maxFR) || ...
                (autoCorr.peakTo40ms(i) < paramInt.peak2MeanThr1 && ...
                cluList.firingRate(i) > paramInt.maxFR1)) && ...
                Clu.SpkWidthR(i) < paramInt.maxSpkWidthR)
            autoCorr.isInterneuron(i) = 1;
            
            if(cluList.refracViolPercent(i) >= refracViolPercentThre || ...
                    (cluList.mahalDist(i) <= mahalDistThre && ...
                    cluList.mahalDist(i) >= 0)|| ...
                    cluList.firingRate(i) < paramInt.minFRInt || ...
                    autoCorr.refract(i) < 2)
                autoCorr.isInterneuron(i) = -1;
            end
        end
    end
    
    autoCorr.isPyrneuron = (cluList.firingRate > minFR ...
        & cluList.refracViolPercent < refracViolPercentThre...
        & cluList.mahalDist > mahalDistThre ...
        & Clu.SpkWidthR > paramInt.maxSpkWidthR ...
        & autoCorr.isInterneuron == 0);
    
    maxRefractViol = 1;
    maxCenterMax = -300;
    maxCenterMax1 = -200;
    refViol = cluList.refracViolPercent;
    centerMax = cluList.centerMax;
    leftMax = cluList.leftMax;
    indGoodCenterMax = centerMax-leftMax < maxCenterMax;
    indGoodCenterMaxLow = centerMax-leftMax < maxCenterMax1;
    indGoodRef = refViol < maxRefractViol;
    
    isHighAmp = indGoodRef & indGoodCenterMax == 1; 
    isHighAmpLow = indGoodRef & indGoodCenterMaxLow == 1; 
    autoCorr.relDepthNeuHDef = depthNeu.relDepthNeuHDef;
    autoCorr.isSpikeHighAmp = isHighAmp;
    autoCorr.isSpikeHighAmp200 = isHighAmpLow;
    
    fullPath = [path fileNameInfo];
    save(fullPath,'autoCorr','paramInt','-append')
end
