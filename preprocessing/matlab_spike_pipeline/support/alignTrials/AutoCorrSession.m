function AutoCorrSession(path,fileName, onlyRun, mazeSess)
% extract autocorrelation features for the selected session
%
% e.g.: detectInt('./','A111-20150301-01_DataStructure_mazeSection1_TrialType1')
% 
% by Yingxue, 08/24/2017

    %%%%%%%% check arguments
    if nargin<4
        disp('At least two arguments are needed for this function.');
        return;
    elseif nargin > 4
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
    fullPath = [path fileNameThetaMod];
    if(exist(fullPath) == 0)
        disp('_ThetaMod_Run.mat file does not exist.');
        return;
    end
    load(fullPath,'thetaModSess'); 
            
    fileNameCCG = [fileName '_CCG_Run' num2str(onlyRun) '.mat'];
    fullPath = [path fileNameCCG];
    if(exist(fullPath) == 0)
        disp('CCG_Run.mat file does not exist.');
        return;
    end
    load(fullPath,'CCGSess'); 
    
    fileNameRec = [fileName '.mat'];
    
    fullPath = [path fileNameRec];
    if(exist(fullPath) == 0)
        disp('File does not exist.');
        return;
    end
    load(fullPath,'cluList');    
    
    if(mazeSess == 0)
        mazeSess = 1;
    end
    ccgT = CCGSess{mazeSess}.ccgT;
    numNeurons = length(cluList.all);    
    indMean = find(abs(ccgT) > 100);    
    indPeak = find(abs(ccgT) < 50);
    ind40ms = find(abs(ccgT) >= 40 & abs(ccgT) <= 50);
    peakTimes = ccgT(indPeak);
    
    % generate Gaussian kernel
    std = 1.5;
    paramT.gaussFilt = gaussFilter(12*std, std);
    lenGaussKernel = length(paramT.gaussFilt);
    normFactor = sum(paramT.gaussFilt);
    paramT.gaussFilt = paramT.gaussFilt./normFactor;
    
    autoCorr = struct('indMean', indMean',...  % indices used to calculate mean autocorrelogram
                  'indPeak', indPeak',...  % indices used to calculate peak autocorrelogram
                  'ind40ms', ind40ms',...  % indices at 20 ms of the autocorelogram
                  'mean',zeros(1,numNeurons),... % mean of autocorrelogram
                  'peakAmp',zeros(1,numNeurons),... % peak amplitude of autocorrelogram
                  'peakTime',zeros(1,numNeurons),... % time of the peak of autocorrelogram
                  'peakToMean',zeros(1,numNeurons),... % peak to mean ratio of autocorrelogram
                  'peakTo40ms',zeros(1,numNeurons),... % peak to 40 ms ratio of autocorrelogram
                  ...
                  'realMean',zeros(1,numNeurons),... % mean of autocorrelogram
                  'realPeakAmp',zeros(1,numNeurons),... % peak amplitude of autocorrelogram
                  'realPeakTime',zeros(1,numNeurons),... % time of the peak of autocorrelogram
                  'realPeakToMean',zeros(1,numNeurons),... % peak to mean ratio of autocorrelogram
                  'realPeakTo40ms',zeros(1,numNeurons),... % peak to 40 ms ratio of autocorrelogram
                  ...
                  'thetaMod',zeros(1,numNeurons),... % theta modulation calculated based on ACG spectrum
                  'thetaModInd',zeros(1,numNeurons)); % theta modulation index
        
    for i = 1:numNeurons
        autocor = conv(CCGSess{mazeSess}.ccgVal(:,i,i),paramT.gaussFilt);
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
        
        autoCorr.realMean(i) = mean(CCGSess{mazeSess}.ccgVal(indMean,i,i));
        [autoCorr.realPeakAmp(i), indPeakCorr] = max(CCGSess{mazeSess}.ccgVal(indPeak,i,i));
        autoCorr.realAmp40ms(i) = mean(CCGSess{mazeSess}.ccgVal(ind40ms,i,i));
        autoCorr.realPeakTime(i) = abs(peakTimes(indPeakCorr(1)));
        
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
        
        autoCorr.thetaMod(i) = thetaModSess{mazeSess}.thetaMod(i);
        autoCorr.thetaModInd(i) = thetaModSess{mazeSess}.thetaModInd(i);
    end
    
    fileNameAutoCorr = [fileName '_AutoCorr_Run' num2str(onlyRun) '_mSess' num2str(mazeSess) '.mat'];
    fullPath = [path fileNameAutoCorr];
    save(fullPath,'autoCorr')
end
