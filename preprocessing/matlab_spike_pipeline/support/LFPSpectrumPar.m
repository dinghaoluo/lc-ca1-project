function LFPSpectrumStruct = LFPSpectrumPar(eeg,indLapList,...
    thetaFqRange,gammaFqRange,numSamples,sampleFq,whitenOrd,indSpeed,figureState)
%  Calculate the averaged LFP spectrum over trials -- Parallelized version
%
% Inputs:
% eeg:                  eeg structure including eeg at each sampling point of each valid trial
% indLapList:           the array containing the index numbers of the valid trials (Caution: since theta strcuture records the theta traces from all the valid trials in 
%                       the task, indLapList here includes the indices of all to be selected cells in the theta structure, rather than referring back to the 
%                       original trial no.)
% thetaFqRange:         the range of theta frequency band ([lowerbound upperbound])
% gammaFqRange:         the range of gamma frequency band ([lowerbound upperbound])
% numSamples:           number of samples
% sampleFq:             sampling frequency
% whitenOrd:            >0: whiten the signal before spectrum analysis, whitenOrd is the order of AR model 
%                       0: do not whiten the signal
% figureState:          0: figure off
%                       2: figure on
%
% Outputs: LFPSpectrumStruct with following fields
% indLapList: 
% LFPSpetgramPerTrial:  LFP spectrum of each trial
% freqLFP:              frequencies at which the PSD is estimated 
% indexFreqTheta:       index of theta frquencies
% normLFPSpectrumTheta: normalized spectrum in theta band
% halfAmpThetaBand:     the bandwidth in which the amplitude is above half
%                       of the max amplitude within theta band
% indexFreqGamma:       index of gamma frquencies
% normLFPSpectrumGamma: normailzed spectrum in gamma band

LFPSpectrumStruct = [];
numTrials = length(indLapList);
if(numTrials == 0)
    return;
end

numSamplesPadConst = 524288;

if(~isempty(numSamplesPadConst))
    if(numSamplesPadConst > numSamples)
        numSamplesPad = numSamplesPadConst;
    else
        numSamplesPad = max(2^(nextpow2(numSamples)),numSamples);
    end
else
    numSamplesPad = max(2^(nextpow2(numSamples)),numSamples);
end

if(mod(numSamplesPad,2) == 0) % length of the spectrum vector is sampleFq/2/stepFreq
    lenSpectrum = numSamplesPad/2+1;
else
    lenSpectrum = (numSamplesPad+1)/2; 
end

stepFreq = sampleFq/numSamplesPad;
freqLFP = 0:stepFreq:(lenSpectrum-1)*stepFreq;
indexFreqTheta = find(freqLFP <= thetaFqRange(2) ...
                        & freqLFP >= thetaFqRange(1));
indexFreqGamma = find(freqLFP <= gammaFqRange(2) ...
                        & freqLFP >= gammaFqRange(1));

LFPSpectrumStruct = struct('indLapList',indLapList,...
                           'freqLFP',freqLFP',...
                           'whitenOrd',whitenOrd,...
                           'indexFreqTheta',indexFreqTheta,...
                           ...
                           'LFPSpectrumPerTrial',zeros(numTrials,lenSpectrum),...
                           'indMaxThetaBand',zeros(1,numTrials),...
                           'maxThetaFreq',zeros(1,numTrials),...
                           'halfAmpThetaBand',zeros(1,numTrials),...
                           ...
                           'LFPSpectrumPerTrialFL',zeros(numTrials,lenSpectrum),...
                           'indMaxThetaBandFL',zeros(1,numTrials),...
                           'maxThetaFreqFL',zeros(1,numTrials),...
                           'halfAmpThetaBandFL',zeros(1,numTrials));
                       
eegTmp = eeg(indLapList);
disp('Trial ');
tStart = tic;
LFPSpectrumPerTrial = zeros(numTrials,lenSpectrum);
halfAmpThetaBand = zeros(1,numTrials);
indMaxThetaBand = zeros(1,numTrials);
maxThetaFreq = zeros(1,numTrials);

LFPSpectrumPerTrialFL = zeros(numTrials,lenSpectrum);
halfAmpThetaBandFL = zeros(1,numTrials);
indMaxThetaBandFL = zeros(1,numTrials);
maxThetaFreqFL = zeros(1,numTrials);

parfor i = 1:numTrials  % calculate the spectrum for each trial and then average over trials
    fprintf('%d  ', i);
    if(whitenOrd > 0) % whiten the signal 
        whitenEEG = WhiteningNoOrdEst(...
                            eegTmp{i}(indSpeed{i}), whitenOrd, 1);
        whitenEEGFullLen = WhiteningNoOrdEst(...
                            eegTmp{i}, whitenOrd, 1); % use all the data
    else
        whitenEEG = eegTmp{i}(indSpeed{i})-mean(eegTmp{i}(indSpeed{i}));
        whitenEEGFullLen = eegTmp{i}-mean(eegTmp{i});
    end
    
    % spectrum with constant data length
    dataEEG = zeros(1,numSamplesPad);
    dataEEG(1:length(whitenEEG)) = whitenEEG;
    [LFPSpectrumTmp,freqArr] = pmtm(dataEEG,4,numSamplesPad,sampleFq);
    LFPSpectrumPerTrial(i,:) = LFPSpectrumTmp';
    thetaBandLFP = LFPSpectrumTmp(indexFreqTheta);
    maxThetaBand = max(thetaBandLFP);
    indMaxThetaBand(i) = find(thetaBandLFP == maxThetaBand,1);
    maxThetaFreq(i) = indexFreqTheta(indMaxThetaBand(i))...
                                        *stepFreq;
    indBelowHalfAmp = find(thetaBandLFP < maxThetaBand/2);
    belowInd = find(indBelowHalfAmp < indMaxThetaBand(i));
    aboveInd = find(indBelowHalfAmp > indMaxThetaBand(i));
    if(isempty(belowInd))
        belowInd = 1;
    end
    if(isempty(aboveInd))
        aboveInd = length(indBelowHalfAmp);
    end
    halfAmpThetaBand(i) = (indBelowHalfAmp(aboveInd(1)) ...
        - indBelowHalfAmp(belowInd(end)))*stepFreq;
    
    % spectrum with the full length of data
    dataEEG = zeros(1,numSamplesPad);
    dataEEG(1:length(whitenEEGFullLen)) = whitenEEGFullLen;
    [LFPSpectrumTmp,freqArr] = pmtm(dataEEG,4,numSamplesPad,sampleFq);
    LFPSpectrumPerTrialFL(i,:) = LFPSpectrumTmp';
    thetaBandLFP = LFPSpectrumTmp(indexFreqTheta);
    maxThetaBand = max(thetaBandLFP);
    indMaxThetaBandFL(i) = find(thetaBandLFP == maxThetaBand,1);
    maxThetaFreqFL(i) = indexFreqTheta(indMaxThetaBandFL(i))...
                                        *stepFreq;
    indBelowHalfAmp = find(thetaBandLFP < maxThetaBand/2);
    belowInd = find(indBelowHalfAmp < indMaxThetaBand(i));
    aboveInd = find(indBelowHalfAmp > indMaxThetaBand(i));
    if(isempty(belowInd))
        belowInd = 1;
    end
    if(isempty(aboveInd))
        aboveInd = length(indBelowHalfAmp);
    end
    halfAmpThetaBandFL(i) = (indBelowHalfAmp(aboveInd(1)) ...
        - indBelowHalfAmp(belowInd(end)))*stepFreq;
    
%     if(figureState == 2)
%         h = plot(LFPSpectrumStruct.freqLFP(indexFreqTheta),...
%               10*log10(LFPSpectrumStruct.LFPSpectrumPerTrial(i,indexFreqTheta)),...
%              'k');
%         set(h,'LineWidth',1.5);
%         hold on
%         h = plot(LFPSpectrumStruct.freqLFP(indexFreqTheta(indMaxThetaBand)),...
%               10*log10(LFPSpectrumStruct.LFPSpectrumPerTrial...
%                           (i,indexFreqTheta(indMaxThetaBand))),'ro');
%         set(h,'LineWidth',2.0);
%         h = plot(LFPSpectrumStruct.freqLFP...
%                       (indexFreqTheta(indBelowHalfAmp(aboveInd(1)))),...
%                     10*log10(LFPSpectrumStruct.LFPSpectrumPerTrial...
%                       (i,indexFreqTheta(indBelowHalfAmp(aboveInd(1))))),'ro');
%         set(h,'LineWidth',2.0);
%         h = plot(LFPSpectrumStruct.freqLFP...
%                       (indexFreqTheta(indBelowHalfAmp(belowInd(end)))),...
%                     10*log10(LFPSpectrumStruct.LFPSpectrumPerTrial...
%                       (i,indexFreqTheta(indBelowHalfAmp(belowInd(end))))),'ro');
%         set(h,'LineWidth',2.0);
%         hold off
%         set(gca,'FontSize',14.0,'Box','on')
%         xlabel('Frequency (Hz)')
%         ylabel('LFP spectrum');
%         title(['Trial no, ' num2str(i)]);
%         pause(1);
%     end
end
tLapse = toc(tStart);
disp(['End of LFP spectrum calculation, total calculation time: ', num2str(tLapse)]);

LFPSpectrumStruct.LFPSpectrumPerTrial = LFPSpectrumPerTrial;
LFPSpectrumStruct.halfAmpThetaBand = halfAmpThetaBand;
LFPSpectrumStruct.indMaxThetaBand = indMaxThetaBand;
LFPSpectrumStruct.maxThetaFreq = maxThetaFreq;

LFPSpectrumStruct.LFPSpectrumPerTrialFL = LFPSpectrumPerTrialFL;
LFPSpectrumStruct.halfAmpThetaBandFL = halfAmpThetaBandFL;
LFPSpectrumStruct.indMaxThetaBandFL = indMaxThetaBandFL;
LFPSpectrumStruct.maxThetaFreqFL = maxThetaFreqFL;
