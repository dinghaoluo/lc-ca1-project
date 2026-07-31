function ThetaModulationAlignedRun(path,fileName,onlyRun,mazeSess,figState)
% calculate theta modulation of each neuron
%
% by Yingxue 08/25/2017

     
    %%%%%%%%% check arguments
    if nargin<2
        disp('At least three arguments are needed for this function.');
        return;
    elseif nargin == 2
        onlyRun = 1;
        mazeSess = 1;
        figState = 0;
    elseif nargin == 3
        mazeSess = 1;
        figState = 0;
    elseif nargin == 4
        figState = 0;
    elseif nargin > 5
        disp('Too many input arguments');        
        return;
    end
        
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    fileNameThetaMod = [fileName '_ThetaModAlignedRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    fileNameCCG = [fileName '_CCGAlignedRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    fileNameOrig = [fileName '.mat'];
    
    fullPath = [path fileNameOrig]; 
    if(exist(fullPath) == 0)
        disp('The file does not exist');
        return;
    end
    load(fullPath,'cluList');
    totClu = length(cluList.all);
    
    fullPath = [path fileNameCCG];
    if(exist(fullPath) == 0)
        disp(['Extended file does not exist. Please run function',...
              ' CalCCGAlignedRun first. onlyRun has to be equal to 1.']);
        return;
    end
    load(fullPath);
    
    fullPath = [path fileName '_PeakFRAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];    
    if(exist(fullPath) == 0)
        disp('The _PeakFRAligned file does not exist');
        return;
    end
    load(fullPath,'pulseMeth');
                       
    %%%%%%%%% initialize constants
    GlobalConst;
   
    % generate Gaussian kernel for filtering ACG
    std = 5;
    paramT.gaussFilt = gaussFilter(12*std, std);
    lenGaussKernel = length(paramT.gaussFilt);
    normFactor = sum(paramT.gaussFilt);
    paramT.gaussFilt = paramT.gaussFilt./normFactor;
    
    paramT.thetaRange = [6 12];
    paramT.deltaRange = [1.5 4];
    paramT.betaRange = [13 30];
    paramT.gammaRange = [30 60];
    paramT.highgammaRange = [60 150];
    
    % generate Gaussian kernel for filtering ACG spetrum
    std = 15;
    paramT.gaussFiltSpect = gaussFilter(12*std, std);
    lenGaussKernelSpect = length(paramT.gaussFiltSpect);
    normFactorSpect = sum(paramT.gaussFiltSpect);
    paramT.gaussFiltSpect = paramT.gaussFiltSpect./normFactorSpect;
    
    thetaMod = struct('ACGSpectrum',[],... % ACG spectrum
                      'ACGSpectrumFreq',[],... % ACG spectrum frequency
                      'ACGSpectrumDeltaPeak',[],... % ACG spectrum delta peak
                      'ACGSpectrumThetaPeak',[],... % ACG spectrum theta peak
                      'ACGSpectrumBetaPeak',[],... % ACG spectrum beta peak
                      'ACGSpectrumGammaPeak',[],... % ACG spectrum gamma peak
                      'ACGSpectrumHighGammaPeak',[],... % ACG spectrum high gamma peak
                      ...
                      'filACGSpectrumDeltaPeak',[],... % filtered ACG spectrum delta peak
                      'filACGSpectrumThetaPeak',[],... % filtered ACG spectrum theta peak
                      'filACGSpectrumBetaPeak',[],... % filtered ACG spectrum beta peak
                      'filACGSpectrumGammaPeak',[],... % filtered ACG spectrum gamma peak
                      'filACGSpectrumHighGammaPeak',[],... % filtered ACG spectrum high gamma peak
                      ...
                      'filteredACGSpectrum',[],... % filtered ACG spectrum
                      'thetaMod',zeros(1,totClu),... % theta modulation calculated based on ACG spectrum
                      'filteredACG',[],... % smoothed ACG
                      'trough', zeros(1,totClu),... % ACG first trough (method 2)
                      'peak', zeros(1,totClu),... % ACG first peak (method 2)
                      'thetaModInd', zeros(1,totClu),... % theta modulation index (method 2)
                      ...
                      'trough3', zeros(1,totClu),... % ACG first trough (method 3)
                      'peak3', zeros(1,totClu),... % ACG first peak (method 3)
                      'troughT3', zeros(1,totClu),... % index of ACG first trough (method 3)
                      'peakT3', zeros(1,totClu),... % index of ACG first peak (method 3)
                      'thetaModInd3', zeros(1,totClu)); % theta modulation index (method 3)
                  
    %% calculate theta modulation for non-stimulated good trials
    disp('Calculate theta modulation for non-stimulated good trials')
    thetaModNonStimGood = thetaMod;
    for j = 1:totClu      
        [thetaModNonStimGood.ACGSpectrum(j,:),...
            thetaModNonStimGood.ACGSpectrumFreq(j,:),...
            thetaModNonStimGood.thetaMod(j)] = ...
            thetaMod1(CCGNonStimGood.ccgVal(:,j,j),CCGNonStimGood.ccgT);

        %% calculate peaks of the ACG spectrum
        % calculate the peak within delta frequencies of ACG spectrum 
        thetaModNonStimGood.ACGSpectrumDeltaPeak(j) = ...
            peakSpectrum(paramT.deltaRange(1),paramT.deltaRange(2),...
                        thetaModNonStimGood.ACGSpectrumFreq(j,:),...
                        thetaModNonStimGood.ACGSpectrum(j,:));

        % calculate the peak within theta frequencies of ACG spectrum 
        thetaModNonStimGood.ACGSpectrumThetaPeak(j) = ...
            peakSpectrum(paramT.thetaRange(1),paramT.thetaRange(2),...
                        thetaModNonStimGood.ACGSpectrumFreq(j,:),...
                        thetaModNonStimGood.ACGSpectrum(j,:));

        % calculate the peak within beta frequencies of ACG spectrum 
        thetaModNonStimGood.ACGSpectrumBetaPeak(j) = ...
            peakSpectrum(paramT.betaRange(1),paramT.betaRange(2),...
                        thetaModNonStimGood.ACGSpectrumFreq(j,:),...
                        thetaModNonStimGood.ACGSpectrum(j,:));

        % calculate the peak within gamma frequencies of ACG spectrum 
        thetaModNonStimGood.ACGSpectrumGammaPeak(j) = ...
            peakSpectrum(paramT.gammaRange(1),paramT.gammaRange(2),...
                        thetaModNonStimGood.ACGSpectrumFreq(j,:),...
                        thetaModNonStimGood.ACGSpectrum(j,:));

        % calculate the peak within high gamma frequencies of ACG spectrum 
        thetaModNonStimGood.ACGSpectrumHighGammaPeak(j) = ...
            peakSpectrum(paramT.highgammaRange(1),paramT.highgammaRange(2),...
                        thetaModNonStimGood.ACGSpectrumFreq(j,:),...
                        thetaModNonStimGood.ACGSpectrum(j,:));


        %% calculate peaks of the filtered ACG spectrum

        filteredACGSpectrum = conv(thetaModNonStimGood.ACGSpectrum(j,:),...
                            paramT.gaussFiltSpect);
        if(mod(lenGaussKernelSpect,2) == 0)
            filteredACGSpectrum = ...
                filteredACGSpectrum(floor(lenGaussKernelSpect/2)+1:...
                    (end-floor(lenGaussKernelSpect/2))+1); 
        else
            filteredACGSpectrum = ...
                filteredACGSpectrum(floor(lenGaussKernelSpect/2)+1:...
                    (end-floor(lenGaussKernelSpect/2))); 
        end
        thetaModNonStimGood.filteredACGSpectrum(j,:) = ...
                filteredACGSpectrum;

        % calculate the peak within delta frequencies of the filtered
        % ACG spectrum
        thetaModNonStimGood.filACGSpectrumDeltaPeak(j) = ...
            peakSpectrum(paramT.deltaRange(1),paramT.deltaRange(2),...
                        thetaModNonStimGood.ACGSpectrumFreq(j,:),...
                        filteredACGSpectrum);

        % calculate the peak within theta frequencies of the filtered
        % ACG spectrum
        thetaModNonStimGood.filACGSpectrumThetaPeak(j) = ...
            peakSpectrum(paramT.thetaRange(1),paramT.thetaRange(2),...
                        thetaModNonStimGood.ACGSpectrumFreq(j,:),...
                        filteredACGSpectrum);

        % calculate the peak within beta frequencies of the filtered
        % ACG spectrum
        thetaModNonStimGood.filACGSpectrumBetaPeak(j) = ...
            peakSpectrum(paramT.betaRange(1),paramT.betaRange(2),...
                        thetaModNonStimGood.ACGSpectrumFreq(j,:),...
                        filteredACGSpectrum);

        % calculate the peak within gamma frequencies of the filtered
        % ACG spectrum
        thetaModNonStimGood.filACGSpectrumGammaPeak(j) = ...
            peakSpectrum(paramT.gammaRange(1),paramT.gammaRange(2),...
                        thetaModNonStimGood.ACGSpectrumFreq(j,:),...
                        filteredACGSpectrum);

        % calculate the peak within high gamma frequencies of the filtered
        % ACG spectrum
        thetaModNonStimGood.filACGSpectrumHighGammaPeak(j) = ...
            peakSpectrum(paramT.highgammaRange(1),paramT.highgammaRange(2),...
                        thetaModNonStimGood.ACGSpectrumFreq(j,:),...
                        filteredACGSpectrum);

        filteredACG = conv(CCGNonStimGood.ccgVal(:,j,j),paramT.gaussFilt);
        if(mod(lenGaussKernel,2) == 0)
            filteredACG = ...
                filteredACG(floor(lenGaussKernel/2)+1:...
                    (end-floor(lenGaussKernel/2))+1); 
        else
            filteredACG = ...
                filteredACG(floor(lenGaussKernel/2)+1:...
                    (end-floor(lenGaussKernel/2))); 
        end
        [thetaModNonStimGood.trough(j),thetaModNonStimGood.peak(j),...
            thetaModNonStimGood.thetaModInd(j)] = ...
            thetaMod2(filteredACG,CCGNonStimGood.ccgT);
        thetaModNonStimGood.filteredACG(j,:) = filteredACG;

        [thetaModNonStimGood.trough3(j),thetaModNonStimGood.peak3(j),...
            thetaModNonStimGood.troughT3(j),thetaModNonStimGood.peakT3(j),...
            thetaModNonStimGood.thetaModInd3(j)] ...
            = thetaMod3(filteredACG,CCGNonStimGood.ccgT);
    end
    
    %% calculate theta modulation for non-stimulated bad trials
    disp('Calculate theta modulation for non-stimulated bad trials')
    thetaModNonStimBad = thetaMod;
    for j = 1:totClu      
        [thetaModNonStimBad.ACGSpectrum(j,:),...
            thetaModNonStimBad.ACGSpectrumFreq(j,:),...
            thetaModNonStimBad.thetaMod(j)] = ...
            thetaMod1(CCGNonStimBad.ccgVal(:,j,j),CCGNonStimBad.ccgT);
        
        filteredACG = conv(CCGNonStimBad.ccgVal(:,j,j),paramT.gaussFilt);
        if(mod(lenGaussKernel,2) == 0)
            filteredACG = ...
                filteredACG(floor(lenGaussKernel/2)+1:...
                    (end-floor(lenGaussKernel/2))+1); 
        else
            filteredACG = ...
                filteredACG(floor(lenGaussKernel/2)+1:...
                    (end-floor(lenGaussKernel/2))); 
        end
        
        [thetaModNonStimBad.trough(j),thetaModNonStimBad.peak(j),...
            thetaModNonStimBad.thetaModInd(j)] = ...
            thetaMod2(filteredACG,CCGNonStimBad.ccgT);
        thetaModNonStimBad.filteredACG(j,:) = filteredACG;

        [thetaModNonStimBad.trough3(j),thetaModNonStimBad.peak3(j),...
            thetaModNonStimBad.troughT3(j),thetaModNonStimBad.peakT3(j),...
            thetaModNonStimBad.thetaModInd3(j)] ...
            = thetaMod3(filteredACG,CCGNonStimBad.ccgT);        
    end
    
    %% added by Yingxue on 2/15/2021
    %% calculate theta modulation for stimulated trials 
    thetaModStim = [];
    thetaModStimCtrl = [];
    for i = 1:length(pulseMeth)
        disp('Calculate theta modulation for stimulated trials')
        thetaModStim{i} = thetaMod;
        for j = 1:totClu      
            [thetaModStim{i}.ACGSpectrum(j,:),...
                thetaModStim{i}.ACGSpectrumFreq(j,:),...
                thetaModStim{i}.thetaMod(j)] = ...
                thetaMod1(CCGStim{i}.ccgVal(:,j,j),CCGStim{i}.ccgT);

            filteredACG = conv(CCGStim{i}.ccgVal(:,j,j),paramT.gaussFilt);
            if(mod(lenGaussKernel,2) == 0)
                filteredACG = ...
                    filteredACG(floor(lenGaussKernel/2)+1:...
                        (end-floor(lenGaussKernel/2))+1); 
            else
                filteredACG = ...
                    filteredACG(floor(lenGaussKernel/2)+1:...
                        (end-floor(lenGaussKernel/2))); 
            end

            [thetaModStim{i}.trough(j),thetaModStim{i}.peak(j),...
                thetaModStim{i}.thetaModInd(j)] = ...
                thetaMod2(filteredACG,CCGStim{i}.ccgT);
            thetaModStim{i}.filteredACG(j,:) = filteredACG;

            [thetaModStim{i}.trough3(j),thetaModStim{i}.peak3(j),...
                thetaModStim{i}.troughT3(j),thetaModStim{i}.peakT3(j),...
                thetaModStim{i}.thetaModInd3(j)] ...
                = thetaMod3(filteredACG,CCGStim{i}.ccgT);
        end
        
        disp('Calculate theta modulation for control trials during stimulation')
        thetaModStimCtrl{i} = thetaMod;
        for j = 1:totClu      
            [thetaModStimCtrl{i}.ACGSpectrum(j,:),...
                thetaModStimCtrl{i}.ACGSpectrumFreq(j,:),...
                thetaModStimCtrl{i}.thetaMod(j)] = ...
                thetaMod1(CCGStimCtrl{i}.ccgVal(:,j,j),CCGStimCtrl{i}.ccgT);

            filteredACG = conv(CCGStimCtrl{i}.ccgVal(:,j,j),paramT.gaussFilt);
            if(mod(lenGaussKernel,2) == 0)
                filteredACG = ...
                    filteredACG(floor(lenGaussKernel/2)+1:...
                        (end-floor(lenGaussKernel/2))+1); 
            else
                filteredACG = ...
                    filteredACG(floor(lenGaussKernel/2)+1:...
                        (end-floor(lenGaussKernel/2))); 
            end

            [thetaModStimCtrl{i}.trough(j),thetaModStimCtrl{i}.peak(j),...
                thetaModStimCtrl{i}.thetaModInd(j)] = ...
                thetaMod2(filteredACG,CCGStimCtrl{i}.ccgT);
            thetaModStimCtrl{i}.filteredACG(j,:) = filteredACG;

            [thetaModStimCtrl{i}.trough3(j),thetaModStimCtrl{i}.peak3(j),...
                thetaModStimCtrl{i}.troughT3(j),thetaModStimCtrl{i}.peakT3(j),...
                thetaModStimCtrl{i}.thetaModInd3(j)] ...
                = thetaMod3(filteredACG,CCGStimCtrl{i}.ccgT);
        end
    end
        
    fullPath = [path fileNameThetaMod];
    save(fullPath, 'thetaModNonStimGood','thetaModNonStimBad',...
        'thetaModStim','thetaModStimCtrl','paramT');
    
    if(figState ~= 0)
        numSess = length(mazeSess);
        meanThetaMod = zeros(1,numSess);
        stdThetaMod = zeros(1,numSess);
        exc = autoCorr.isPyrneuron;
        for i = 1:numSess
            thetaMod = thetaModNonStimGood.thetaMod(exc);
            valid = isnan(thetaMod) == 0;
            meanThetaMod(i) = mean(thetaMod(valid));
            stdThetaMod(i) = std(thetaMod(valid));
        end
        barPlot(1:length(mazeSess),meanThetaMod,stdThetaMod,...
                'Session','Mean theta modulation',...
                'Mean theta modulation per session');
    end
                   
    clear mydata;

end

function peakFreq = peakSpectrum(startFreq,stopFreq,spectrumFreq,...
                                spectrum)
    indFreq = find(spectrumFreq >= startFreq...
                & spectrumFreq <= stopFreq);
    [~,indPeak] = max(spectrum(indFreq));
    peakFreq = spectrumFreq(indFreq(indPeak(1)));
end
