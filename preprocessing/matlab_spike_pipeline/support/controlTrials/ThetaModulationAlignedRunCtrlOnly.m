function ThetaModulationAlignedRunCtrlOnly(path,fileName,onlyRun,mazeSess,figState)
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
    fileNameThetaModCtrl = [fileName '_ThetaModAlignedRunCtrl_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    fileNameThetaMod = [fileName '_ThetaModAlignedRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    fileNameCCG = [fileName '_CCGAlignedRunCtrl_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
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
              ' CalCCGAlignedRunCtrlOnly first. onlyRun has to be equal to 1.']);
        return;
    end
    load(fullPath,'CCGNonStim');
    
    fullPath = [path fileNameThetaMod];
    if(exist(fullPath) == 0)
        disp('The _ThetaModAlignedRun file does not exist');
        return;
    end
    load(fullPath,'paramT');
                           
    %%%%%%%%% initialize constants
    GlobalConst;
   
    % generate Gaussian kernel for filtering ACG
    
    thetaModNonStim = struct('ACGSpectrum',[],... % ACG spectrum
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
                  
    %% calculate theta modulation for non-stimulated trials
    disp('Calculate theta modulation for non-stimulated trials')
    for j = 1:totClu      
        [thetaModNonStim.ACGSpectrum(j,:),...
            thetaModNonStim.ACGSpectrumFreq(j,:),...
            thetaModNonStim.thetaMod(j)] = ...
            thetaMod1(CCGNonStim.ccgVal(:,j,j),CCGNonStim.ccgT);

        %% calculate peaks of the ACG spectrum
        % calculate the peak within delta frequencies of ACG spectrum 
        thetaModNonStim.ACGSpectrumDeltaPeak(j) = ...
            peakSpectrum(paramT.deltaRange(1),paramT.deltaRange(2),...
                        thetaModNonStim.ACGSpectrumFreq(j,:),...
                        thetaModNonStim.ACGSpectrum(j,:));

        % calculate the peak within theta frequencies of ACG spectrum 
        thetaModNonStim.ACGSpectrumThetaPeak(j) = ...
            peakSpectrum(paramT.thetaRange(1),paramT.thetaRange(2),...
                        thetaModNonStim.ACGSpectrumFreq(j,:),...
                        thetaModNonStim.ACGSpectrum(j,:));

        % calculate the peak within beta frequencies of ACG spectrum 
        thetaModNonStim.ACGSpectrumBetaPeak(j) = ...
            peakSpectrum(paramT.betaRange(1),paramT.betaRange(2),...
                        thetaModNonStim.ACGSpectrumFreq(j,:),...
                        thetaModNonStim.ACGSpectrum(j,:));

        % calculate the peak within gamma frequencies of ACG spectrum 
        thetaModNonStim.ACGSpectrumGammaPeak(j) = ...
            peakSpectrum(paramT.gammaRange(1),paramT.gammaRange(2),...
                        thetaModNonStim.ACGSpectrumFreq(j,:),...
                        thetaModNonStim.ACGSpectrum(j,:));

        % calculate the peak within high gamma frequencies of ACG spectrum 
        thetaModNonStim.ACGSpectrumHighGammaPeak(j) = ...
            peakSpectrum(paramT.highgammaRange(1),paramT.highgammaRange(2),...
                        thetaModNonStim.ACGSpectrumFreq(j,:),...
                        thetaModNonStim.ACGSpectrum(j,:));


        %% calculate peaks of the filtered ACG spectrum
        lenGaussKernelSpect = length(paramT.gaussFiltSpect);
        filteredACGSpectrum = conv(thetaModNonStim.ACGSpectrum(j,:),...
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
        thetaModNonStim.filteredACGSpectrum(j,:) = ...
                filteredACGSpectrum;

        % calculate the peak within delta frequencies of the filtered
        % ACG spectrum
        thetaModNonStim.filACGSpectrumDeltaPeak(j) = ...
            peakSpectrum(paramT.deltaRange(1),paramT.deltaRange(2),...
                        thetaModNonStim.ACGSpectrumFreq(j,:),...
                        filteredACGSpectrum);

        % calculate the peak within theta frequencies of the filtered
        % ACG spectrum
        thetaModNonStim.filACGSpectrumThetaPeak(j) = ...
            peakSpectrum(paramT.thetaRange(1),paramT.thetaRange(2),...
                        thetaModNonStim.ACGSpectrumFreq(j,:),...
                        filteredACGSpectrum);

        % calculate the peak within beta frequencies of the filtered
        % ACG spectrum
        thetaModNonStim.filACGSpectrumBetaPeak(j) = ...
            peakSpectrum(paramT.betaRange(1),paramT.betaRange(2),...
                        thetaModNonStim.ACGSpectrumFreq(j,:),...
                        filteredACGSpectrum);

        % calculate the peak within gamma frequencies of the filtered
        % ACG spectrum
        thetaModNonStim.filACGSpectrumGammaPeak(j) = ...
            peakSpectrum(paramT.gammaRange(1),paramT.gammaRange(2),...
                        thetaModNonStim.ACGSpectrumFreq(j,:),...
                        filteredACGSpectrum);

        % calculate the peak within high gamma frequencies of the filtered
        % ACG spectrum
        thetaModNonStim.filACGSpectrumHighGammaPeak(j) = ...
            peakSpectrum(paramT.highgammaRange(1),paramT.highgammaRange(2),...
                        thetaModNonStim.ACGSpectrumFreq(j,:),...
                        filteredACGSpectrum);

        filteredACG = conv(CCGNonStim.ccgVal(:,j,j),paramT.gaussFilt);
        lenGaussKernel = length(paramT.gaussFilt);
        if(mod(lenGaussKernel,2) == 0)
            filteredACG = ...
                filteredACG(floor(lenGaussKernel/2)+1:...
                    (end-floor(lenGaussKernel/2))+1); 
        else
            filteredACG = ...
                filteredACG(floor(lenGaussKernel/2)+1:...
                    (end-floor(lenGaussKernel/2))); 
        end
        [thetaModNonStim.trough(j),thetaModNonStim.peak(j),...
            thetaModNonStim.thetaModInd(j)] = ...
            thetaMod2(filteredACG,CCGNonStim.ccgT);
        thetaModNonStim.filteredACG(j,:) = filteredACG;

        [thetaModNonStim.trough3(j),thetaModNonStim.peak3(j),...
            thetaModNonStim.troughT3(j),thetaModNonStim.peakT3(j),...
            thetaModNonStim.thetaModInd3(j)] ...
            = thetaMod3(filteredACG,CCGNonStim.ccgT);
    end    
    
    fullPath = [path fileNameThetaModCtrl];
    save(fullPath, 'thetaModNonStim','paramT');
                      
    clear mydata;

end

function peakFreq = peakSpectrum(startFreq,stopFreq,spectrumFreq,...
                                spectrum)
    indFreq = find(spectrumFreq >= startFreq...
                & spectrumFreq <= stopFreq);
    [~,indPeak] = max(spectrum(indFreq));
    peakFreq = spectrumFreq(indFreq(indPeak(1)));
end
