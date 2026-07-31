function ThetaModulation(path,fileName,onlyRun,figState)
% calculate theta modulation of each neuron
%
% by Yingxue 08/25/2017

     
    %%%%%%%%% check arguments
    if nargin<2
        disp('At least three arguments are needed for this function.');
        return;
    elseif nargin == 2
        onlyRun = 1;
        figState = 0;
    elseif nargin == 3
        figState = 0;
    elseif nargin > 4
        disp('Too many input arguments');        
        return;
    end
        
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    fileNameThetaMod = [fileName '_ThetaMod_Run' num2str(onlyRun) '.mat'];
    fileNameCCG = [fileName '_CCG_Run' num2str(onlyRun) '.mat'];
    fileName = [fileName '.mat'];
    
    fullPath = [path fileName]; 
    if(exist(fullPath) == 0)
        disp('The file does not exist');
        return;
    end
    load(fullPath,'cluList');
    totClu = length(cluList.all);
    
    fullPath = [path fileNameCCG];
    if(exist(fullPath) == 0)
        disp(['Extended file does not exist. Please run function',...
              ' CalCCG first. onlyRun has to be equal to 1.']);
        return;
    end
    load(fullPath);
                   
    %%%%%%%%% initialize constants
    GlobalConst;
    fileNameInfo = [fileName(1:end-4) '_Info.mat'];
        
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo_smTr(path,fileName);
    end
    load(fullPath);
    mazeSess = beh.mazeSessAll;
    
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
                      'filteredACG',zeros(totClu,length(CCGSess{1}.ccgT)),... % smoothed ACG
                      'trough', zeros(1,totClu),... % ACG first trough (method 2)
                      'peak', zeros(1,totClu),... % ACG first peak (method 2)
                      'thetaModInd', zeros(1,totClu),... % theta modulation index (method 2)
                      ...
                      'trough3', zeros(1,totClu),... % ACG first trough (method 3)
                      'peak3', zeros(1,totClu),... % ACG first peak (method 3)
                      'troughT3', zeros(1,totClu),... % index of ACG first trough (method 3)
                      'peakT3', zeros(1,totClu),... % index of ACG first peak (method 3)
                      'thetaModInd3', zeros(1,totClu)); % theta modulation index (method 3)
                  
    % calculate theta modulation for each subsession
    disp('Calculate theta modulation for each subsession')
    thetaModSess = cell(1,length(mazeSess));
    for i = 1:length(mazeSess)
        disp(['Session ' num2str(i)]);
        thetaModSess{i} = thetaMod;
        for j = 1:totClu      
            [thetaModSess{i}.ACGSpectrum(j,:),...
                thetaModSess{i}.ACGSpectrumFreq(j,:),...
                thetaModSess{i}.thetaMod(j)] = ...
                thetaMod1(CCGSess{i}.ccgVal(:,j,j),CCGSess{i}.ccgT);
            
            %% calculate peaks of the ACG spectrum
            % calculate the peak within delta frequencies of ACG spectrum 
            thetaModSess{i}.ACGSpectrumDeltaPeak(j) = ...
                peakSpectrum(paramT.deltaRange(1),paramT.deltaRange(2),...
                            thetaModSess{i}.ACGSpectrumFreq(j,:),...
                            thetaModSess{i}.ACGSpectrum(j,:));
                        
            % calculate the peak within theta frequencies of ACG spectrum 
            thetaModSess{i}.ACGSpectrumThetaPeak(j) = ...
                peakSpectrum(paramT.thetaRange(1),paramT.thetaRange(2),...
                            thetaModSess{i}.ACGSpectrumFreq(j,:),...
                            thetaModSess{i}.ACGSpectrum(j,:));
                        
            % calculate the peak within beta frequencies of ACG spectrum 
            thetaModSess{i}.ACGSpectrumBetaPeak(j) = ...
                peakSpectrum(paramT.betaRange(1),paramT.betaRange(2),...
                            thetaModSess{i}.ACGSpectrumFreq(j,:),...
                            thetaModSess{i}.ACGSpectrum(j,:));
                        
            % calculate the peak within gamma frequencies of ACG spectrum 
            thetaModSess{i}.ACGSpectrumGammaPeak(j) = ...
                peakSpectrum(paramT.gammaRange(1),paramT.gammaRange(2),...
                            thetaModSess{i}.ACGSpectrumFreq(j,:),...
                            thetaModSess{i}.ACGSpectrum(j,:));
                        
            % calculate the peak within high gamma frequencies of ACG spectrum 
            thetaModSess{i}.ACGSpectrumHighGammaPeak(j) = ...
                peakSpectrum(paramT.highgammaRange(1),paramT.highgammaRange(2),...
                            thetaModSess{i}.ACGSpectrumFreq(j,:),...
                            thetaModSess{i}.ACGSpectrum(j,:));
            
            
            %% calculate peaks of the filtered ACG spectrum
            
            filteredACGSpectrum = conv(thetaModSess{i}.ACGSpectrum(j,:),...
                                paramT.gaussFiltSpect);
            if(mod(lenGaussKernel,2) == 0)
                filteredACGSpectrum = ...
                    filteredACGSpectrum(floor(lenGaussKernelSpect/2)+1:...
                        (end-floor(lenGaussKernelSpect/2))+1); 
            else
                filteredACGSpectrum = ...
                    filteredACGSpectrum(floor(lenGaussKernelSpect/2)+1:...
                        (end-floor(lenGaussKernelSpect/2))); 
            end
            thetaModSess{i}.filteredACGSpectrum(j,:) = ...
                    filteredACGSpectrum;
            
            % calculate the peak within delta frequencies of the filtered
            % ACG spectrum
            thetaModSess{i}.filACGSpectrumDeltaPeak(j) = ...
                peakSpectrum(paramT.deltaRange(1),paramT.deltaRange(2),...
                            thetaModSess{i}.ACGSpectrumFreq(j,:),...
                            filteredACGSpectrum);
                        
            % calculate the peak within theta frequencies of the filtered
            % ACG spectrum
            thetaModSess{i}.filACGSpectrumThetaPeak(j) = ...
                peakSpectrum(paramT.thetaRange(1),paramT.thetaRange(2),...
                            thetaModSess{i}.ACGSpectrumFreq(j,:),...
                            filteredACGSpectrum);
                        
            % calculate the peak within beta frequencies of the filtered
            % ACG spectrum
            thetaModSess{i}.filACGSpectrumBetaPeak(j) = ...
                peakSpectrum(paramT.betaRange(1),paramT.betaRange(2),...
                            thetaModSess{i}.ACGSpectrumFreq(j,:),...
                            filteredACGSpectrum);
                        
            % calculate the peak within gamma frequencies of the filtered
            % ACG spectrum
            thetaModSess{i}.filACGSpectrumGammaPeak(j) = ...
                peakSpectrum(paramT.gammaRange(1),paramT.gammaRange(2),...
                            thetaModSess{i}.ACGSpectrumFreq(j,:),...
                            filteredACGSpectrum);
                        
            % calculate the peak within high gamma frequencies of the filtered
            % ACG spectrum
            thetaModSess{i}.filACGSpectrumHighGammaPeak(j) = ...
                peakSpectrum(paramT.highgammaRange(1),paramT.highgammaRange(2),...
                            thetaModSess{i}.ACGSpectrumFreq(j,:),...
                            filteredACGSpectrum);
            
            filteredACG = conv(CCGSess{i}.ccgVal(:,j,j),paramT.gaussFilt);
            if(mod(lenGaussKernel,2) == 0)
                filteredACG = ...
                    filteredACG(floor(lenGaussKernel/2)+1:...
                        (end-floor(lenGaussKernel/2))+1); 
            else
                filteredACG = ...
                    filteredACG(floor(lenGaussKernel/2)+1:...
                        (end-floor(lenGaussKernel/2))); 
            end
            [thetaModSess{i}.trough(j),thetaModSess{i}.peak(j),...
                thetaModSess{i}.thetaModInd(j)] = ...
                thetaMod2(filteredACG,CCGSess{i}.ccgT);
            thetaModSess{i}.filteredACG(j,:) = filteredACG;
            
            %% added by Yingxue on 4/7/2020
            [thetaModSess{i}.trough3(j),thetaModSess{i}.peak3(j),...
                thetaModSess{i}.troughT3(j),thetaModSess{i}.peakT3(j),...
                thetaModSess{i}.thetaModInd3(j)] ...
                = thetaMod3(filteredACG,CCGSess{i}.ccgT);
            
        end
    end
    
    fullPath = [path fileNameThetaMod];
    save(fullPath, 'thetaModSess','paramT');
    
    if(figState ~= 0)
        numSess = length(mazeSess);
        meanThetaMod = zeros(1,numSess);
        stdThetaMod = zeros(1,numSess);
        exc = autoCorr.isPyrneuron;
        for i = 1:numSess
            thetaMod = thetaModSess{i}.thetaMod(exc);
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
