function ThetaModulation_CtrlOnly(path,fileName,onlyRun,mazeSess)
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
    elseif nargin == 3
        mazeSess = 1;
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
    fileNameThetaModCtrl = [fileName '_ThetaMod_Ctrl_Run' num2str(onlyRun) '_mazeSess' num2str(mazeSess) '.mat'];
    fileNameCCG = [fileName '_CCG_Ctrl_Run' num2str(onlyRun) '_mazeSess' num2str(mazeSess) '.mat'];
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
              ' CalCCG_GoodTr first. onlyRun has to be equal to 1.']);
        return;
    end
    load(fullPath);
    
    fullPath = [path fileNameThetaMod];
    if(exist(fullPath) == 0)
        disp(['_ThetaMod file does not exist. Please run function',...
              ' ThetaModulation first. ']);
        return;
    end
    load(fullPath,'paramT');
                   
    %%%%%%%%% initialize constants
    GlobalConst;
    fileNameInfo = [fileName(1:end-4) '_Info.mat'];
        
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo_smTr(path,fileName);
    end
    load(fullPath);
   
    thetaModSessCtrl = struct('ACGSpectrum',[],... % ACG spectrum
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
                      'filteredACG',zeros(totClu,length(CCGSessCtrl.ccgT)),... % smoothed ACG
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
    disp('Calculate theta modulation for ctrl trials')
    if(~isempty(CCGSessCtrl))
            thetaModSessCtrl = calThetaMod(totClu,thetaModSessCtrl,CCGSessCtrl,paramT);
    end
        
    fullPath = [path fileNameThetaModCtrl];
    save(fullPath, 'thetaModSessCtrl','paramT');
             
    clear mydata;

end

function peakFreq = peakSpectrum(startFreq,stopFreq,spectrumFreq,...
                                spectrum)
    indFreq = find(spectrumFreq >= startFreq...
                & spectrumFreq <= stopFreq);
    [~,indPeak] = max(spectrum(indFreq));
    peakFreq = spectrumFreq(indFreq(indPeak(1)));
end

function thetaModSess = calThetaMod(totClu,thetaModSess,CCGSess,paramT)
    if(isempty(CCGSess.ccgVal))
        return;
    end
    for j = 1:totClu  
        [thetaModSess.ACGSpectrum(j,:),...
            thetaModSess.ACGSpectrumFreq(j,:),...
            thetaModSess.thetaMod(j)] = ...
            thetaMod1(CCGSess.ccgVal(:,j,j),CCGSess.ccgT);

        %% calculate peaks of the ACG spectrum
        % calculate the peak within delta frequencies of ACG spectrum 
        thetaModSess.ACGSpectrumDeltaPeak(j) = ...
            peakSpectrum(paramT.deltaRange(1),paramT.deltaRange(2),...
                    thetaModSess.ACGSpectrumFreq(j,:),...
                    thetaModSess.ACGSpectrum(j,:));

        % calculate the peak within theta frequencies of ACG spectrum 
        thetaModSess.ACGSpectrumThetaPeak(j) = ...
            peakSpectrum(paramT.thetaRange(1),paramT.thetaRange(2),...
                    thetaModSess.ACGSpectrumFreq(j,:),...
                    thetaModSess.ACGSpectrum(j,:));
                
        % calculate the peak within beta frequencies of ACG spectrum 
        thetaModSess.ACGSpectrumBetaPeak(j) = ...
            peakSpectrum(paramT.betaRange(1),paramT.betaRange(2),...
                        thetaModSess.ACGSpectrumFreq(j,:),...
                        thetaModSess.ACGSpectrum(j,:));
                    
        % calculate the peak within gamma frequencies of ACG spectrum 
        thetaModSess.ACGSpectrumGammaPeak(j) = ...
            peakSpectrum(paramT.gammaRange(1),paramT.gammaRange(2),...
                        thetaModSess.ACGSpectrumFreq(j,:),...
                        thetaModSess.ACGSpectrum(j,:));    
                    
        % calculate the peak within high gamma frequencies of ACG spectrum 
        thetaModSess.ACGSpectrumHighGammaPeak(j) = ...
                peakSpectrum(paramT.highgammaRange(1),paramT.highgammaRange(2),...
                            thetaModSess.ACGSpectrumFreq(j,:),...
                            thetaModSess.ACGSpectrum(j,:));
                        
        %% calculate peaks of the filtered ACG spectrum         
        filteredACGSpectrum = conv(thetaModSess.ACGSpectrum(j,:),...
                            paramT.gaussFiltSpect);
        lenGaussKernelSpect = length(paramT.gaussFiltSpect);
        if(mod(lenGaussKernelSpect,2) == 0)
            filteredACGSpectrum = ...
                filteredACGSpectrum(floor(lenGaussKernelSpect/2)+1:...
                    (end-floor(lenGaussKernelSpect/2))+1); 
        else
            filteredACGSpectrum = ...
                filteredACGSpectrum(floor(lenGaussKernelSpect/2)+1:...
                    (end-floor(lenGaussKernelSpect/2))); 
        end
        thetaModSess.filteredACGSpectrum(j,:) = ...
                filteredACGSpectrum;             
            
        % calculate the peak within delta frequencies of the filtered
        % ACG spectrum
        thetaModSess.filACGSpectrumDeltaPeak(j) = ...
            peakSpectrum(paramT.deltaRange(1),paramT.deltaRange(2),...
                        thetaModSess.ACGSpectrumFreq(j,:),...
                        thetaModSess.filteredACGSpectrum(j,:));    
                    
        % calculate the peak within theta frequencies of the filtered
        % ACG spectrum
        thetaModSess.filACGSpectrumThetaPeak(j) = ...
            peakSpectrum(paramT.thetaRange(1),paramT.thetaRange(2),...
                        thetaModSess.ACGSpectrumFreq(j,:),...
                        thetaModSess.filteredACGSpectrum(j,:));     
                    
        % calculate the peak within beta frequencies of the filtered
        % ACG spectrum
        thetaModSess.filACGSpectrumBetaPeak(j) = ...
            peakSpectrum(paramT.betaRange(1),paramT.betaRange(2),...
                        thetaModSess.ACGSpectrumFreq(j,:),...
                        thetaModSess.filteredACGSpectrum(j,:));  
                    
        % calculate the peak within gamma frequencies of the filtered
        % ACG spectrum
        thetaModSess.filACGSpectrumGammaPeak(j) = ...
            peakSpectrum(paramT.gammaRange(1),paramT.gammaRange(2),...
                        thetaModSess.ACGSpectrumFreq(j,:),...
                        thetaModSess.filteredACGSpectrum(j,:));   
                    
        % calculate the peak within high gamma frequencies of the filtered
        % ACG spectrum
        thetaModSess.filACGSpectrumHighGammaPeak(j) = ...
            peakSpectrum(paramT.highgammaRange(1),paramT.highgammaRange(2),...
                        thetaModSess.ACGSpectrumFreq(j,:),...
                        thetaModSess.filteredACGSpectrum(j,:));   
                    
        filteredACG = conv(CCGSess.ccgVal(:,j,j),paramT.gaussFilt);
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
        [thetaModSess.trough(j),thetaModSess.peak(j),...
            thetaModSess.thetaModInd(j)] = ...
            thetaMod2(filteredACG,CCGSess.ccgT);
        thetaModSess.filteredACG(j,:) = filteredACG;        
        
        [thetaModSess.trough3(j),thetaModSess.peak3(j),...
            thetaModSess.troughT3(j),thetaModSess.peakT3(j),...
            thetaModSess.thetaModInd3(j)] ...
            = thetaMod3(thetaModSess.filteredACG(j,:),...
            CCGSess.ccgT);
    end
end
