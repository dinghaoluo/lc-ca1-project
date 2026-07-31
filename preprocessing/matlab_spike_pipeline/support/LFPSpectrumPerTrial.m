function LFPSpectrumPerTrial(path,fileName,figureState)
% calculate the delta oscillation frequency over left and right trials
% note:
% 1. using trials.eeg to compute the LFP spectrum
% 2. the spectrum of eeg and unit firing rate profile are calculated for
% individual trials before averaging
% path:         the path of the recording file
% fileName:     name of the recording file
% timeBin:      2SD of the Gaussian filter used to obtain the firing rate
%               profile (in second), default value is 0.005 s
% figureState:  0: figure off
%               1: plot the LFP spectrum for individual eeg trace
%               2: plot the LFP spectrum for individual eeg trace and
%               individual trial
%
% Example:
% LFPSpectrumPerTial
% ('../','i01_maze06_MS.005_DataStructure_mazeSection13_TrialType1_whlDirCW',1)

    %%%%%%%%% check arguments
    if nargin<2
        disp('At least two arguments are needed for this function.');
        return;
    elseif nargin == 2
        figureState = 0;
    elseif nargin > 3
        disp('Too many input arguments.');
        return;
    end

    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(isempty(indexFileName))
        fileNameOscFreq = [fileName '_LFPSpectrum.mat'];
        fileName = [fileName '.mat'];
    else
        fileNameOscFreq = [fileName(1:indexFileName(end)-1)...
                            '_LFPSpectrum.mat'];
    end

    LFPSpectrumStruct = [];
    thetaPAmp = [];
    thetaTAmp = [];

    fullPath = [path fileName];
    if(exist(fullPath) == 0)
        disp('File does not exist.');
        return;
    end
    load(fullPath,'trials');

    %%%%%%%%% initialize constants
    GlobalConst;
    fileNameInfo = [fileName(1:end-4) '_Info.mat'];

    fullPath = [path fileNameInfo];
    if(~exist(fullPath))
        BasicInfo(path,fileName);
    end
    load(fullPath);
    lapList = 1:beh.numTrials;
    
    %%%%%%%%% param structure
    paramOscFreq = struct(...
        'thetaBand', [lowerCrossCorrTheta,upperCrossCorrTheta],...
        'gammaBand', [lowerCrossCorrGamma upperCrossCorrGamma]);

    %%%%%%%%% prepare figure
    % prepare figure
    if(figureState == 1)
        [figNew,pos] = CreateFig();
        set(0,'Units','pixels')
        figTitle = 'LFP spectrum';
        set(figure(figNew),'OuterPosition',pos,'Name',figTitle);
    end
    
    fullPath = [path fileNameOscFreq];
    
    %%%%%%%%% calculate the variance of theta peak amplitude
    disp('Calculate peak theta amplitude:')
    thetaPeakAmpArr = getRecField(trials,'thetaPeak_tAmpl',lapList);
    thetaTroughAmpArr = getRecField(trials,'thetaTrough_tAmpl',lapList);
    speed = getRecField(trials,'speed',lapList);
    indSpeed = [];
    numSample = zeros(1,length(lapList));
    for i = 1:length(lapList)
        indSpeed{i} = find(speed{lapList(i)} > minSpeed);
        numSample(i) = length(indSpeed{lapList(i)});
    end
    maxNumSample = max(numSample);
   
    [thetaPAmp.MeanArr, thetaPAmp.StdArr, ...
        thetaPAmp.Mean, thetaPAmp.Std] = ...
            thetaPeakAmp(thetaPeakAmpArr,indSpeed,lapList);       
        % average peak theta amplitude 
    [thetaTAmp.MeanArr, thetaTAmp.StdArr, ...
        thetaTAmp.Mean, thetaTAmp.Std] = ...
            thetaPeakAmp(thetaTroughAmpArr,indSpeed,lapList);       
        % average trough theta amplitude
    save(fullPath,'paramOscFreq','thetaPAmp', 'thetaTAmp','-v7.3');
    
    %%%%%%%%% caculate LFP spectrum
    disp('Calculate LFP spectrum:');
    eeg = getRecField(trials,'eeg',lapList); 
        % extract eeg from the data structure
   
    RMSthetaFreqFiltRange = [6 10];
    if(~isempty(lapList))
        % calcualte both, whitened and non-whitened LFP power spectra
        for whiteSpect = 1
            switch whiteSpect
                case 0; pOpt = 0;
                        LFPSpectrumStruct = LFPSpectrumPar(eeg,lapList,...
                        paramOscFreq.thetaBand,paramOscFreq.gammaBand,...
                        maxNumSample,sampleFq,pOpt,indSpeed,figureState);
                        save(fullPath,'LFPSpectrumStruct','-append');
                case 1; pOptArr = zeros(1,length(lapList));
                        for i = 1:length(lapList)
                            [sbc, fpe, pOptArr(i)] = arest(eeg{i}, 1, 100);  
                        end
                        pOpt = floor(mean(pOptArr));
                        LFPSpectrumStructW = LFPSpectrumPar(eeg,lapList,...
                        paramOscFreq.thetaBand,paramOscFreq.gammaBand,...
                        maxNumSample,sampleFq,pOpt,indSpeed,figureState);
                        save(fullPath,'LFPSpectrumStructW','-append');
                    % using the first trial of the LFP data to estimate 
                    % the order of AR model
            end
        end
        % calculate theta RMS (as an estimate of theta amplitude that is
        % independent of peak/valley detection or phas estimate
        thetaAmpRMSmV = getThetaRms(eeg, sampleFq, ...
            RMSthetaFreqFiltRange, lapList, indSpeed);
        save(fullPath,'thetaAmpRMSmV', '-append');
    end
    
    if(figureState == 1)
        subplot(2,1,1)
        h = plot(LFPSpectrumStruct.freqLFP(LFPSpectrumStruct.indexFreqTheta)',...
            LFPSpectrumStruct.LFPSpectrum(:,LFPSpectrumStruct.indexFreqTheta'));
        set(h,'LineWidth',2);
        subplot(2,1,2)
         h = plot(LFPSpectrumStruct.freqLFP(LFPSpectrumStruct.indexFreqTheta)',...
            LFPSpectrumStruct.LFPSpectrum(:,LFPSpectrumStruct.indexFreqTheta'));
        set(h,'LineWidth',2);
    end    
end
