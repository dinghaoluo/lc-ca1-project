function [filteredSpikeArray,filteredSpikeArrayNormT] = ...
    ConvSpikeTrainTimePar_smTr2P(path,fileName,timeBin,onlyRun,saveFlag)
% Smoothed the spike trains for individual neurons over individual trials
% (over distance instead of time)
% -- Parallelized version
% path:         the path of the recording file
% fileName:     name of the recording file
% timeBin:      2SD of the Gaussian filter used to obtain the firing rate
%               profile (in second), default value is 200ms
% onlyRun:      1: only consider the time period when the animal is running 
% saveFlag:     0: do not save the result to a file 
%               otherwise: save the result to a file (default)
%
% Return:
% filteredSpikeArrayGo: a structure with N cells, where N = number of trials.
%                       Each cell i contains an array of size numNeurons x 
%                       Nsamples(i), with each row the smoothed firing 
%                       profile of each neuron   
%
% Example:
% [filteredSpikeArrayGo] = ConvSpikeTrainTimePar_smTr('./',
%           'A111-20150301-01_DataStructure_mazeSection1_TrialType1',0.2,1,1)

    %%%%%%%%% check arguments
    if nargin<2
        disp('At least two arguments are needed for this function.');
        return;
    elseif nargin == 2
        timeBin = 0.2; % s
        onlyRun = 1;
        saveFlag = 1;
    elseif nargin == 3
        onlyRun = 1;
        saveFlag = 1;
    elseif nargin == 4
        saveFlag = 1;
    elseif nargin > 5
        disp('Too many input arguments.');
        return;
    end
    
    GlobalConst2P;

    
    %%%%%%%%% load recording file
    indexFileName = strfind(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    timebinStr = num2str(timeBin*1000);
    ind = strfind(timebinStr,'.');
    if(~isempty(ind))
        timebinStr(ind) = 'p';
    end
    fileNameConv = [fileName '_convSpikesTime' timebinStr ...
                    'ms_Run' num2str(onlyRun) '.mat'];
    fileNameExt = [fileName '_ext.mat'];
    fileName = [fileName '.mat'];
    
    fullPath = [path fileName];
    if(exist(fullPath) == 0)
        disp('File does not exist.');
        return;
    end
    load(fullPath);
    
    if(onlyRun == 1)
        fullPath = [path fileNameExt];
        if(exist(fullPath) == 0)
            disp(['Extended file does not exist'....
                'Please run function SpikeDuringRun first']);
            return;
        end
        load(fullPath);
    end
    
    %%%%%%%%% initialize constants
    
    std = floor(timeBin*sampleFq); 
    % keep it the same as in the function ConvSpikeTrain_Aligned
    
    % generate Gaussian kernel
    paramC.gaussFilt = gaussFilter2P(12*std, std);
    lenGaussKernel = length(paramC.gaussFilt);
    normFactor = sum(paramC.gaussFilt);
    paramC.gaussFilt = paramC.gaussFilt./normFactor;
    
    filteredSpikeArray = [];
    filteredSpikeArrayNormTNormAmp = [];
    dFFArray = [];
    dFFArrayNormAmp = [];
    
    numTrials = length(trials);
    trialLenArr = [];
    for i = 1:numTrials
        trialLenArr = [trialLenArr trials{i}.Nsamples];
    end
    trialLen = max(trialLenArr);
    trialLen = min([trialLenT*sampleFq,trialLen]);
    numNeurons = size(trials{1}.spikes,2);
    
    tStart = tic;
    for i = 1:length(trials)
        % structure modified by Yingxue on 11/14/2021
        if(onlyRun == 0)
            if(isempty(trials{i}))
                filteredSpikeArray{i} = [];
                filteredSpikeArrayNormT{i} = [];
                filteredSpikeArrayNormTNormAmp{i} = [];
                continue;
            end
            spikes = trials{i}.spikesSM;
            dFF = trials{i}.dFFGF; % changed on 9/6/2023 from dFFSM to dFFGF
        else
            if(isempty(trialsExt{i}))
                filteredSpikeArray{i} = [];
                filteredSpikeArrayNormT{i} = [];
                filteredSpikeArrayNormTNormAmp{i} = [];
                continue;
            end
            spikes = trialsExt{i}.spikesSM;
            dFF = trialsExt{i}.dFFGF; % changed on 9/6/2023 from dFFSM to dFFGF
        end
        %%%%%%
        
%         disp(['Trial ' num2str(i)]);
        filteredSpikeArrayTmp = zeros(numNeurons,trialLen);
        filteredSpikeArrayNormAmpTmp = zeros(numNeurons,trialLen);
        dFFArrayTmp = zeros(numNeurons,trialLen);
        dFFArrayNormAmpTmp = zeros(numNeurons,trialLen);
        
        for j = 1:numNeurons
            filteredSpikeTmp = conv(spikes(:,j)',paramC.gaussFilt);
            filteredSpikeTmp = ...
                filteredSpikeTmp(floor(lenGaussKernel/2)+1:...
                (end-lenGaussKernel+floor(lenGaussKernel/2)+1)); 
                % cut the convolution result to be the same length as the original data
            if(length(filteredSpikeTmp) > trialLen)
                filteredSpikeTmp = filteredSpikeTmp(1:trialLen); 
            end
            filteredSpikeArrayTmp(j,1:length(filteredSpikeTmp)) = filteredSpikeTmp;  
            maxAmp = max(filteredSpikeArrayTmp(j,:));
            filteredSpikeArrayNormAmpTmp(j,:) = ...
                filteredSpikeArrayTmp(j,:)/maxAmp; 
            
            trLen = min(length(filteredSpikeTmp),trialLen);
            dFFArrayTmp(j,1:trLen) = dFF(1:trLen,j)';
            dFFArrayNormAmpTmp(j,1:trLen) = ...
                (dFF(1:trLen,j)'-min(dFFArrayTmp(j,:)))/...
                (max(dFFArrayTmp(j,:)) - min(dFFArrayTmp(j,:)));
        end
       
        filteredSpikeArray{i} = filteredSpikeArrayTmp*sampleFq; % added *sampleFq on 1/7/2022 ;
        filteredSpikeArrayNormAmp{i} = filteredSpikeArrayNormAmpTmp; 
        dFFArray{i} = dFFArrayTmp;
        dFFArrayNormAmp{i} = dFFArrayNormAmpTmp;
    end

    tLapse = toc(tStart);
    disp(['End of convolution calculation, total calculation time: ', ...
            num2str(tLapse)]);
    
    if(saveFlag ~= 0)
        fullPath = [path fileNameConv];
        save(fullPath, 'filteredSpikeArray',...
                       'filteredSpikeArrayNormAmp',...
                       'dFFArray',...
                       'dFFArrayNormAmp',...
                       'paramC','-v7.3'); 
    end

    clear mydata;
    clear all;

end
