function [filteredSpikeArray,filteredSpikeArrayNormT] = ...
    ConvSpikeTrainDistPar_smTr2P(path,fileName,spaceBin,onlyRun,saveFlag)
% Smoothed the spike trains for individual neurons over individual trials
% (over distance instead of time)
% -- Parallelized version
% path:         the path of the recording file
% fileName:     name of the recording file
% spaceBin:      2SD of the Gaussian filter used to obtain the firing rate
%               profile (in second), default value is 70 pixels
% sampFq:       sample frequency
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
% [filteredSpikeArrayGo] = ConvSpikeTrainDistPar_smTr('./',
%           'A111-20150301-01_DataStructure_mazeSection1_TrialType1',10,1,1)

    %%%%%%%%% check arguments
    if nargin<2
        disp('At least two arguments are needed for this function.');
        return;
    elseif nargin == 2
        spaceBin = 2; % cm
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
    fileNameConv = [fileName '_convSpikesDist' num2str(spaceBin) ...
                    'mm_Run' num2str(onlyRun) '.mat'];
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
    fileNameInfo = [fileName(1:end-4) '_Info.mat'];     
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo(path,fileName);
    end
    load(fullPath);
    
    std = floor(spaceBin/spaceMergeBin/2);
    
    % generate Gaussian kernel
    paramC.gaussFilt = gaussFilter2P(12*std, std);
    lenGaussKernel = length(paramC.gaussFilt);
    normFactor = sum(paramC.gaussFilt);
    paramC.gaussFilt = paramC.gaussFilt./normFactor;
    
    filteredSpikeArray = [];
    filteredSpikeArrayNormT = [];
    
    tracks = unique(beh.trackLen);
    for i = 1:length(tracks)
        if(spaceMergeBin ~= 0)
            paramC.spaceSteps{i} = [0:spaceMergeBin:tracks(i)];
        else
            paramC.spaceSteps{i} = [0:tracks(i)];
        end
        paramC.spaceBin = spaceBin;
    end
    
    spikesDist = [];
        
    tStart = tic;
    for i = 1:length(lapList)
        if(onlyRun == 0)
            if(isempty(trials{i}))
                filteredSpikeArray{i} = [];
                filteredSpikeArrayNormT{i} = [];
                filteredSpikeArrayNormTNormAmp{i} = [];
                continue;
            end
            Spikes = trials{i}.spikesSM;
            xMM = trials{i}.xMM;
        else
            if(isempty(trialsExt{i}))
                filteredSpikeArray{i} = [];
                filteredSpikeArrayNormT{i} = [];
                filteredSpikeArrayNormTNormAmp{i} = [];
                continue;
            end
            Spikes = trialsExt{i}.spikesSM;
            %%%% corrected a bug on 4/2/2019, should only consider the
            %%%% distance that is over the speed threshold when
            %%%% calculating the timePerBin
            xMM = trialsExt{i}.xMM;
        end
        
        indTrack = find(tracks == beh.trackLen(i));
        [spikesPerBinTmp,timePerBinTmp] = ...
            spikeTime2Dist(Spikes,xMM,...
                    paramC.spaceSteps{indTrack});
        
        nsamples = length(paramC.spaceSteps{indTrack});
        spikesDist.spikesPerBin{i} = spikesPerBinTmp;
        spikesDist.timePerBin{i} = timePerBinTmp/sampleFq;
        spikesDist.avgTimePerBin{i} = sum(timePerBinTmp)/sampleFq/(beh.trackLen(i)/spaceMergeBin); % added time estimation per space bin on 1/7/2022 by Yingxue
        
%         disp(['Trial ' num2str(i)]);
        if(i == 66)
            a = 1;
        end
        filteredSpikeArrayTmp = zeros(rec.numNeurons,nsamples);
        filteredSpikeArrayNormTTmp = zeros(rec.numNeurons,nsamples);
        filteredSpikeArrayNormTNormAmpTmp = zeros(rec.numNeurons,nsamples);
        for j = 1:rec.numNeurons
%             fprintf('%d  ', j);
            spikeArray = zeros(1,nsamples);
            spikeArrayNorm = zeros(1,nsamples);
            if(sum(spikesPerBinTmp(j,:)))
                spikeArray = spikesPerBinTmp(j,:)/spikesDist.avgTimePerBin{i}; % changed by Yingxue on 1/7/2022, added /spikesDist.avgTimePerBin{i}
                spikeArrayNorm = spikesPerBinTmp(j,:)./spikesDist.timePerBin{i}; % changed by Yingxue on 2/8/2022 from timePerBinTmp to spikesDist.timePerBin{i}
                spikeArrayNorm(isnan(spikeArrayNorm)) = 0;
                spikeArray1 = [spikeArray(nsamples-lenGaussKernel+1:nsamples)...
                        spikeArray spikeArray(1:lenGaussKernel)];
                spikeArrayNorm1 = [spikeArrayNorm(nsamples-lenGaussKernel+1:nsamples)...
                        spikeArrayNorm spikeArrayNorm(1:lenGaussKernel)];
                filteredSpikeTmp = conv(spikeArray1,paramC.gaussFilt);
                filteredSpikeArrayTmp(j,:) = ...
                    filteredSpikeTmp(floor(lenGaussKernel/2)+lenGaussKernel+1:...
                        (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1)); 
                    % cut the convolution result to be the same length 
                    % as the original data
                filteredSpikeTmp = conv(spikeArrayNorm1,paramC.gaussFilt);
                filteredSpikeArrayNormTTmp(j,:) = ...
                    filteredSpikeTmp(floor(lenGaussKernel/2)+lenGaussKernel+1:...
                        (end-2*lenGaussKernel+floor(lenGaussKernel/2)+1));
                filteredSpikeArrayNormTNormAmpTmp(j,:) = ...
                    filteredSpikeArrayNormTTmp(j,:)...
                    /max(filteredSpikeArrayNormTTmp(j,:));
            end
        end
        filteredSpikeArray{i} = filteredSpikeArrayTmp;
        filteredSpikeArrayNormT{i} = filteredSpikeArrayNormTTmp;
        filteredSpikeArrayNormTNormAmp{i} = filteredSpikeArrayNormTNormAmpTmp;
    end

    tLapse = toc(tStart);
    disp(['End of convolution calculation, total calculation time: ', ...
            num2str(tLapse)]);
    
    if(saveFlag ~= 0)
        fullPath = [path fileNameConv];
        save(fullPath, 'filteredSpikeArray','filteredSpikeArrayNormT',...
                       'filteredSpikeArrayNormTNormAmp',...
                       'spikesDist','paramC','-v7.3'); 
    end

    clear mydata;
    clear all;

end

%%%% corrected a bug on 4/2/2019, should only consider the
%%%% distance that is over the speed threshold when calculating the time
%%%% per bin
function [spikesPerBin,timePerBin] = ...
                            spikeTime2Dist(spikes,dist, spaceSteps)
                        
    numBins = length(spaceSteps);
    step = spaceSteps(2) - spaceSteps(1);

    timePerBin = hist(dist(dist<=spaceSteps(end)+step/2),spaceSteps);
%     for i = 1:numBins
%         ind = find(dist >= spaceSteps(i)-step/2 & dist < spaceSteps(i)+step/2);
%         if(~isempty(ind))
%             timePerBin(i) = length(ind);
%         else
%             timePerBin(i) = 1;
%         end
%     end

    numNeurons = size(spikes,2);
    spikesPerBin = zeros(numNeurons,numBins);
    accumTime = 1;
    for i = 1:numBins
        if(timePerBin(i) == 0)
            continue;
        end
        spikePerBinTmp = spikes(accumTime:accumTime+timePerBin(i)-1,:)';
        spikesPerBin(:,i) = sum(spikePerBinTmp,2);
        accumTime = accumTime + timePerBin(i);
    end
    
end