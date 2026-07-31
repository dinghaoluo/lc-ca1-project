function [filteredSpikeArray,filteredSpikeArrayNormT] = ...
    ConvSpikeTrainDistPar_smTr(path,fileName,spaceBin,onlyRun,saveFlag)
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
    
    GlobalConst;
    
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
    
    if(onlyRun == 1) % this line is added by Yingxue on 11/14/2021
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
    paramC.gaussFilt = gaussFilter(12*std, std);
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
        % structure modified by Yingxue on 11/14/2021
        if(onlyRun == 0)
            if(isempty(trials{i}))
                filteredSpikeArray{i} = [];
                filteredSpikeArrayNormT{i} = [];
                filteredSpikeArrayNormTNormAmp{i} = [];
                continue;
            end
            Spikes = trials{i}.spikes;
            xMM = trials{i}.xMM;
        else
            if(isempty(trialsExt{i}))
                filteredSpikeArray{i} = [];
                filteredSpikeArrayNormT{i} = [];
                filteredSpikeArrayNormTNormAmp{i} = [];
                continue;
            end
            Spikes = trialsExt{i}.spikes;
            %%%% corrected a bug on 4/2/2019, should only consider the
            %%%% distance that is over the speed threshold when
            %%%% calculating the timePerBin
            xMM = trialsExt{i}.xMM;
        end
        %
        
        indTrack = find(tracks == beh.trackLen(i));
        [spikesInd,spikesDistTmp,timePerBinTmp] = ...
            spikeTime2Dist(Spikes,trials{i}.xMM, xMM,...
                    paramC.spaceSteps{indTrack});
        
        nsamples = length(paramC.spaceSteps{indTrack});
        spikesDist.distPerSpike{i} = spikesDistTmp;
        spikesDist.timePerBin{i} = timePerBinTmp/sampleFq;
        spikesDist.avgTimePerBin{i} = sum(timePerBinTmp)/sampleFq/beh.trackLen(i); % added time estimation per space bin on 1/3/2022 by Yingxue
        
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
            if(~isempty(spikesDistTmp{j}))
                for m = 1:length(spikesDistTmp{j})
                    ind = spikesInd{j}(m);
                    spikeArray(ind) = spikeArray(ind) + 1; 
                    spikeArrayNorm(ind) = spikeArrayNorm(ind) ...
                        + 1/spikesDist.timePerBin{i}(ind); % changed from spikesDist.timePerBin{i}(ind) to timePerBinTmp(ind) on 11/7/2021 by Yingxue;
                                                           % changed back to spikesDist.timePerBin{i}(ind) on 1/3/2022 by Yingxue                             
                end
                spikeArray = spikeArray/spikesDist.avgTimePerBin{i}; % added on 1/3/2022 by Yingxue
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
function [spikesInd,spikesDist,timePerBin] = ...
                            spikeTime2Dist(spikes,dist,distRun, spaceSteps)
                        
    numBins = length(spaceSteps);
    step = spaceSteps(2) - spaceSteps(1);
    timePerBin = zeros(1,numBins);
    for i = 1:numBins
        ind = find(distRun >= spaceSteps(i)-step/2 & distRun < spaceSteps(i)+step/2);
        if(~isempty(ind))
            timePerBin(i) = length(ind);
        else
            timePerBin(i) = 1;
        end
    end

    numNeurons = length(spikes);
    spikesInd = cell(1,numNeurons);
    spikesDist = cell(1,numNeurons);
    for i = 1:numNeurons % number of neurons
        for j = 1:length(spikes{i}) % number of spikes per neuron per trial
            curSpike = dist(spikes{i}(j));
            tmp = spaceSteps - curSpike;
            %%%%% added on 4/8/2019 by Yingxue
            %%% remove all the spikes that occurred after the trackLen
            if(curSpike > spaceSteps(end) + step/2)
                continue;
            end
            %%%%%
            ind = find(tmp >= 0,1);
            if(tmp(ind) == 0 | ind == 1)
                spikesInd{i} = [spikesInd{i} ind];
                spikesDist{i} = [spikesDist{i} spaceSteps(ind)];
            else
                if(tmp(ind) < abs(tmp(ind-1)))
                    spikesInd{i} = [spikesInd{i} ind];
                    spikesDist{i} = [spikesDist{i} spaceSteps(ind)];
                else
                    spikesInd{i} = [spikesInd{i} ind-1];
                    spikesDist{i} = [spikesDist{i} spaceSteps(ind-1)];
                end
            end
        end
    end
end