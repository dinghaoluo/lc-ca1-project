function pcaConvSpikesDist(path, fileName, spaceBin, onlyRun)
% calculate pca of the convolution of spikes over distance
    
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    
    fileNameConv = [fileName '_convSpikesDist' num2str(spaceBin)...
                        'mm_Run' num2str(onlyRun) '.mat'];                    
    fileNameInfo = [fileName '_Info.mat'];
    fileNameFR = [fileName '_FR_Run' num2str(onlyRun) '.mat'];
    fileNameRec = [fileName '.mat'];
    
    fullPath = [path fileNameConv];
    if(exist(fullPath) == 0)
        disp(['The firing profile file does not exist. Please call ',...
                'function "ConvSpikeTrainDistParVR" first.']);
        return;
    end
    load(fullPath,'filteredSpikeArrayNormT','paramC');
    
    fullPath = [path fileNameInfo];
    if(exist(fullPath,'file') == 0)
        BasicInfo(path,fileName);
    end
    load(fullPath);
    mazeSess = beh.mazeSessAll;
    numSamples = zeros(1,length(unique(beh.trackLen)));
    for i = 1:length(numSamples)
        numSamples(i) = length(paramC.spaceSteps{i});
    end
    maxNumSamples = max(numSamples);
    
    fullPath = [path fileNameFR];
    if(exist(fullPath) == 0)
        disp(['The neuron depth file does not exist. Please call ',...
                'function "GetNeuRelativeDepth" first.']);
        return;
    end
    load(fullPath,'mFRStruct');
    
    fullPath = [path fileNameRec];
    if(exist(fullPath) == 0)
        disp('The recording file does not exist.');
        disp(fullPath)
        return;
    end
    load(fullPath,'trials','cluList');
    
    GlobalConst;
    
    pyrNeurons = find(rec.isIntern == 0 & ...
        mFRStruct.mFR > minFR);
    
    % order trials based on the first lick position
    mazeSessTr = beh.mazeSess(beh.indGoodLap);
    firstLickPos = zeros(1,length(beh.indGoodLap));
    for i = beh.indGoodLap
        lickDist = trials{i}.xMM(trials{i}.lickLfpInd);
        ind = find(lickDist > 300, 1, 'first'); 
        firstLickPos(i) = mean(lickDist(ind));
    end
    indGoodLapOrdered = [];
    for i = 1
        ind = find(mazeSessTr == i);
        [~, indGoodLapOrderedTmp] = sort(firstLickPos(ind));
        indGoodLapOrdered = [indGoodLapOrdered; ind(indGoodLapOrderedTmp)];
    end

    filteredSpikeArray = zeros(length(pyrNeurons),...
        length(indGoodLapOrdered)*size(filteredSpikeArrayNormT{1},2));
    for i = 1:length(pyrNeurons)
        filteredSpikeArrayNeu = zeros(length(indGoodLapOrdered),...
            size(filteredSpikeArrayNormT{1},2));
        for j = 1:length(indGoodLapOrdered)
            filteredSpikeArrayNeu(j,:) = ...
                filteredSpikeArrayNormT{indGoodLapOrdered(j)}(pyrNeurons(i),:);
        end
        filteredSpikeArray(i,:) = reshape(filteredSpikeArrayNeu',1,[]);
    end
            
%     for i = 1:length(beh.indGoodLap)
%         filteredSpikeArray(i,:) = ...
%             filteredSpikeArrayNormT{beh.indGoodLap(i)}(NeuronNo,:);
%     end
    filteredSpikeArray = ...
        bsxfun(@minus,filteredSpikeArray,mean(filteredSpikeArray));
    
    [coeff,score,latent,tsquared,explained,mu] = ...
        pca(filteredSpikeArray);
    
    proj1 = filteredSpikeArray*coeff(:,1);
    [~,indProj1] = sort(proj1);
    figure;
    imagesc(filteredSpikeArray(indProj1,:))
    
    proj2 = filteredSpikeArray*coeff(:,2);
    [~,indProj2] = sort(proj2);
    figure
    imagesc(filteredSpikeArray(indProj2,:))
    
    proj3 = filteredSpikeArray*coeff(:,3);
    [~,indProj3] = sort(proj3);
    figure
    imagesc(filteredSpikeArray(indProj3,:))
    
    proj4 = filteredSpikeArray*coeff(:,4);
    [~,indProj4] = sort(proj4);
    figure
    imagesc(filteredSpikeArray(indProj4,:))
    
    z = linkage(projArr(:,1:15),'ward');
    c = cluster(z,'maxclust',10);
    dendrogram(z,0)
    
end
