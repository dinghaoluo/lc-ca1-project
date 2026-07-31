function plotSequenceFRMap(path, fileName, spaceBin, onlyRun)
    
    fileNameInfo = [fileName '_Info.mat'];
    fileNameConv = [fileName '_convSpikesDist' num2str(spaceBin) ...
                    'mm_Run' num2str(onlyRun) '.mat'];
    fileNameSpInfo = [fileName '_SpInfo_Run' num2str(onlyRun) '.mat'];
    
%     neuronArr = [2,6,22,23,30,34,39,44,45,48,49,50,58,62,64,66,68];
        
    if(onlyRun == 1)
        fullPath = [path fileNameConv];
        if(exist(fullPath) == 0)
            disp(['Extended file does not exist. Please run function',...
                  ' ConvSpikeTrainDistPar first']);
            return;
        end
        load(fullPath);
    end
    
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo(path,fileName);
    end
    load(fullPath);
%     mazeSess = beh.mazeSessAll';
%     mazeType = beh.mazeSess;

    fullPath = [path fileNameSpInfo];
    if(exist(fullPath) == 0)
        disp(['Please run GetFRMapInfo function first to calculate spatial' ...
              'information']);
    end
    load(fullPath);
    
    mazeSess = beh.mazeSessAll';
    mazeType = beh.mazeSess;
    
    neuronArr1 = [];
    spInfoArr1 = [];
    neuronArr = [];
    for i = 1:length(beh.mazeSessAll)
        neuronArrTmp = find(spatialInfoSess{i}.spatialInfo > 0.2 &...
            spatialInfoSess{i}.sparsity > 1.15 & ...
            spatialInfoSess{i}.meanFR > 0.25); 
        neuronArr1 = [neuronArr1 neuronArrTmp];
        spInfoArr1 = [spInfoArr1 spatialInfoSess{i}.spatialInfo(neuronArrTmp)];
    end
    neuronArrU = unique(neuronArr1);
    for i = neuronArrU
        if(sum(neuronArr1 == i) > length(mazeSess)/2)
            neuronArr = [neuronArr i];
        else
            if(sum(spInfoArr1(neuronArr1 == i) > 0.5) >= 1)
                neuronArr = [neuronArr i];
            end
        end
    end
    
    GlobalConst; 
    
    filSpikeArrayNormTNormAmpPerNeuron = cell(1,3);
    for sess = mazeSess
        indTr = find(mazeType == sess);
        indTr = indTr';
        filSpikeArrayNormTNormAmpPerNeuron{sess} = ...
            zeros(length(neuronArr),size(filteredSpikeArrayNormTNormAmp{1},2));
        for i = 1:length(neuronArr)
            for tr = indTr
                filSpikeArrayNormTNormAmpPerNeuron{sess}(i,:) = ...
                    filSpikeArrayNormTNormAmpPerNeuron{sess}(i,:) + ...
                    filteredSpikeArrayNormTNormAmp{tr}(neuronArr(i),:);
            end
            filSpikeArrayNormTNormAmpPerNeuron{sess}(i,:) = ...
                filSpikeArrayNormTNormAmpPerNeuron{sess}(i,:)/...
                max(filSpikeArrayNormTNormAmpPerNeuron{sess}(i,:));
        end
        
        if(sess == 2)
            [maxFiltered,indMaxFiltered] = ...
                max(filSpikeArrayNormTNormAmpPerNeuron{sess},[],2);
            [~,indSorted] = sort(indMaxFiltered);
        end
    end
    
    for sess = mazeSess
        filSpikeArrayNormTNormAmpPerNeuronSorted{sess} = ...
            filSpikeArrayNormTNormAmpPerNeuron{sess}(indSorted,:);
        figure
        imagesc(paramC.spaceSteps{1}, 1:length(neuronArr),...
            filSpikeArrayNormTNormAmpPerNeuronSorted{sess},[0 1]); 
    end
    
end
