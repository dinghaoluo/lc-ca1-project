function plotFieldOrderedByLick(path,fileName,spaceBin,methodTheta,onlyRun)
% plot field according to licking time within a trial
% sortMethod: 1: based on the first lick position
%             2: based on the running distance of the first segment
%             3: based on the running speed of the first segment

    alignDist = 500; % align the first licks to 500 mm
    
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    
    fileNameConv = [fileName '_convSpikesDist' num2str(spaceBin)...
                        'mm_Run' num2str(onlyRun) '.mat'];                    
    fileNameSpInfo = [fileName '_SpInfo_Run' num2str(onlyRun) '.mat'];                    
    fileNameInfo = [fileName '_Info.mat'];
    fileNameFR = [fileName '_FR_Run' num2str(onlyRun) '.mat'];
    fileNameDepth = [fileName '_Depth.mat'];
    fileNameSpeed = [fileName '_runSpeed.mat'];
    fileNameRec = [fileName '.mat'];
    
    fullPath = [path fileNameConv];
    if(exist(fullPath) == 0)
        disp(['The firing profile file does not exist. Please call ',...
                'function "ConvSpikeTrainDistParVR" first.']);
        return;
    end
    load(fullPath,'filteredSpikeArrayNormT','paramC');
    
    fullPath = [path fileNameSpInfo];
    if(exist(fullPath) == 0)
        disp(['The spatial information file does not exist. Please call ',...
                'function "GetFRMapInfo" first.']);
        return;
    end
        
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
    
    fullPath = [path fileNameDepth];
    if(exist(fullPath) == 0)
        disp(['The mean firing rate file does not exist. Please call ',...
                'function "MeanFiringRateVR" first.']);
        return;
    end
    load(fullPath);
    
    fullPath = [path fileNameSpeed];
    if(exist(fullPath) == 0)
        disp(['The run speed file does not exist. Please call ',...
                'function "RunSpeed" first.']);
        return;
    end
    load(fullPath,'runSegments');
    
    fullPath = [path fileNameRec];
    if(exist(fullPath) == 0)
        disp('The recording file does not exist.');
        disp(fullPath)
        return;
    end
    load(fullPath,'trials','cluList');
    
    indNeurons = 1:rec.numNeurons;
    indGoodLap = beh.indGoodLap;
    %% sort trials based on the first lick position
    mazeSessTr = beh.mazeSess(indGoodLap);
    
    % find first lick position
    firstLickPos = zeros(1,length(indGoodLap));
    for i = indGoodLap
        lickDist = trials{i}.xMM(trials{i}.lickLfpInd);
        ind = find(lickDist > 300, 1, 'first'); 
        firstLickPos(i) = mean(lickDist(ind));
    end
        
%         ind = [1:119]';
%         [firstLickPosOrdered{1}, indGoodLapOrdered{1}] = ...
%                 sort(firstLickPos(ind));
%         indGoodLapOrdered{1} = ind(indGoodLapOrdered{1});

    % order trials by where the lick occurs
    for i = mazeSess'
        ind = find(mazeSessTr == i);
        [firstLickPosOrdered{i}, indGoodLapOrdered{i}] = ...
            sort(firstLickPos(ind));
        indGoodLapOrdered{i} = ind(indGoodLapOrdered{i});
    end
    
    filteredSpikeArrayByLick = ...
        alignFilteredSpArrayToLick(indGoodLapOrdered,firstLickPos,...
                                    filteredSpikeArrayNormT,alignDist);
    
    count = 0;
    for i = indNeurons
        count = count + 1;

        if(mod(count-1,16) == 0)
            [figNew,pos] = CreateFig();
            set(0,'Units','pixels') 
            figTitle = 'All the trials';
            set(figure(figNew),'OuterPosition',...
                [pos(1) pos(2) pos(3)*2 pos(4)*2.2],'Name',figTitle)
        end
        
        subplot(4,4,mod(count-1,16)+1)
        
        figTitle = ['Neu ' num2str(i) '(' num2str(cluList.shank(i))...
                    ' ' num2str(cluList.localClu(i)) ') D' ...
                    num2str(depthNeu.relDepthNeuHDef(i))];
                
        plotFRProfIndNeuronIndTrial(gca,filteredSpikeArrayByLick,i,...
                   indGoodLapOrdered,alignDist,figTitle);
    end 
end

function plotFRProfIndNeuronIndTrial...
            (handle,filteredSpikeArray,neuronNo,indLapSess,alignDist,figTitle)
    numSamples = size(filteredSpikeArray{1},2);
    numTrTot = 0;
    numTrSess = [];
    lickPosArr = [];
    for sess = 1:length(indLapSess)
        numTrTot = numTrTot + length(indLapSess{sess});
        numTrSess = [numTrSess numTrTot];
    end
    FRProfilePerTrial = zeros(numTrTot,numSamples);
    
    numTr = 0;
    
    for sess = 1:length(indLapSess)
        for i = 1:length(indLapSess{sess})
            numTr = numTr+1;
            FRProfilePerTrial(numTr,:) = filteredSpikeArray{neuronNo}(indLapSess{sess}(i),:); 
            %./max(filteredSpikeArray{indLaps(i)}(neuronNo,:));
        end
    end  
    FRProfilePerTrial = FRProfilePerTrial/max(FRProfilePerTrial(:));
    h = imagesc(0:numSamples-1,1:numTr,FRProfilePerTrial);
    hold on
    for i = 1:length(numTrSess)-1
        h = plot([0 numSamples],numTrSess(i)*ones(1,2),'r');
        set(h,'LineWidth',1);
    end
    h = plot([alignDist alignDist],[0, numTrTot],'r');
    set(h,'LineWidth',1);
     
    set(gca,'FontSize',8.0,'Box','on','XLim',[0 numSamples-1],'YLim',[0 numTr]);
    xlabel('Dist (mm)');
    ylabel('Trial No.');
    title(figTitle);
end