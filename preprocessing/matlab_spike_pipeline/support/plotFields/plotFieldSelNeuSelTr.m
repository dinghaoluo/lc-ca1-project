function plotFieldSelNeuSelTr(path,fileName,spaceBin,methodTheta,onlyRun,sortMethod, indNeurons, cluBehNoGood, cluBehNoBad)
% plot field according to licking time within a trial
% sortMethod: 1: based on the first lick position
%             2: based on the running distance of the first segment
%             3: based on the running speed of the first segment
% e.g. plotFieldSelNeuSelTr('./','A012-20190221-01_DataStructure_mazeSection1_TrialType1',20,1,1,13,[16 27 30 38 65 85 99 108 134], [1 5], [0 2 3 4])
   
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    
    fileNameConv = [fileName '_convSpikesDist' num2str(spaceBin)...
                        'mm_Run' num2str(onlyRun) '.mat'];                               
    fileNameInfo = [fileName '_Info.mat'];
    fileNameRun = [fileName '_ext.mat'];
    fileNameFR = [fileName '_FR_Run' num2str(onlyRun) '.mat'];
    fileNameDepth = [fileName '_Depth.mat'];
    fileNameSpeed = [fileName '_runSpeed.mat'];
    fileNameRec = [fileName '.mat'];
    fileNameCluster = [fileName '_kmeanBeh.mat'];
    
    fullPath = [path fileNameConv];
    if(exist(fullPath) == 0)
        disp(['The firing profile file does not exist. Please call ',...
                'function "ConvSpikeTrainDistParVR" first.']);
        return;
    end
    load(fullPath,'filteredSpikeArrayNormT','paramC');
    
%     fullPath = [path fileNameSpInfo];
%     if(exist(fullPath) == 0)
%         disp(['The spatial information file does not exist. Please call ',...
%                 'function "GetFRMapInfo" first.']);
%         return;
%     end
        
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
    
    fullPath = [path fileNameRun];
    if(exist(fullPath) == 0)
        disp(['The spikes during run file does not exist. Please call ',...
                'function "SpikeDuringRun" first.']);
        return;
    end
    load(fullPath);
    
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
    
%     fullPath = [path fileNameSpeed];
%     if(exist(fullPath) == 0)
%         disp(['The run speed file does not exist. Please call ',...
%                 'function "RunSpeed" first.']);
%         return;
%     end
%     load(fullPath,'runSegments','speedSpectro');
    
    fullPath = [path fileNameRec];
    if(exist(fullPath) == 0)
        disp('The recording file does not exist.');
        disp(fullPath)
        return;
    end
    load(fullPath,'trials','cluList');
    
%     fullPath = [path fileNameCluster];
%     if(exist(fullPath) == 0)
%         disp('The behavior cluster file does not exist.');
%         disp(fullPath)
%         return;
%     end
%     load(fullPath);
    
    indGoodLap = beh.indGoodLap;
    %% sort trials 
    mazeSessTr = beh.mazeSess(indGoodLap);
    if(sortMethod == 0)
        % original field
%         for i = mazeSess'
        for i = 1:3
            ind = find(mazeSessTr == i);
            ind = beh.indGoodLap(ind);
            orderedParam{i} = zeros(1,length(ind));
            indGoodLapOrdered{i} = ind;
        end
        
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

    %         plotFRProfIndNeuronIndTrial(gca,filteredSpikeArrayNormT,i,...
    %                    orderedParam,indGoodLapOrdered,figTitle);
            plotSpikesIndNeuronIndTrial(gca,trialsExt,i,orderedParam,...
                indGoodLapOrdered,figTitle);

        end
    elseif(sortMethod == 13) % plotting the selected trials from a behavior cluster 
        indLapCluGood{1} = [];
        orderedParam{1} = [];
        for n = 1:length(cluBehNoGood)
            indClu = find(kmeanClu == cluBehNoGood(n));
            orderedParam{1} = [orderedParam{1} zeros(1,length(indClu))];
            indLapCluGood{1} = [indLapCluGood{1} indLaps(indClu)];
            indLapCluGood{1} = sort(indLapCluGood{1});
        end    
        
        count = 0;
        for i = indNeurons
            count = count + 1;

            if(mod(count-1,16) == 0)
                [figNew,pos] = CreateFig();
                set(0,'Units','pixels') 
                figTitle = 'All the trials';
                set(figure(figNew),'OuterPosition',...
                    [pos(1) pos(2) 987 808],'Name',figTitle)
            end

            subplot(4,4,mod(count-1,16)+1)

            figTitle = ['Neu ' num2str(i) '(' num2str(cluList.shank(i))...
                        ' ' num2str(cluList.localClu(i)) ') D' ...
                        num2str(depthNeu.relDepthNeuHDef(i))];

    %         plotFRProfIndNeuronIndTrial(gca,filteredSpikeArrayNormT,i,...
    %                    orderedParam,indGoodLapOrdered,figTitle);
            plotSpikesIndNeuronIndTrial(gca,trialsExt,i,orderedParam,...
                indLapCluGood,figTitle);
        end
        
        indLapCluBad{1} = [];
        orderedParam{1} = [];
        for n = 1:length(cluBehNoBad)
            indClu = find(kmeanClu == cluBehNoBad(n));
            orderedParam{1} = [orderedParam{1} zeros(1,length(indClu))];
            indLapCluBad{1} = [indLapCluBad{1} indLaps(indClu)];
        end    
        indLapCluBad{1} = sort(indLapCluBad{1});
        
        count = 0;
        for i = indNeurons
            count = count + 1;

            if(mod(count-1,16) == 0)
                [figNew,pos] = CreateFig();
                set(0,'Units','pixels') 
                figTitle = 'All the trials';
                set(figure(figNew),'OuterPosition',...
                    [pos(1) pos(2) 987 808],'Name',figTitle)
            end

            subplot(4,4,mod(count-1,16)+1)

            figTitle = ['Neu ' num2str(i) '(' num2str(cluList.shank(i))...
                        ' ' num2str(cluList.localClu(i)) ') D' ...
                        num2str(depthNeu.relDepthNeuHDef(i))];

    %         plotFRProfIndNeuronIndTrial(gca,filteredSpikeArrayNormT,i,...
    %                    orderedParam,indGoodLapOrdered,figTitle);
            plotSpikesIndNeuronIndTrial(gca,trialsExt,i,orderedParam,...
                indLapCluBad,figTitle);
        end
    end
end

function plotSpikesIndNeuronIndTrial(handle,spikes,neuronNo,paramSess,indLapSess,figTitle)
    numTrTot = 0;
    numTrSess = [];
    paramArr = [];
    for sess = 1:length(paramSess)
        numTrTot = numTrTot + length(paramSess{sess});
        numTrSess = [numTrSess numTrTot];
        paramArr = [paramArr paramSess{sess}];
    end
    
    numTr = 0;
    hold on;
    for sess = 1:length(paramSess)
        for i = 1:length(indLapSess{sess})
            numTr = numTr+1;
            h = plot(spikes{indLapSess{sess}(i)}.spikesMM{neuronNo},...
                numTr*ones(1,length(spikes{indLapSess{sess}(i)}.spikesMM{neuronNo})),...
                'k.');
            set(h,'MarkerSize',5.5);
        end
    end
    for i = 1:length(numTrSess)-1
        h = plot([0 1800],numTrSess(i)*ones(1,2),'r');
        set(h,'LineWidth',1);
    end
    set(gca,'XLim',[0 1800],'YLim',[0 numTrTot],'Ydir','reverse');
    xlabel('Dist (mm)');
    ylabel('Trial No.');
    title(figTitle);
end

function plotFRProfIndNeuronIndTrial...
            (handle,filteredSpikeArrayNormT,neuronNo,paramSess,indLapSess,figTitle)
    numSamples = size(filteredSpikeArrayNormT{1},2);
    numTrTot = 0;
    numTrSess = [];
    paramArr = [];
    for sess = 1:length(paramSess)
        numTrTot = numTrTot + length(paramSess{sess});
        numTrSess = [numTrSess numTrTot];
        paramArr = [paramArr paramSess{sess}];
    end
    FRProfilePerTrial = zeros(numTrTot,numSamples);
    
    numTr = 0;
    
    for sess = 1:length(paramSess)
        for i = 1:length(indLapSess{sess})
            numTr = numTr+1;
            FRProfilePerTrial(numTr,:) = filteredSpikeArrayNormT{indLapSess{sess}(i)}(neuronNo,:); 
            %./max(filteredSpikeArrayNormT{indLaps(i)}(neuronNo,:));
        end
    end  
    FRProfilePerTrial = FRProfilePerTrial/max(FRProfilePerTrial(:));
    h = imagesc(0:numSamples-1,1:numTr,FRProfilePerTrial);
    hold on
    for i = 1:length(numTrSess)-1
        h = plot([0 numSamples],numTrSess(i)*ones(1,2),'r');
        set(h,'LineWidth',1);
    end
    for i = 1:length(paramArr)
        if(~isnan(paramArr(i)))
            h = plot([paramArr(i) paramArr(i)],[i-1, i],'r');
            set(h,'LineWidth',1);
        end
    end
            
    set(gca,'FontSize',8.0,'Box','on','XLim',[0 numSamples-1],'YLim',[0 numTr]);
    xlabel('Dist (mm)');
    ylabel('Trial No.');
    title(figTitle);
end