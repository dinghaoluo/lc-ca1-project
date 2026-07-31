function FieldDetection_CtrlOnly(path,fileName,onlyRun,spaceBin,figureState,mazeSess)
% Firing field detection based on correlation distance and spatial information
% path:         the path of the recording file
% fileName:     name of the recording file
% onlyRun:      1: only consider the time period when the animal is running 
%
% Example:
% FieldDetection_GoodTr('./','A111-20150301-01_DataStructure_mazeSection1_TrialType1',1,2)

    %%%%%%%%% check arguments
    if nargin<2
        disp('At least two arguments are needed for this function.');
        return;
    elseif nargin == 2
        figureState = 0;
        spaceBin = 20; % mm
        onlyRun = 1;
        mazeSess = 1;
    elseif nargin == 3
        figureState = 0;   
        spaceBin = 20; % mm
        mazeSess = 1;
    elseif nargin == 4
        figureState = 0; 
        mazeSess = 1;
    elseif nargin == 5
        mazeSess = 1;
    elseif nargin > 6
        disp('Too many input arguments.');
        return;
    end
    
    %%%%%%%%% initialize constants
    GlobalConst;
        
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    
    fileNameFW = [fileName '_FieldSpCorr_Ctrl_Run' num2str(onlyRun) ...
        '_mazeSess' num2str(mazeSess) '.mat'];
    
    fileNameConv = [fileName '_convSpikesDist' num2str(spaceBin)...
                        'mm_Run' num2str(onlyRun) '.mat'];
    fileNamePeakFR = [fileName '_PeakFRCtrl' num2str(spaceBin)...
                        'mm_Run' num2str(onlyRun) '_mazeSess' num2str(mazeSess) '.mat'];
    fileNameFR = [fileName '_FR_Ctrl_Run' num2str(onlyRun) ...
             '_mazeSess' num2str(mazeSess) '.mat'];
    fileNameSpInfo = [fileName '_SpInfo_Ctrl_Run' num2str(onlyRun)...
             '_mazeSess' num2str(mazeSess) '.mat'];
    fileNameCorr = [fileName '_meanSpikesCorrDist_Ctrl_Run' num2str(onlyRun)...
             '_mazeSess' num2str(mazeSess) '.mat'];
    fileNameInfo = [fileName '_Info.mat'];
    fileName = [fileName '.mat'];
     
    fullPath = [path fileName];
    if(exist(fullPath) == 0)
        disp('The recording file does not exist.');
        disp(fullPath)
        return;
    end
    load(fullPath,'cluList');
    
    fullPath = [path fileNameConv];
    if(exist(fullPath) == 0)
        disp(['The firing profile file does not exist. Please call ',...
                'function "ConvSpikeTrainDistPar" first.']);
        return;
    end
    load(fullPath,'filteredSpikeArrayNormT','paramC');
    
    fullPath = [path fileNamePeakFR];
    if(exist(fullPath) == 0)
        disp(['The peak firing rate file does not exist. Please call ',...
                'function "PeakFiringRate_smTrCtrlOnly" first.']);
        return;
    end
    load(fullPath,'pFRStructNormTSessCtrl');
    
    fullPath = [path fileNameFR];
    if(exist(fullPath) == 0)
        disp(['The mean firing rate file does not exist. Please call ',...
                'function "MeanFiringRateCtrlOnly" first.']);
        return;
    end
    load(fullPath,'mFRStructSessCtrl');
    
    fullPath = [path fileNameSpInfo];
    if(exist(fullPath) == 0)
        disp(['The spatial information file does not exist. Please call ',...
                'function "GetFRMapInfo_CtrlOnly" first.']);
        return;
    end
    load(fullPath,'spatialInfoSessCtrl');
    
    fullPath = [path fileNameCorr];
    if(exist(fullPath) == 0)
        disp(['The mean correlation in distance file does not exist. Please call ',...
                'function "meanNeuronSpikeCorrDist_CtrlOnly" first.']);
        return;
    end
    load(fullPath,'meanCorrDistCtrl');
    
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo_smTr(path,fileName);
    end
    load(fullPath);
    
    paramF.minNumTr = 15;
    paramF.percNumActiveTr = 0.5;
    paramF.minCorr = 0.15;
    paramF.minSpInfo = 0.25;
    paramF.minCorrHighSpInfo = 0.09;
    paramF.minHighSpInfo = 0.7;
    paramF.minInstFR = 0.15;
    paramF.lowThreFieldMeanInstFR = 0.08;
    paramF.maxFieldWidth = 1500;
    paramF.reboundCheckRegion = 300; % there is no clear peak within ReboundCheckRegion (mm) region before or after the field
    paramF.maxReboundMean = 0.2; % rebound mean inst FR should be smaller than maxReboundMean * meanInstFR
    paramF.reboundHeight = 0.4; % rebound height should be smaller than reboundHeight*peakFieldFR
    
    numSamples = zeros(1,length(unique(beh.trackLen)));
    for i = 1:length(numSamples)
        numSamples(i) = length(paramC.spaceSteps{i});
    end
    
    if(~isempty(meanCorrDistCtrl))
            fieldSpCorrSessCtrl= FieldSpInfoCorr(...
                pFRStructNormTSessCtrl,spatialInfoSessCtrl.adaptSpatialInfo,...
                meanCorrDistCtrl,paramF);
    end
   
    fullPath = [path fileNameFW];
    save(fullPath, 'fieldSpCorrSessCtrl');
    fieldSpCorrSessCtrl
    
    if(figureState == 2)
       %%% all the trials
       if(~isempty(fieldSpCorrSessCtrl))
       %% good trials
            indLaps = mFRStructSessCtrl.indLapList;
            count = 0;
            for i = 1:rec.numNeurons  
                if(mFRStructSessCtrl.mFR(i) > 0.05) %minFR
                   strNumField = '';
                   indSessBorder = [];
                   totTrialNo = 0;
                   fieldInfoTmp = getFieldInfoIndNeuron(i,fieldSpCorrSessCtrl); 
                    if(~isempty(fieldInfoTmp))
                        strNumField = [strNumField, ' ', num2str(size(fieldInfoTmp,1))];
                    else
                        strNumField = [strNumField, ' 0'];
                    end
                    % get the field information
                    count = count + 1;

                    if(mod(count-1,16) == 0)
                        [figNew,pos] = CreateFig();
                        set(0,'Units','pixels') 
                        figTitle = 'All the trials';
                        set(figure(figNew),'OuterPosition',...
                            [pos(1) pos(2)-500 pos(3)*2 pos(4)*2.2],'Name',figTitle)
                    end

                    subplot(4,4,mod(count-1,16)+1)

                    figTitle = ['Neu ' num2str(i) '(' num2str(cluList.shank(i))...
                        ' ' num2str(cluList.localClu(i)) ') Field' strNumField ' G Sp' ...
                        num2str(spatialInfoSessCtrl.adaptSpatialInfo(i),'%.2f')...
                        'Co' num2str(meanCorrDistCtrl.meanNZ(i),'%.2f')];                
                    plotFRProfIndNeuronIndTrial(gca,filteredSpikeArrayNormT,...
                        numSamples(1),i,indLaps,figTitle,1);
                else
                    disp(['Firng rate of neuron ' num2str(i) ' is too low: ' ...
                        num2str(mFRStructSessCtrl.mFR(i)) ' Hz']);
                    continue;
                end
            end
        end
        
    end
    
end

function plotFRProfIndNeuronIndTrial...
            (handle,filteredSpikeArrayNormT,numSamples,neuronNo,indLaps,figTitle,indSessBorder)
    numTr = length(indLaps);
    FRProfilePerTrial = zeros(numTr,numSamples);
    for i = 1:numTr
        if(isempty(filteredSpikeArrayNormT{i}))
            continue;
        end
%         disp(['Tr no. ' num2str(i)]);
        szArr = size(filteredSpikeArrayNormT{i},2);
        if(szArr > numSamples)
            szArr = numSamples;
        end
%         szArr = size(filteredSpikeArrayNormT{indLaps(i)},2);
        FRProfilePerTrial(i,1:szArr) = filteredSpikeArrayNormT{indLaps(i)}(neuronNo,1:numSamples); 
            %./max(filteredSpikeArrayNormT{indLaps(i)}(neuronNo,:));
    end
%     FRProfilePerTrial = FRProfilePerTrial/max(FRProfilePerTrial(:));
    h = imagesc(0:numSamples-1,1:numTr,FRProfilePerTrial);
    if(~isempty(indSessBorder))
        hold on
        for i = 1:length(indSessBorder)
            h = plot([0 numSamples],indSessBorder(i)*ones(1,2),'r');
            set(h,'LineWidth',1);
        end
    end    
            
    set(gca,'FontSize',8.0,'Box','on','XLim',[0 numSamples-1],'YLim',[0 numTr]);
    xlabel('Dist (mm)');
    ylabel('Trial No.');
    title(figTitle);
end
