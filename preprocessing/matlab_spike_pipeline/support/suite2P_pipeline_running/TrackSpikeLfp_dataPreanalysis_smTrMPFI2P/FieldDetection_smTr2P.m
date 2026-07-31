function FieldDetection_smTr2P(path,fileName,onlyRun,figureState,spaceBin,mazeSess,intervalDSpInfo)
% Firing field detection based on correlation distance and spatial information
% path:         the path of the recording file
% fileName:     name of the recording file
% onlyRun:      1: only consider the time period when the animal is running 
%
% Example:
% FieldDetectionAligned('./','A111-20150301-01_DataStructure_mazeSection1_TrialType1',1,2)

    %%%%%%%%% check arguments
    if nargin<2
        disp('At least two arguments are needed for this function.');
        return;
    elseif nargin == 2
        figureState = 0;
        onlyRun = 1;
        mazeSess = 1;
        spaceBin = 20;
        intervalDSpInfo = 1800;
    elseif nargin == 3
        figureState = 0;   
        mazeSess = 1;
        spaceBin = 20;
        intervalDSpInfo = 1800;
    elseif nargin == 4
        mazeSess = 1;
        spaceBin = 20;
        intervalDSpInfo = 1800;
    elseif nargin == 5
        spaceBin = 20;
        intervalDSpInfo = 1800;
    elseif nargin == 6
        intervalDSpInfo = 1800;
    elseif nargin > 7
        disp('Too many input arguments.');
        return;
    end
    
    %%%%%%%%% initialize constants
    GlobalConst2P;
        
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    
    fileNameFW = [fileName '_FieldSpCorr_Run' num2str(onlyRun) '.mat'];
    
    fileNameConv = [fileName '_convSpikesDist' num2str(spaceBin) ...
                    'mm_Run' num2str(onlyRun) '.mat'];
    fileNamePeakFR = [fileName '_PeakFR' num2str(spaceBin) ...
                      'mm_Run' num2str(onlyRun) '.mat'];    
    fileNameFR = [fileName '_FR_Run' num2str(onlyRun) '.mat'];
    fileNameSpInfo = [fileName '_SpInfoAligned_msess' num2str(mazeSess) ...
                        '_Run' num2str(onlyRun) '.mat'];
    fileNameSpInfo = [fileName '_SpInfo_Run' num2str(onlyRun) '.mat'];
    fileNameCorr = [fileName '_meanSpikesCorrDist_Run' num2str(onlyRun) '_intD' ...
            num2str(intervalDSpInfo) '.mat'];
    fileNameInfo = [fileName '_Info.mat'];
    fileName = [fileName '.mat'];
         
    fullPath = [path fileNameConv];
    if(exist(fullPath) == 0)
        disp(['The firing profile file does not exist. Please call ',...
                'function "ConvSpikeTrainDistPar_smTr2P" first.']);
        return;
    end
    load(fullPath,'filteredSpikeArrayNormT','paramC');
    
    fullPath = [path fileNamePeakFR];
    if(exist(fullPath) == 0)
        disp(['The peak firing rate file does not exist. Please call ',...
                'function "PeakFiringRate_smTr2P" first.']);
        return;
    end
    load(fullPath,'pFRStructNormT');
    
    fullPath = [path fileNameFR];
    if(exist(fullPath) == 0)
        disp(['The mean firing rate file does not exist. Please call ',...
                'function "MeanFiringRate2P" first.']);
        return;
    end
    load(fullPath,'mFRStructSM');
    
    fullPath = [path fileNameSpInfo];
    if(exist(fullPath) == 0)
        disp(['The spatial information file does not exist. Please call ',...
                'function "GetFRMapInfo2P" first.']);
        return;
    end
    load(fullPath,'spatialInfoSess');
    spatialInfoSess = spatialInfoSess{mazeSess};
    
    fullPath = [path fileNameCorr];
    if(exist(fullPath) == 0)
        disp(['The mean correlation in time file does not exist. Please call ',...
                'function "meanNeuronSpikeCorrDist_smTr2P" first.']);
        return;
    end
    load(fullPath,'meanCorrDist');
    
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo_smTr(path,fileName);
    end
    load(fullPath);
    
    fullPath = [path fileName]; 
    if(exist(fullPath) == 0)
        disp('The recording file does not exist');
        return;
    end
    load(fullPath,'cluList');
        
    trackLen = unique(beh.trackLen(beh.mazeSess == mazeSess));
    paramF.minNumTr = 15;
    paramF.percNumActiveTr = 0.25;
    paramF.minCorr = 0.1;
    paramF.minSpInfo = 0.9;
    paramF.minHighSpInfoCorr = 0.038;
    paramF.minHighSpInfo =1.3;
    paramF.minCorrHigh = 0.1;
    paramF.minCorrHighSpInfo = 0.6;
    paramF.minInstFR = 0.01;    
    paramF.lowThreFieldMeanInstFR = 0.1;
    paramF.timeBin = timeStep;
    paramF.maxFieldWidth = trackLen*0.7;
    paramF.reboundCheckRegion = trackLen*0.3; % there is no clear peak within ReboundCheckRegion (mm) region before or after the field
    paramF.maxReboundMean = 0.7; % rebound mean inst FR should be smaller than maxReboundMean * meanInstFR
    paramF.reboundHeight = 0.4; % rebound height should be smaller than reboundHeight*peakFieldFR
    paramF.intervalTSpInfo = intervalDSpInfo;
    paramF.maxFieldWidth1 = trackLen*0.7; 
    
    fieldSpCorr = [];
        
    if(~isempty(meanCorrDist))
        meanCorrDist.meanGoodNZ = meanCorrDist.meanNZ;
        meanCorrDist.nGoodNonZeroTr = meanCorrDist.nNonZeroTr;
        fieldSpCorr = FieldSpInfoCorr2P(...
            pFRStructNormT,spatialInfoSess.adaptSpatialInfo,...
            meanCorrDist,paramF);
        fieldSpCorr.FW = fieldSpCorr.FW*paramF.timeBin;
    end
    
    fullPath = [path fileNameFW];
    save(fullPath, 'fieldSpCorr','paramF');
    fieldSpCorr
    
    if(figureState == 2)
       %%% all the trials
       neuronNo = 1:length(mFRStructSM.mFR);
       figNo = 1;
       if(~isempty(fieldSpCorr))
       %% good trials
            count = 0;
            indGoodLap = mFRStructSM.indLapList;
            [~,indGoodLap] = setdiff(indGoodLap,1:startTrNo);
            for i = neuronNo 
                if(mFRStructSM.mFR(i) > 0.01) %minFR
                   strNumField = '';
                   totTrialNo = 0;
                   fieldInfoTmp = getFieldInfoIndNeuron2P(i,fieldSpCorr); 
                    if(~isempty(fieldInfoTmp))
                        strNumField = [strNumField, ' ', num2str(size(fieldInfoTmp,1))];
                    else
                        strNumField = [strNumField, ' 0'];
                    end
                    
                            % get the field information
                    count = count + 1;

                    if(mod(count-1,16) == 0)
                        if(count > 1)
                            fullpath = [path fileName(1:end-4) '_DistRun' num2str(onlyRun) '_' num2str(figNo)];
                            print('-painters', '-dpdf', fullpath, '-r600')
                            savefig([fullpath '.fig']);
                            figNo = figNo+1;
                        end
    
                        [figNew,pos] = CreateFig2P();
                        set(0,'Units','pixels') 
                        figTitle = 'All the trials';
                        set(figure(figNew),'OuterPosition',...
                            [pos(1) pos(2)-300 pos(3)*1.6 pos(4)*1.6],'Name',figTitle)
                    end

                    subplot(4,4,mod(count-1,16)+1)

                    figTitle = ['Neu ' num2str(i) '(' num2str(cluList.localClu(i)-1) ') Field' strNumField ' G Sp' ...
                        num2str(spatialInfoSess.adaptSpatialInfo(i),'%.2f')...
                        'Co' num2str(meanCorrDist.meanGoodNZ(i),'%.2f')];  
                    
                    filteredSpikeArray = pFRStructNormT.filteredSpArrAll{i};
                    plotFRProfIndNeuronIndTrial(gca,...
                        filteredSpikeArray./max(filteredSpikeArray,[],2),...
                        intervalDSpInfo,indGoodLap,figTitle,...
                        paramC.spaceSteps{mazeSess});
                else
                    disp(['Firng rate of neuron ' num2str(i) ' is too low: ' ...
                        num2str(mFRStructSM.mFR(i)) ' Hz']);
                    continue;
                end
            end
            
            fullpath = [path fileName(1:end-4) '_DistRun' num2str(onlyRun) '_' num2str(figNo)];
            print('-painters', '-dpdf', fullpath, '-r600')
            savefig([fullpath '.fig']);
       end        
    end    
end

function plotFRProfIndNeuronIndTrial...
            (handle,filteredSpikeArrayT,intervalT,indLaps,figTitle,timeStep)
    numTr = length(indLaps);
    h = imagesc(timeStep,1:numTr,filteredSpikeArrayT(indLaps,:));
            
    set(gca,'FontSize',8.0,'Box','on','YLim',[0 numTr],'XLim',[0 intervalT]);
    xlabel('Dist (mm)');
    ylabel('Trial No.');
    title(figTitle);
end
