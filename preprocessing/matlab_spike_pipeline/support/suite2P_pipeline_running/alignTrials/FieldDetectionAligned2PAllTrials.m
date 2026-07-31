function FieldDetectionAligned2PAllTrials(path,fileName,onlyRun,figureState,mazeSess,intervalTSpInfo)
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
        intervalTSpInfo = 10;
    elseif nargin == 3
        figureState = 0;   
        mazeSess = 1;
        intervalTSpInfo = 10;
    elseif nargin == 4
        mazeSess = 1;
        intervalTSpInfo = 10;
    elseif nargin == 5
        intervalTSpInfo = 10;
    elseif nargin > 6
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
    
    fileNameFW = [fileName '_FieldSpCorrAlignedAllTr_Run' num2str(mazeSess) ...
                        '_Run' num2str(onlyRun) '.mat'];
    
    fileNameConv = [fileName '_convSpikesAlignedRun_msess' num2str(mazeSess) ...
                        '_Run' num2str(onlyRun) '.mat'];
    fileNamePeakFR = [fileName '_PeakFRAlignedRun_msess' num2str(mazeSess) ...
                        '_Run' num2str(onlyRun) '.mat'];
    fileNamePeakFRCue = [fileName '_PeakFRAlignedCue_msess' num2str(mazeSess) ...
                        '_Run' num2str(onlyRun) '.mat'];
    fileNameFR = [fileName '_FRAlignedRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    fileNameSpInfo = [fileName '_SpInfoAligned_msess' num2str(mazeSess) ...
                        '_Run' num2str(onlyRun) '.mat'];
    fileNameSpInfoCue = [fileName '_SpInfoAlignedCue_msess' num2str(mazeSess) ...
                        '_Run' num2str(onlyRun) '.mat'];
    fileNameCorr = [fileName '_meanSpikesCorrTAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalTSpInfo) '.mat'];
    fileNameInfo = [fileName '_Info.mat'];
    fileNameBeh = [fileName '_behPar_msess' num2str(mazeSess) '.mat'];
    fileName = [fileName '.mat'];
         
    fullPath = [path fileNameConv];
    if(exist(fullPath) == 0)
        disp(['The firing profile file does not exist. Please call ',...
                'function "ConvSpikeTrain_Aligned" first.']);
        return;
    end
    load(fullPath,'filteredSpikeArrayRun','paramC');
    
    fullPath = [path fileNamePeakFR];
    if(exist(fullPath) == 0)
        disp(['The peak firing rate file does not exist. Please call ',...
                'function "PeakFiringRate_Aligned" first.']);
        return;
    end
    load(fullPath,'pFRNonStimGoodStruct','trialNoNonStimGood');
    
    fullPath = [path fileNamePeakFRCue];
    if(exist(fullPath) == 0)
        disp(['The peak firing rate file does not exist. Please call ',...
                'function "PeakFiringRate_Aligned" first.']);
        return;
    end
    load(fullPath,'pFRNonStimGoodStructCue');
    
    fullPath = [path fileNameFR];
    if(exist(fullPath) == 0)
        disp(['The mean firing rate file does not exist. Please call ',...
                'function "MeanFiringRateAligned" first.']);
        return;
    end
    load(fullPath,'mFRStructNonStimGood','mFRStructNonStimBad');
    
    fullPath = [path fileNameSpInfo];
    if(exist(fullPath) == 0)
        disp(['The spatial information file does not exist. Please call ',...
                'function "GetFRMapInfo_Aligned" first.']);
        return;
    end
    load(fullPath,'spatialInfoSessNonStimGood');
    spatialInfoSessNonStimGoodRun = spatialInfoSessNonStimGood;
    
    fullPath = [path fileNameSpInfoCue];
    if(exist(fullPath) == 0)
        disp(['The spatial information file does not exist. Please call ',...
                'function "GetFRMapInfo_AlignedCue" first.']);
        return;
    end
    load(fullPath,'spatialInfoSessNonStimGood');
    spatialInfoSessNonStimGoodCue = spatialInfoSessNonStimGood;
    
    fullPath = [path fileNameCorr];
    if(exist(fullPath) == 0)
        disp(['The mean correlation in time file does not exist. Please call ',...
                'function "meanNeuronSpikeCorrT" first.']);
        return;
    end
    load(fullPath,'meanCorrTRun','meanCorrTCue');
    
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo_smTr(path,fileName);
    end
    load(fullPath);
    
    fullPath = [path fileNameBeh]; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    fullPath = [path fileName]; 
    if(exist(fullPath) == 0)
        disp('The recording file does not exist');
        return;
    end
    load(fullPath,'cluList');
            
    paramF.minNumTr = 10;
    paramF.percNumActiveTr = 0.15;
    paramF.minCorr = 0.14;
    paramF.minSpInfo = 1.0;
    paramF.minHighSpInfoCorr = 0.08;
    paramF.minHighSpInfo = 1.27;
    paramF.minCorrHigh = 0.17;
    paramF.minCorrHighSpInfo = 0.4;
    paramF.minInstFR = 0.005;    
    paramF.lowThreFieldMeanInstFR = 0.1;
    paramF.timeBin = timeStep;
    paramF.maxFieldWidth = 4.5/paramF.timeBin;
    paramF.reboundCheckRegion = 0.5/paramF.timeBin; % there is no clear peak within ReboundCheckRegion (mm) region before or after the field
    paramF.maxReboundMean = 0.7; % rebound mean inst FR should be smaller than maxReboundMean * meanInstFR
    paramF.reboundHeight = 0.5; % rebound height should be smaller than reboundHeight*peakFieldFR
    paramF.intervalTSpInfo = intervalTSpInfo/paramF.timeBin;
    paramF.maxFieldWidth1 = mean(behPar.numSamplesRun(mFRStructNonStimGood.indLapList))*0.75/sampleFq/paramF.timeBin; 
    
    fieldSpCorrSessNonStimGood = [];
    
    if(length(mFRStructNonStimGood.indLapList) ~= length(trialNoNonStimGood))
        disp('Number of trials does not match between mFR and pFR. Please check');
    end
    
    if(~isempty(meanCorrTRun))
        meanCorrTRun.meanGoodNZ = meanCorrTRun.meanNZ;
        meanCorrTRun.nGoodNonZeroTr = meanCorrTRun.nNonZeroTr;
        fieldSpCorrSessNonStimGood = FieldSpInfoCorrAligned2P(...
            pFRNonStimGoodStruct,spatialInfoSessNonStimGoodRun.adaptSpatialInfo,...
            meanCorrTRun,paramF);
        
        fieldSpCorrSessNonStimGood.FW = fieldSpCorrSessNonStimGood.FW*paramF.timeBin;
    end
    
    if(~isempty(meanCorrTCue))
        meanCorrTCue.meanGoodNZ = meanCorrTCue.meanNZ;
        meanCorrTCue.nGoodNonZeroTr = meanCorrTCue.nNonZeroTr;
        fieldSpCorrSessNonStimGoodCue = FieldSpInfoCorrAligned2P(...
            pFRNonStimGoodStructCue,spatialInfoSessNonStimGoodCue.adaptSpatialInfo,...
            meanCorrTCue,paramF); %% changed from pFRNonStimGoodStruct to pFRNonStimCueGoodStruct by Yingxue on 1/21/2022
        fieldSpCorrSessNonStimGoodCue.FW = fieldSpCorrSessNonStimGoodCue.FW*paramF.timeBin;
    end
    
    fullPath = [path fileNameFW];
    save(fullPath, 'fieldSpCorrSessNonStimGood','fieldSpCorrSessNonStimGoodCue','paramF');
    fieldSpCorrSessNonStimGood
    fieldSpCorrSessNonStimGoodCue
    
    if(figureState == 2)
       %%% all the trials
       neuronNo = 1:length(mFRStructNonStimGood.mFR);
       figNo = 1;
       if(~isempty(fieldSpCorrSessNonStimGood))
       %% good trials
            count = 0;
            indGoodLap = mFRStructNonStimGood.indLapList;
            indGoodLap = setdiff(indGoodLap,1:startTrNo);
            indSessBorder = length(indGoodLap);
            indBadLap = [];
            if(isfield(mFRStructNonStimBad,'indLapList'))
                indBadLap = mFRStructNonStimBad.indLapList;
                indBadLap = setdiff(indBadLap,1:startTrNo);
            end
            indGoodLap = [indGoodLap; indBadLap];
            for i = neuronNo 
                if(mFRStructNonStimGood.mFR(i) > 0.01) %minFR
                   strNumField = '';
                   totTrialNo = 0;
                   fieldInfoTmp = getFieldInfoIndNeuron2P(i,fieldSpCorrSessNonStimGood); 
                    if(~isempty(fieldInfoTmp))
                        strNumField = [strNumField, ' ', num2str(size(fieldInfoTmp,1))];
                    else
                        strNumField = [strNumField, ' 0'];
                    end
                   
                    count = count + 1;

                    if(mod(count-1,16) == 0)
                        if(count > 1)
                            fullpath = [path fileName(1:end-4) '_Run' num2str(onlyRun) '_' num2str(figNo)];
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
                        num2str(spatialInfoSessNonStimGoodRun.adaptSpatialInfo(i),'%.2f')...
                        'Co' num2str(meanCorrTRun.meanGoodNZ(i),'%.2f')];                
                    plotFRProfIndNeuronIndTrial(gca,...
                        filteredSpikeArrayRun{i}./max(filteredSpikeArrayRun{i},[],2),...
                        intervalTSpInfo,indGoodLap',figTitle,indSessBorder,...
                        0:timeStep:timeStep*size(filteredSpikeArrayRun{i},2));
                else
                    disp(['Firng rate of neuron ' num2str(i) ' is too low: ' ...
                        num2str(mFRStructNonStimGood.mFR(i)) ' Hz']);
                    continue;
                end
            end
            
            fullpath = [path fileName(1:end-4) '_Run' num2str(onlyRun) '_' num2str(figNo)];
            print('-painters', '-dpdf', fullpath, '-r600')
            savefig([fullpath '.fig']);
        end
        
    end
    
    a=1;
    
end

function plotFRProfIndNeuronIndTrial...
            (handle,filteredSpikeArrayT,intervalT,indLaps,figTitle,indSessBorder,timeStep)
    numTr = length(indLaps);
    numTotalTr = size(filteredSpikeArrayT,1);
    indLapsTmp = [indLaps setdiff(1:numTotalTr,indLaps)];
    h = imagesc(timeStep,1:numTotalTr,filteredSpikeArrayT(indLapsTmp,:));
    if(~isempty(indSessBorder))
        hold on
        for i = 1:length(indSessBorder)
            h = plot([0 intervalT],indSessBorder(i)*ones(1,2),'r');
            set(h,'LineWidth',1);
        end
    end    
            
    set(gca,'FontSize',8.0,'Box','on','YLim',[0 numTotalTr],'XLim',[0 intervalT]);
    xlabel('Time (s)');
    ylabel('Trial No.');
    title(figTitle);
end
