function FieldDetectionAligned(path,fileName,onlyRun,figureState,mazeSess,intervalTSpInfo)
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
    GlobalConst;
        
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    
    fileNameFW = [fileName '_FieldSpCorrAligned_Run' num2str(mazeSess) ...
                        '_Run' num2str(onlyRun) '.mat'];
    
    fileNameConv = [fileName '_convSpikesAligned_msess' num2str(mazeSess) ...
                        '_Run' num2str(onlyRun) '.mat'];
    fileNamePeakFR = [fileName '_PeakFRAligned_msess' num2str(mazeSess) ...
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
                'function "ConvSpikeTrain_Aligned" first.']);
        return;
    end
    load(fullPath,'filteredSpikeArrayRun','filteredSpikeArrayCue','paramC','timeStep');
    
    fullPath = [path fileNamePeakFR];
    if(exist(fullPath) == 0)
        disp(['The peak firing rate file does not exist. Please call ',...
                'function "PeakFiringRate_Aligned" first.']);
        return;
    end
    load(fullPath,'pFRNonStimGoodStruct','pFRNonStimCueGoodStruct','trialNoNonStimGood');
    
    fullPath = [path fileNameFR];
    if(exist(fullPath) == 0)
        disp(['The mean firing rate file does not exist. Please call ',...
                'function "MeanFiringRateAligned" first.']);
        return;
    end
    load(fullPath,'mFRStructNonStimGood');
    
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
        
    paramF.minNumTr = 15;
    paramF.percNumActiveTr = 0.4;
    paramF.minCorr = 0.12;
    paramF.minSpInfo = 1;
    paramF.minCorrHighSpInfo = 0.09;
    paramF.minHighSpInfo = 2.5;
    paramF.minInstFR = 0.09;    
    paramF.lowThreFieldMeanInstFR = 0.1;
    paramF.timeBin = diff(timeStep(1:2))/1250;
    paramF.maxFieldWidth = 4.5/paramF.timeBin;
    paramF.reboundCheckRegion = 1.7/paramF.timeBin; % there is no clear peak within ReboundCheckRegion (mm) region before or after the field
    paramF.maxReboundMean = 0.18; % rebound mean inst FR should be smaller than maxReboundMean * meanInstFR
    paramF.reboundHeight = 0.45; % rebound height should be smaller than reboundHeight*peakFieldFR
    paramF.intervalTSpInfo = intervalTSpInfo/paramF.timeBin;
    paramF.maxFieldWidth1 = mean(behPar.numSamplesRun(mFRStructNonStimGood.indLapList))*0.7/sampleFq/paramF.timeBin; 
    
    fieldSpCorrSessNonStimGood = [];
    
    if(length(mFRStructNonStimGood.indLapList) ~= length(trialNoNonStimGood))
        disp('Number of trials does not match between mFR and pFR. Please check');
    end
    
    if(~isempty(meanCorrTRun))
        fieldSpCorrSessNonStimGood = FieldSpInfoCorrAligned(...
            pFRNonStimGoodStruct,spatialInfoSessNonStimGoodRun.adaptSpatialInfo,...
            meanCorrTRun,paramF);
        fieldSpCorrSessNonStimGood.FW = fieldSpCorrSessNonStimGood.FW*paramF.timeBin;
    end
    
    if(~isempty(meanCorrTCue))
        fieldSpCorrSessNonStimGoodCue = FieldSpInfoCorrAligned(...
            pFRNonStimCueGoodStruct,spatialInfoSessNonStimGoodCue.adaptSpatialInfo,...
            meanCorrTCue,paramF); %% changed from pFRNonStimGoodStruct to pFRNonStimCueGoodStruct by Yingxue on 1/21/2022
        fieldSpCorrSessNonStimGoodCue.FW = fieldSpCorrSessNonStimGoodCue.FW*paramF.timeBin;
    end
    
    fullPath = [path fileNameFW];
    save(fullPath, 'fieldSpCorrSessNonStimGood','fieldSpCorrSessNonStimGoodCue','paramF');
    fieldSpCorrSessNonStimGood
    fieldSpCorrSessNonStimGoodCue
    
    if(figureState == 2)
       %%% all the trials
       if(~isempty(fieldSpCorrSessNonStimGood))
       %% good trials
            count = 0;
            numFig = 0;
            for i = 1:rec.numNeurons  
                if(mFRStructNonStimGood.mFR(i) > 0.05) %minFR
                   strNumField = '';
                   totTrialNo = 0;
                   fieldInfoTmp = getFieldInfoIndNeuron(i,fieldSpCorrSessNonStimGood); 
                    if(~isempty(fieldInfoTmp))
                        strNumField = [strNumField, ' ', num2str(size(fieldInfoTmp,1))];
                    else
                        strNumField = [strNumField, ' 0'];
                    end
                    indGoodLap = mFRStructNonStimGood.indLapList;
                    indSessBorder = totTrialNo+length(indGoodLap);
                            % get the field information
                    count = count + 1;

                    if(mod(count-1,16) == 0)
%                         if(numFig > 0)
%                             fileName1 = [path fileNameFW(1:end-4) num2str(numFig)];
%                             saveas(gcf,fileName1);
%                             print('-painters', '-dpng', fileName1, '-r600') 
%                         end
%                         numFig = numFig + 1;
                        [figNew,pos] = CreateFig();
                        set(0,'Units','pixels') 
                        figTitle = 'All the trials';
                        set(figure(figNew),'OuterPosition',...
                            [pos(1) pos(2)-500 pos(3)*2 pos(4)*2.2],'Name',figTitle)
                    end

                    subplot(4,4,mod(count-1,16)+1)

                    figTitle = ['Neu ' num2str(i) '(' num2str(cluList.shank(i))...
                        ' ' num2str(cluList.localClu(i)) ') Field' strNumField ' G Sp' ...
                        num2str(spatialInfoSessNonStimGood.adaptSpatialInfo(i),'%.2f')...
                        'Co' num2str(meanCorrTRun.meanGoodNZ(i),'%.2f')];                
                    plotFRProfIndNeuronIndTrial(gca,filteredSpikeArrayRun{i},...
                        intervalTSpInfo,indGoodLap',figTitle,indSessBorder,timeStep/1000);
%                     plotFRProfIndNeuronIndTrial(gca,filteredSpikeArrayCue{i},...
%                         intervalTSpInfo,indGoodLap',figTitle,indSessBorder,timeStep/1000);
                else
                    disp(['Firng rate of neuron ' num2str(i) ' is too low: ' ...
                        num2str(mFRStructNonStimGood.mFR(i)) ' Hz']);
                    continue;
                end
            end
            
%             fileName1 = [path fileNameFW(1:end-4) num2str(numFig)];
%             saveas(gcf,fileName1);
%             print('-painters', '-dpng', fileName1, '-r600') 
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
