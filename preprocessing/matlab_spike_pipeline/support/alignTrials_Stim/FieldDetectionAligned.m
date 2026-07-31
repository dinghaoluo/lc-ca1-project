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
    fileNameCorr = [fileName '_meanSpikesCorrTAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '_intT' ...
            num2str(intervalTSpInfo) '.mat'];
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
                'function "ConvSpikeTrain_Aligned" first.']);
        return;
    end
    load(fullPath,'filteredSpikeArrayRun','paramC','timeStep');
    
    fullPath = [path fileNamePeakFR];
    if(exist(fullPath) == 0)
        disp(['The peak firing rate file does not exist. Please call ',...
                'function "PeakFiringRate_Aligned" first.']);
        return;
    end
    load(fullPath,'pFRNonStimGoodStruct','trialNoNonStimGood');
    
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
    
    fullPath = [path fileNameCorr];
    if(exist(fullPath) == 0)
        disp(['The mean correlation in time file does not exist. Please call ',...
                'function "meanNeuronSpikeCorrT" first.']);
        return;
    end
    load(fullPath,'meanCorrTRun');
    
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo_smTr(path,fileName);
    end
    load(fullPath);
    
    paramF.minNumTr = 15;
    paramF.percNumActiveTr = 0.5;
    paramF.minCorr = 0.15;
    paramF.minSpInfo = 1;
    paramF.minInstFR = 0.15;
    paramF.lowThreFieldMeanInstFR = 0.05;
    paramF.timeBin = diff(timeStep(1:2))/1000;
    paramF.maxFieldWidth = 5.5/paramF.timeBin;
   
    fieldSpCorrSessNonStimGood = [];
    
    if(~isempty(meanCorrTRun))
        fieldSpCorrSessNonStimGood = FieldSpInfoCorr(...
            pFRNonStimGoodStruct,spatialInfoSessNonStimGood.adaptSpatialInfo,...
            meanCorrTRun,paramF);
        fieldSpCorrSessNonStimGood.FW = fieldSpCorrSessNonStimGood.FW*paramF.timeBin;
    end
    
    fullPath = [path fileNameFW];
    save(fullPath, 'fieldSpCorrSessNonStimGood','paramF');
    
    if(figureState == 2)
       %%% all the trials
       if(~isempty(fieldSpCorrSessNonStimGood))
       %% good trials
            count = 0;
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
                        'Co' num2str(meanCorrTRun.meanNZ(i),'%.2f')];                
                    plotFRProfIndNeuronIndTrial(gca,filteredSpikeArrayRun{i},...
                        intervalTSpInfo,indGoodLap,figTitle,indSessBorder,timeStep/1000);
                else
                    disp(['Firng rate of neuron ' num2str(i) ' is too low: ' ...
                        num2str(mFRStructNonStimGood.mFR(i)) ' Hz']);
                    continue;
                end
            end
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
