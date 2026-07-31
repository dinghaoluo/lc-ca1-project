function probDistrPopActivityCorr(pathAnal,recType,spaceBin,onlyRun)
% for selected neurons with high corr., calculate the population level
% trial to trial correlation
% default analysis path ['Z:\Raphael_tests\mice_expdata\Analysis\'];
% e.g. probDistrPopActivityCorr('Z:\Raphael_tests\mice_expdata\Analysis\',1,20,1)

    RecordingList;
    if(recType == 1) 
        recListPath = listRecordingsPassiveLickPath;
        recListFileName = listRecordingsPassiveLickFileName;
        recListSess = mazeSessionPassiveLick;
    elseif(recType == 2)
        recListPath = listRecordingsActiveLickPath;
        recListFileName = listRecordingsActiveLickFileName;
        recListSess = mazeSessionActiveLick;
    elseif(recType == 3)
        % passive no cue
        recListPath = listRecordingsNoCuePath;
        recListFileName = listRecordingsNoCueFileName;
        recListSess = mazeSessionNoCue;
    end
    recListLen = size(recListPath,1);
    
    load([pathAnal 'recList' num2str(recType) '_activityCorrNeu_Run1.mat']);
    
    recListPopActivityCorr.indPyr = [];
    recListPopActivityCorr.indLaps = [];
    recListPopActivityCorr.corrMatx = [];
    recListPopActivityCorr.meanCorr = zeros(1,recListLen);
    
    for i = 1:recListLen
        path = recListPath(i,:);
        fileName = recListFileName(i,:);
        fullPath = [path fileName ...
            '_convSpikesDist' num2str(spaceBin) ...
                    'mm_Run' num2str(onlyRun) '.mat'];
        if(exist(fullPath,'file') == 0)
            disp(['The firing profile file does not exist. Try to run the',...
                    'function again with fileState = 0.']);
            return;
        end
        load(fullPath,'filteredSpikeArrayNormT');
        
%         if recType == 2
%             fileNameCluster = [fileName '_kmeanBeh.mat'];
%             fullPath = [path fileNameCluster];
%             if(exist(fullPath) == 0)
%                 disp('The behavior cluster file does not exist.');
%                 disp(fullPath)
%                 return;
%             end
%             load(fullPath);
%             
%             [indCluGood,indCluBad] = HandSelectedCluFromBehKMean(fileName);
%             
%             indTr = indLaps(kmeanClu == indCluGood);
%         else
        fileNameInfo = [fileName '_Info.mat'];     
        fullPath = [path fileNameInfo];
        if(exist(fullPath) == 0)
            BasicInfo(path,fileName);
        end
        load(fullPath,'beh');
        fullPath = [path fileName '_behPar_msess' num2str(recListSess(i)) '.mat']; 
        if(exist(fullPath) == 0)
            disp('The _behPar file does not exist');
            return;
        end
        load(fullPath);
        indTr = find(beh.mazeSess == recListSess(i));
        indTr = indTr(behPar.indTrBadBeh == 0);
        
        recListPopActivityCorr.indLaps{i} = indTr;
        
        indGoodNeu = recListActivityCorr.indPyrGoodNeu(...
            recListActivityCorr.indRecGoodNeu == i);
        recListPopActivityCorr.indPyr{i} = indGoodNeu; 
        
        if(length(indGoodNeu) < 3)
            recListPopActivityCorr.corrMatx{i} = [];
            continue;
        end
        
        filteredSpikeArrayTmp = zeros(length(indTr),length(indTr));
        
        count = 0;
        for j = 1:length(indTr)
            for k = j+1:length(indTr)
                count = count+1;
                filteredSpikeArrayTmp(j,k) = ...
                    corr2(filteredSpikeArrayNormT{indTr(j)}(indGoodNeu,:),...
                        filteredSpikeArrayNormT{indTr(k)}(indGoodNeu,:)); 
                if(isnan(filteredSpikeArrayTmp(j,k)))
                    filteredSpikeArrayTmp(j,k) = 0;
                end
            end
        end
        recListPopActivityCorr.corrMatx{i} = filteredSpikeArrayTmp;
        recListPopActivityCorr.meanCorr(i) = ...
            sum(sum(triu(filteredSpikeArrayTmp,1)))/count; 
    end
    
    fullPath = [pathAnal 'popActivityCorr_RecList' num2str(recType) '_Run'...
        num2str(onlyRun) '.mat'];
    save(fullPath,'recListPopActivityCorr');
end
