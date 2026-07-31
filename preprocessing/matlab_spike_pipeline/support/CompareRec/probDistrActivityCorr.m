function probDistrActivityCorr(pathAnal,recType,onlyRun)
% plot the distribution of neuronal activity correlation over all the trials
% for each neuron
% default analysis path ['Z:\Raphael_tests\mice_expdata\Analysis\'];
% e.g. probDistrActivityCorr('Z:\Raphael_tests\mice_expdata\Analysis\',1,1)

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
    
    param.minFR = 0.3;
    param.maxFR = 5;
    param.minCorr = 0.015;
    param.percActiveTr = 0.4;
    
    recListActivityCorr.percGoodPyrNeuPerRec = zeros(1,recListLen);
    recListActivityCorr.meanCorrActiveTrPerNeuPerRec = zeros(1,recListLen);
    recListActivityCorr.stdCorrActiveTrPerNeuPerRec = zeros(1,recListLen);
    recListActivityCorr.meanCorrActiveTrGoodNeuPerRec = zeros(1,recListLen);
    recListActivityCorr.stdCorrActiveTrGoodNeuPerRec = zeros(1,recListLen);
    recListActivityCorr.indRec = [];
    recListActivityCorr.indPyrNeu = [];
    recListActivityCorr.indRecGoodNeu = [];
    recListActivityCorr.indPyrGoodNeu = [];
    recListActivityCorr.meanCorrActiveTrPerNeu = [];
    recListActivityCorr.meanCorrAllTrPerNeu = [];
    recListActivityCorr.meanCorrGoodTrPerNeu = [];
    recListActivityCorr.meanCorrTActiveTrPerNeu = [];
    recListActivityCorr.meanCorrTAllTrPerNeu = [];
    recListActivityCorr.meanCorrTGoodTrPerNeu = [];
    recListActivityCorr.meanFRActiveTrPerNeu = [];
    recListActivityCorr.spatialInfoPerNeu = [];
    recListActivityCorr.adaptSpatialInfoPerNeu = [];
    recListActivityCorr.sparsityPerNeu = [];
    recListActivityCorr.SNRPerNeu = [];
    
    for i = 1:recListLen
        if(recType <= 3)
            fullPath = [recListPath(i,:) recListFileName(i,:) ...
                '_behPar_msess' num2str(recListSess(i)) '.mat']; 
            if(exist(fullPath) == 0)
                disp('The _behPar file does not exist');
                return;
            end
            load(fullPath);
            
            %% load activity correlation in distance file
            fullPath = [recListPath(i,:) recListFileName(i,:) ...
                    '_NeuActivityCorr_Run' num2str(onlyRun) '.mat'];
            if(exist(fullPath,'file') == 2)
                load(fullPath,'activityCorr');    
            end
            if(iscell(activityCorr))
                activityCorr = activityCorr{1};
            end
            
            %% indices of good trials
            indBadTr = behPar.totRunLenT > 13 | behPar.numRun > 10 | ...
                behPar.totStopLenT > 2 | behPar.rewarded == -1;
                % remove the nonstopping criteria 
            indGoodTr = find(indBadTr == 0);
            [indGoodTr,indIndGoodTr] = intersect(activityCorr.indLaps,indGoodTr);
            numGoodTr = length(indGoodTr);
            totNumCorr = (numGoodTr*numGoodTr-numGoodTr)/2;
            recListActivityCorr.indGoodTr = indGoodTr;
            
            %% average correlation dist of good trials
            recListActivityCorr.meanCorrGoodTrPerNeuPerRec{i} = ...
                zeros(1,length((activityCorr.indPyrNeu)));
            for n = 1:length(activityCorr.indPyrNeu)
                corrMatxPerNeuTmp = activityCorr.corrMatxPerNeu{n}...
                    (indIndGoodTr,indIndGoodTr);
                recListActivityCorr.meanCorrGoodTrPerNeuPerRec{i}(n) = ...
                    sum(sum(triu(corrMatxPerNeuTmp,1)))/totNumCorr;
            end
            
            %% load activity correlation in time file
            fullPath = [recListPath(i,:) recListFileName(i,:) ...
                    '_NeuActivityCorrT_Run' num2str(onlyRun) '.mat'];
            if(exist(fullPath,'file') == 2)
                load(fullPath,'activityCorrT');    
            end
            if(iscell(activityCorrT))
                activityCorrT = activityCorrT{1};
            end
            %% average correlation Time of good trials
            recListActivityCorr.meanCorrTGoodTrPerNeuPerRec{i} = ...
                zeros(1,length((activityCorr.indPyrNeu)));
            for n = 1:length(activityCorrT.indPyrNeu)
                corrMatxPerNeuTmp = activityCorrT.corrMatxPerNeu{n}...
                    (indIndGoodTr,indIndGoodTr);
                recListActivityCorr.meanCorrTGoodTrPerNeuPerRec{i}(n) = ...
                    sum(sum(triu(corrMatxPerNeuTmp,1)))/totNumCorr;
            end
            
%             fullPath = [recListPath(i,:) recListFileName(i,:) ...
%                     '_NeuActivitySimT_Run' num2str(onlyRun) '.mat'];
%             if(exist(fullPath,'file') == 2)
%                 load(fullPath,'spikeTrainSimT');    
%             end
            
            fullPath = [recListPath(i,:) recListFileName(i,:) ...
                    '_SpInfo_Run' num2str(onlyRun) '.mat'];
            if(exist(fullPath,'file') == 2)
                load(fullPath,'spatialInfoSess');     
            end
            spatialInfo = spatialInfoSess{1};
        else  %% need to modify this part for corrT and other measures
            fullPath = [recListPath(i,:) recListFileName(i,:) ...
                    '_ActivityCorrKluBeh_Run' num2str(onlyRun) '.mat'];
            if(exist(fullPath,'file') == 2)
                load(fullPath,'behBestCluActivityCorr');    
            end
            activityCorr = behBestCluActivityCorr;
            
            fullPath = [recListPath(i,:) recListFileName(i,:) ...
                    '_SpInfoKluBeh_Run' num2str(onlyRun) '.mat'];
            if(exist(fullPath,'file') == 2)
                load(fullPath,'spatialInfoKluBehGood');    
            end
            spatialInfo = spatialInfoKluBehGood;
        end
        
        indNeu = activityCorr.indPyrNeu;
        if(recType <= 3)
            meanFR = spatialInfo.meanFR(indNeu);
        else
            meanFR = spatialInfo.meanFR;
        end
        indGoodFR = meanFR > param.minFR & meanFR < param.maxFR;
        indGoodNeu = meanFR > param.minFR & meanFR < param.maxFR & ...
            activityCorr.meanCorrActiveTrPerNeu >= param.minCorr & ...
            activityCorr.activeTr/length(activityCorr.indLaps) > param.percActiveTr;
        
        recListActivityCorr.percGoodPyrNeuPerRec(i) = ...
            sum(indGoodNeu)/sum(indGoodFR);
        recListActivityCorr.meanCorrActiveTrPerNeuPerRec(i) = ...
            mean(activityCorr.meanCorrActiveTrPerNeu(indGoodFR));
        recListActivityCorr.stdCorrActiveTrPerNeuPerRec(i) = ...
            std(activityCorr.meanCorrActiveTrPerNeu(indGoodFR));
        
        if(sum(indGoodNeu) > 0)
            recListActivityCorr.meanCorrActiveTrGoodNeuPerRec(i) = ...
                mean(activityCorr.meanCorrActiveTrPerNeu(indGoodNeu));
            recListActivityCorr.stdCorrActiveTrGoodNeuPerRec(i) = ...
                std(activityCorr.meanCorrActiveTrPerNeu(indGoodNeu));
        end
        
        recListActivityCorr.indRec = ...
            [recListActivityCorr.indRec ones(1,sum(indGoodFR))*i];        
        recListActivityCorr.indPyrNeu = ...
            [recListActivityCorr.indPyrNeu indNeu(indGoodFR)];
        
        recListActivityCorr.indRecGoodNeu = ...
            [recListActivityCorr.indRecGoodNeu ones(1,sum(indGoodNeu))*i];     
        recListActivityCorr.indPyrGoodNeu = ...
            [recListActivityCorr.indPyrGoodNeu indNeu(indGoodNeu)];
        
        recListActivityCorr.meanCorrActiveTrPerNeu = ...
            [recListActivityCorr.meanCorrActiveTrPerNeu ...
            activityCorr.meanCorrActiveTrPerNeu(indGoodFR)];
        recListActivityCorr.meanCorrAllTrPerNeu = ...
            [recListActivityCorr.meanCorrAllTrPerNeu ...
            activityCorr.meanCorrAllTrPerNeu(indGoodFR)];
        recListActivityCorr.meanCorrGoodTrPerNeu = ...
            [recListActivityCorr.meanCorrGoodTrPerNeu ...
            recListActivityCorr.meanCorrGoodTrPerNeuPerRec{i}(indGoodFR)];
        
        %% corrT
        recListActivityCorr.meanCorrTActiveTrPerNeuPerRec(i) = ...
            mean(activityCorrT.meanCorrActiveTrPerNeu(indGoodFR));
        recListActivityCorr.meanCorrTActiveTrPerNeu = ...
            [recListActivityCorr.meanCorrTActiveTrPerNeu ...
            activityCorrT.meanCorrActiveTrPerNeu(indGoodFR)];
        recListActivityCorr.meanCorrTAllTrPerNeu = ...
            [recListActivityCorr.meanCorrTAllTrPerNeu ...
            activityCorrT.meanCorrAllTrPerNeu(indGoodFR)];
        recListActivityCorr.meanCorrTGoodTrPerNeu = ...
            [recListActivityCorr.meanCorrTGoodTrPerNeu ...
            recListActivityCorr.meanCorrTGoodTrPerNeuPerRec{i}(indGoodFR)];
        
%         %% simT
%         recListActivityCorr.meanSimTActiveTrPerNeuPerRec(i) = ...
%             mean(spikeTrainSimT.meanCorrActiveTrPerNeu(indGoodFR));
%         recListActivityCorr.meanCorrTActiveTrPerNeu = ...
%             [recListActivityCorr.meanCorrTActiveTrPerNeu ...
%             activityCorrT.meanCorrActiveTrPerNeu(indGoodFR)];
%         recListActivityCorr.meanCorrTAllTrPerNeu = ...
%             [recListActivityCorr.meanCorrTAllTrPerNeu ...
%             activityCorrT.meanCorrAllTrPerNeu(indGoodFR)];
        
        recListActivityCorr.meanFRActiveTrPerNeu = ...
            [recListActivityCorr.meanFRActiveTrPerNeu ...
            meanFR(indGoodFR)];
        
        recListActivityCorr.spatialInfoPerNeu = ...
            [recListActivityCorr.spatialInfoPerNeu ...
            spatialInfo.spatialInfo(indGoodFR)];
        recListActivityCorr.adaptSpatialInfoPerNeu = ...
            [recListActivityCorr.adaptSpatialInfoPerNeu ...
            spatialInfo.adaptSpatialInfo(indGoodFR)];
        recListActivityCorr.sparsityPerNeu = ...
            [recListActivityCorr.sparsityPerNeu ...
            spatialInfo.sparsity(indGoodFR)];
        recListActivityCorr.SNRPerNeu = ...
            [recListActivityCorr.SNRPerNeu ...
            spatialInfo.SNR(indGoodFR)];
        
    end
    
    if(exist(pathAnal,'dir') == 0)
        mkdir(pathAnal);
    end
    fullPath = [pathAnal 'recList' num2str(recType)...
        '_activityCorrNeu_Run' num2str(onlyRun) '.mat'];
    save(fullPath,'recListActivityCorr','param','recListPath','recListFileName');
        
    % activity correlation over dist
    xlabelTmp = 'Corr. neuronal activity over good trials';
    if(recType == 1)
        titleTmp = ['Passive licking (num. neurons ' ...
            num2str(length(recListActivityCorr.indPyrNeu)) ')'];
    elseif(recType == 2)
        titleTmp = ['Active licking (num. neurons ' ...
            num2str(length(recListActivityCorr.indPyrNeu)) ')'];
    elseif(recType == 3)
        titleTmp = ['Passive no cue (num. neurons ' ...
            num2str(length(recListActivityCorr.indPyrNeu)) ')'];
    end
%     plotHist(recListActivityCorr.meanCorrGoodTrPerNeu,-0.025:0.0025:0.35,1,...
%         xlabelTmp,titleTmp);
%     if(recType == 1)
%         print([pathAnal 'CorrActivityPassive'],'-dpdf','-r600');
%     elseif(recType == 2)
%         print([pathAnal 'CorrActivityActive'],'-dpdf','-r600');
%     elseif(recType == 3)
%         print([pathAnal 'CorrActivityPassiveNoCue'],'-dpdf','-r600');
%     end
%     
%     % activity correlation over time
%     xlabelTmp = 'CorrT. neuronal activity over good trials';
% 
%     plotHist(recListActivityCorr.meanCorrTGoodTrPerNeu,-0.025:0.0025:0.35,1,...
%         xlabelTmp,titleTmp);
%     if(recType == 1)
%         print([pathAnal 'CorrTActivityPassive'],'-dpdf','-r600');
%     elseif(recType == 2)
%         print([pathAnal 'CorrTActivityActive'],'-dpdf','-r600');
%     elseif(recType == 3)
%         print([pathAnal 'CorrTActivityPassiveNoCue'],'-dpdf','-r600');
%     end
%     
%     % activity similarity over time
%     xlabelTmp = 'SimT. neuronal activity over trials';
% 
%     plotHist(recListActivityCorr.meanSimTActiveTrPerNeu,-0.025:0.0025:0.35,1,...
%         xlabelTmp,titleTmp);
%     if(recType == 1)
%         print([pathAnal 'SimTActivityPassive'],'-dpdf','-r600');
%     elseif(recType == 2)
%         print([pathAnal 'SimTActivityActive'],'-dpdf','-r600');
%     elseif(recType == 3)
%         print([pathAnal 'SimTActivityPassiveNoCue'],'-dpdf','-r600');
%     end

%     % meanFR
%     xlabelTmp = 'Mean firing rate (Hz)';
% 
%     plotHist(recListActivityCorr.meanFRActiveTrPerNeu,...
%         [param.minFR:0.025:param.maxFR],1,...
%         xlabelTmp,titleTmp);
%     if(recType == 1)
%         print([pathAnal 'MeanFRPassive'],'-dpdf','-r600');
%     elseif(recType == 2)
%         print([pathAnal 'MeanFRActive'],'-dpdf','-r600');
%     elseif(recType == 3)
%         print([pathAnal 'MeanFRPassiveNoCue'],'-dpdf','-r600');
%     end
% 
%     % spatial info
%     xlabelTmp = 'Spatial infomation';
%     
%     plotHist(recListActivityCorr.spatialInfoPerNeu,...
%         [min(recListActivityCorr.spatialInfoPerNeu)-0.1:0.025:...
%             max(recListActivityCorr.spatialInfoPerNeu)+0.1],...
%         1,xlabelTmp,titleTmp);
%     if(recType == 1)
%         print([pathAnal 'SpatInfoPassive'],'-dpdf','-r600');
%     elseif(recType == 2)
%         print([pathAnal 'SpatInfoActive'],'-dpdf','-r600');
%     elseif(recType == 3)
%         print([pathAnal 'SpatInfoPassiveNoCue'],'-dpdf','-r600');
%     end
%     
%     % adapt spatial info
%     xlabelTmp = 'Adapt spatial infomation';
%     
%     plotHist(recListActivityCorr.adaptSpatialInfoPerNeu,...
%         [min(recListActivityCorr.adaptSpatialInfoPerNeu)-0.1:0.04:...
%             max(recListActivityCorr.adaptSpatialInfoPerNeu)+0.1],...
%         1,xlabelTmp,titleTmp);
%     if(recType == 1)
%         print([pathAnal 'AdaptSpatInfoPassive'],'-dpdf','-r600');
%     elseif(recType == 2)
%         print([pathAnal 'AdaptSpatInfoActive'],'-dpdf','-r600');
%     elseif(recType == 3)
%         print([pathAnal 'AdaptSpatInfoPassiveNoCue'],'-dpdf','-r600');
%     end
%     
%     % sparsity
%     xlabelTmp = 'Sparsity';
% 
%     plotHist(recListActivityCorr.sparsityPerNeu,...
%         [min(recListActivityCorr.sparsityPerNeu)-0.1:0.07:...
%             max(recListActivityCorr.sparsityPerNeu)+0.1],...
%         1,xlabelTmp,titleTmp);
%     if(recType == 1)
%         print([pathAnal 'SparsityPassive'],'-dpdf','-r600');
%     elseif(recType == 2)
%         print([pathAnal 'SparsityActive'],'-dpdf','-r600');
%     elseif(recType == 3)
%         print([pathAnal 'SparsityPassiveNoCue'],'-dpdf','-r600');
%     end
%     
%     % SNR
%     xlabelTmp = 'SNR';
%     
%     plotHist(recListActivityCorr.SNRPerNeu,...
%         [min(recListActivityCorr.SNRPerNeu)-0.1:0.1:...
%             max(recListActivityCorr.SNRPerNeu)+0.1],...
%         1,xlabelTmp,titleTmp);
%     if(recType == 1)
%         print([pathAnal 'SNRPassive'],'-dpdf','-r600');
%     elseif(recType == 2)
%         print([pathAnal 'SNRActive'],'-dpdf','-r600');
%     elseif(recType == 3)
%         print([pathAnal 'SNRPassiveNoCue'],'-dpdf','-r600');
%     end
    
end

function plotHist(data,centerpoints,ismean,xl,t)
    fig = figure;
    fig.Renderer = 'Painters';
    [counts,centers] = hist(data,centerpoints);
    prob = counts/length(data);
    h = bar(centers,prob,0.9);
    set(h,'EdgeColor',[0.3 0.3 0.3],'FaceColor',[0.5 0.5 0.5]);
    if(ismean == 0)
        meanData = prctile(data,95); 
    else
        meanData = mean(data(~isnan(data))); 
    end
    hold
    h = plot([meanData meanData],[0 max(prob)+0.02],'r:');
    set(h,'LineWidth',1.5);
    set(gca,'FontSize',14,'YLim',[0 max(prob)+0.02]);
    xlabel(xl);
    ylabel('Prob.')
    title(t);
end
