function activityCorrNeuronsBehCluster(path,fileName,spaceBin,onlyRun,mazeSess)
% calculate the correlation of activity between trials for each neuron
% trials are selected based on the behavior clustering

    %%%%%%%%% check arguments
    if nargin<2
        disp('At least two arguments are needed for this function.');
        return;
     elseif nargin == 2
        spaceBin = 2; % cm
        onlyRun = 1;
        mazeSess = 1;
    elseif nargin == 3
        onlyRun = 1;
        mazeSess = 1;
    elseif nargin == 4
        mazeSess = 1;    
    elseif nargin > 5
        disp('Too many input arguments.');
        return;
    end
    
    fullPath = [path fileName '.mat'];
    if(exist(fullPath,'file') == 0)
        disp('recording file does not exist');
    else
        load(fullPath,'cluList');
    end
    
    fileNameConv = [fileName '_convSpikesDist' num2str(spaceBin) ...
                    'mm_Run' num2str(onlyRun) '.mat'];          
    fullPath = [path fileNameConv];
    if(exist(fullPath,'file') == 0)
        if(fileState == 0)
            disp('File does not exist.');
        else
            disp(['The firing profile file does not exist. Try to run the',...
                    'function again with fileState = 0.']);
        end
        return;
    end
    load(fullPath,'filteredSpikeArrayNormT');
    
    fileNameCluster = [fileName '_kmeanBeh.mat'];
    fullPath = [path fileNameCluster];
    if(exist(fullPath) == 0)
        disp('The behavior cluster file does not exist.');
        disp(fullPath)
        return;
    end
    load(fullPath);
    
    fileNameInfo = [fileName '_Info.mat'];     
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo(path,fileName);
    end
    load(fullPath,'autoCorr');
    
    indPyrNeu = find(autoCorr.isPyrneuron == 1);
%     ind = abs(cluList.centerMax) > 500;

    [indCluGood,indCluBad] = HandSelectedCluFromBehKMean(fileName);
    
    indTr = indLaps(kmeanClu == indCluGood);
    
    behBestCluActivityCorr.indPyrNeu = indPyrNeu;
    behBestCluActivityCorr.indLaps = indTr;
    behBestCluActivityCorr.corrMatxPerNeu = [];
    behBestCluActivityCorr.activeTr = zeros(1,length(indPyrNeu));
    behBestCluActivityCorr.meanCorrActiveTrPerNeu = zeros(1,length(indPyrNeu));
    behBestCluActivityCorr.meanCorrAllTrPerNeu = zeros(1,length(indPyrNeu));
    for n = 1:length(indPyrNeu)
        corrMatxPerNeu = ...
            zeros(length(indTr),length(indTr));
        filteredSpikeArrayTmp = zeros(length(indTr),...
            size(filteredSpikeArrayNormT{1},2));
        for i = 1:length(indTr)
            filteredSpikeArrayTmp(i,:) = ...
                filteredSpikeArrayNormT{indTr(i)}(indPyrNeu(n),:);
        end
        count = 0;
        countAll = 0;
        activeTr = 0;
        for i = 1:length(indTr)
            if(sum(filteredSpikeArrayTmp(i,:)) > 0)
                activeTr = activeTr + 1;
            end
            for j = i+1:length(indTr)
                countAll = countAll + 1;
                if(sum(filteredSpikeArrayTmp(i,:)) == 0 || ...
                        sum(filteredSpikeArrayTmp(j,:)) == 0)
                    corrCoef(1,2) = 0;
                else
                    count = count + 1;
                    corrCoef = corrcoef(filteredSpikeArrayTmp(i,:), ...
                            filteredSpikeArrayTmp(j,:));
                    if(isnan(corrCoef(1,2)))
                        corrCoef(1,2) = 0;
                    end
                end
                corrMatxPerNeu(i,j) = corrCoef(1,2);
                corrMatxPerNeu(j,i) = corrMatxPerNeu(i,j);    
            end
        end
        behBestCluActivityCorr.corrMatxPerNeu{n} = corrMatxPerNeu;
        behBestCluActivityCorr.activeTr(n) = activeTr;
        behBestCluActivityCorr.meanCorrActiveTrPerNeu(n) = ...
            sum(sum(triu(corrMatxPerNeu,1)))/count;
        behBestCluActivityCorr.meanCorrAllTrPerNeu(n) = ...
            sum(sum(triu(corrMatxPerNeu,1)))/countAll;
    end
    indGoodNeu = behBestCluActivityCorr.meanCorrActiveTrPerNeu >= 0.01 & ...
            behBestCluActivityCorr.activeTr/length(indTr) > 0.4;
    behBestCluActivityCorr.indGoodPyrNeu = indPyrNeu(indGoodNeu);
        
    indTr = indLaps(kmeanClu == indCluBad);
    
    behBadCluActivityCorr.indPyrNeu = indPyrNeu;
    behBadCluActivityCorr.indLaps = indTr;
    behBadCluActivityCorr.corrMatxPerNeu = [];
    behBadCluActivityCorr.activeTr = zeros(1,length(indPyrNeu));
    behBadCluActivityCorr.meanCorrActiveTrPerNeu = zeros(1,length(indPyrNeu));
    behBadCluActivityCorr.meanCorrAllTrPerNeu = zeros(1,length(indPyrNeu));
    for n = 1:length(indPyrNeu)
        corrMatxPerNeu = ...
            zeros(length(indTr),length(indTr));
        filteredSpikeArrayTmp = zeros(length(indTr),...
            size(filteredSpikeArrayNormT{1},2));
        for i = 1:length(indTr)
            filteredSpikeArrayTmp(i,:) = ...
                filteredSpikeArrayNormT{indTr(i)}(indPyrNeu(n),:);
        end
        count = 0;
        countAll = 0;
        activeTr = 0;
        for i = 1:length(indTr)
            if(sum(filteredSpikeArrayTmp(i,:)) > 0)
                activeTr = activeTr + 1;
            end
            for j = i+1:length(indTr)
                countAll = countAll + 1;
                if(sum(filteredSpikeArrayTmp(i,:)) == 0 || ...
                        sum(filteredSpikeArrayTmp(j,:)) == 0)
                    corrCoef(1,2) = 0;
                else
                    count = count + 1;
                    corrCoef = corrcoef(filteredSpikeArrayTmp(i,:), ...
                            filteredSpikeArrayTmp(j,:));
                    if(isnan(corrCoef(1,2)))
                        corrCoef(1,2) = 0;
                    end
                end
                corrMatxPerNeu(i,j) = corrCoef(1,2);
                corrMatxPerNeu(j,i) = corrMatxPerNeu(i,j);    
            end
        end
        behBadCluActivityCorr.corrMatxPerNeu{n} = corrMatxPerNeu;
        behBadCluActivityCorr.activeTr(n) = activeTr;
        behBadCluActivityCorr.meanCorrActiveTrPerNeu(n) = ...
            sum(sum(triu(corrMatxPerNeu,1)))/count;
        behBadCluActivityCorr.meanCorrAllTrPerNeu(n) = ...
            sum(sum(triu(corrMatxPerNeu,1)))/countAll;
    end
    
    figure;
    hist(behBestCluActivityCorr.meanCorrActiveTrPerNeu,-0.025:0.01:0.4);
    xlabel('neuron activity corr. across trials')
    ylabel('count')
    title(fileName(1:13))
    if(sum(behBestCluActivityCorr.meanCorrActiveTrPerNeu > 0.4) > 0)
        disp('There are neurons with correlation > 0.4');
    end
    
    indSelPyr = [];
    if(~isempty(strfind(fileName,'A011-20190218')))
        indSelPyr = [3 21 53 54 59 80 95 129 134 140 144 165 172 177 186 206 216 226];
    elseif(~isempty(strfind(fileName,'A011-20190219')))
        indSelPyr = [9 16 19 90 93 130 140 143 182];
    elseif(~isempty(strfind(fileName,'A012-20190224')))
        indSelPyr = [13 15 16 30 31 38 39 55 63 65 66 85 86 87 88 98 103 108 134];
    elseif(~isempty(strfind(fileName,'A012-20190221')))
        indSelPyr = [15 30 31 33 35 40 42 43 52 66 69 70];
    end
    
    if(~isempty(indSelPyr))
        [~,ia,ib] = intersect(indPyrNeu,indSelPyr);
        figure;
        plot(behBestCluActivityCorr.meanCorrActiveTrPerNeu(ia),...
                behBadCluActivityCorr.meanCorrActiveTrPerNeu(ia),'bo');
        hold on;
        maxMeanCorr =  max([behBestCluActivityCorr.meanCorrActiveTrPerNeu(ia),...
                behBadCluActivityCorr.meanCorrActiveTrPerNeu(ia)]);
        plot([0,maxMeanCorr],[0, maxMeanCorr],'r');
        set(gca,'XLim',[0 maxMeanCorr]);
        title(fileName(1:13))
        xlabel('Good behavior')
        ylabel('Bad behavior')
    end
    
    fileNameSpInfo = [fileName '_ActivityCorrKluBeh_Run' num2str(onlyRun) '.mat'];
    fileNameFull = [path fileNameSpInfo];
    save(fileNameFull,'behBestCluActivityCorr','behBadCluActivityCorr')
end

% A012-20190224
% indSelPyr = [13 15 16 30 31 38 39 55 63 65 66 85 86 87 88 98 103 108 134];
% A011-20190219
% indSelPyr = [9 16 19 90 93 130 140 143 182];
% A011-20190218
% indSelPyr = [3 21 53 54 59 80 95 129 134 140 144 165 172 177 186 206 216
% 226];
% A012-20190221
% indSelPyr = [15 30 31 33 35 40 42 43 52 66 69 70];
