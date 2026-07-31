function activityCorrNeurons(path,fileName,spaceBin,onlyRun,mazeSess)
% calculate the correlation of activity between trials for each neuron

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
        disp(['The firing profile file does not exist. Try to run the',...
                    'function again with fileState = 0.']);
        return;
    end
    load(fullPath,'filteredSpikeArrayNormT');
        
    fileNameInfo = [fileName '_Info.mat'];     
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo(path,fileName);
    end
    load(fullPath,'beh','autoCorr');
    
    indPyrNeu = find(autoCorr.isPyrneuron == 1);
%     ind = abs(cluList.centerMax) > 500;
   
    for sess = 1 % only consider no cue and start cue sessions 
        if(~isempty(strfind(fileName,'A001-20180928')) || ...
                ~isempty(strfind(fileName,'A004-20181027')))
            indTr = find(beh.mazeType == 1);
        elseif(~isempty(strfind(fileName,'A002-20181011')))
            indTr = find(beh.mazeType == 7);
        else
            indTr = intersect(find(beh.mazeSess == mazeSess),beh.indGoodLap);
        end
        
%         indTr = find(beh.mazeType == sess);
%         mazeSessTr = beh.mazeSess(indTr);
%         mazeSess = unique(mazeSessTr);        
%         indTr = indTr(mazeSessTr == mazeSess(1)); 
    
        activityCorr{sess}.indPyrNeu = indPyrNeu;
        activityCorr{sess}.indLaps = indTr;
        activityCorr{sess}.corrMatxPerNeu = [];
        activityCorr{sess}.meanCorrActiveTrPerNeu = zeros(1,length(indPyrNeu));
        activityCorr{sess}.meanCorrAllTrPerNeu = zeros(1,length(indPyrNeu));
        activityCorr{sess}.activeTr = zeros(1,length(indPyrNeu));
        for n = 1:length(indPyrNeu)
            corrMatxPerNeu = ...
                zeros(length(indTr),length(indTr));
            filteredSpikeArrayTmp = zeros(length(indTr),...
                size(filteredSpikeArrayNormT{indTr(1)},2));
            for i = 1:length(indTr)
                disp(['i = ' num2str(i) ' n = ' num2str(n)]);
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
                    countAll = countAll+1;
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
            activityCorr{sess}.corrMatxPerNeu{n} = corrMatxPerNeu;
            activityCorr{sess}.activeTr(n) = activeTr;
            activityCorr{sess}.meanCorrActiveTrPerNeu(n) = ...
                sum(sum(triu(corrMatxPerNeu,1)))/count;
            activityCorr{sess}.meanCorrAllTrPerNeu(n) = ...
                sum(sum(triu(corrMatxPerNeu,1)))/countAll;
        end
        indGoodNeu = activityCorr{sess}.meanCorrActiveTrPerNeu >= 0.01 & ...
            activityCorr{sess}.activeTr/length(indTr) > 0.4;
        activityCorr{sess}.indGoodPyrNeu = indPyrNeu(indGoodNeu);
        
        figure;
        hist(activityCorr{sess}.meanCorrActiveTrPerNeu,-0.025:0.01:0.4);
        xlabel('neuron activity corr. across trials')
        ylabel('count')
        title(fileName(1:13))
        if(sum(activityCorr{sess}.meanCorrActiveTrPerNeu > 0.4) > 0)
            disp('There are neurons with correlation > 0.4');
        end
    end
        
    fileNameSpInfo = [fileName '_NeuActivityCorr_Run' num2str(onlyRun) '.mat'];
    fileNameFull = [path fileNameSpInfo];
    save(fileNameFull,'activityCorr');
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
