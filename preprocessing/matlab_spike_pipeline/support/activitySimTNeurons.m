function activitySimTNeurons(path,fileName,timeBin,onlyRun,mazeSess)
% calculate the correlation of activity between trials for each neuron

    %%%%%%%%% check arguments
    if nargin<2
        disp('At least two arguments are needed for this function.');
        return;
     elseif nargin == 2
        timeBin = 0.096; % cm
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
    
    timebinStr = num2str(timeBin*1000);
    ind = strfind(timebinStr,'.');
    if(~isempty(ind))
        timebinStr(ind) = 'p';
    end
    fileNameConv = [fileName '_convSpikesTime' timebinStr ...
                    'ms_Run' num2str(onlyRun) '.mat'];          
    fullPath = [path fileNameConv];
    if(exist(fullPath,'file') == 0)
        disp(['The firing profile file does not exist. Try to run the',...
                    'function again with fileState = 0.']);
        return;
    end
    load(fullPath,'filteredSpikeArray');
        
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
        
        [spikeTrainSimT,nonZeroTr] = calSpikeTrainSimT(length(indPyrNeu),...
            filteredSpikeArray(indPyrNeu),indTr);

    end
        
    fileNameSpInfo = [fileName '_NeuActivitySimT_Run' num2str(onlyRun) '.mat'];
    fileNameFull = [path fileNameSpInfo];
    save(fileNameFull,'spikeTrainSimT','nonZeroTr');
end

function [spikeTrainSimT,nonZeroTr] = calSpikeTrainSimT(neuronNo,filteredSpikeArray,trials)
    parfor n = 1:neuronNo
        spikeArr = filteredSpikeArray{n}(trials,:);
        nonZeroTr{n} = sum(spikeArr') > 0;
        spikeTrainSimTTmp = pdist(spikeArr,'cosine');
        spikeTrainSimTTmp = squareform(spikeTrainSimTTmp);
%         for i = 1:trialNo
%             if(sum(spikeArr(i,:)) ~= 0)
%                 for j = i+1:trialNo
%                     if(sum(spikeArr(j,:)) ~= 0)
%                         spikesCorrTmp(i,j) = ...
%                             corr(spikeArr(i,:)',spikeArr(j,:)','Type','Spearman');
%                     end
%                 end
%             end            
%         end
        spikeTrainSimT{n} = spikeTrainSimTTmp;
    end
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
