function CalCCG_GoodTr(path,fileName,onlyRun)
% Calculate CCG for each subsession
% 
% by Yingxue 8/25/2017

    %%%%%%%%% check arguments
    if nargin<2
        disp('At least three arguments are needed for this function.');
        return;
    elseif nargin == 2
        onlyRun = 1;
    elseif nargin > 3
        disp('Too many input arguments');        
        return;
    end
        
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    fileNameCCG = [fileName '_CCG_GoodTr_Run' num2str(onlyRun) '.mat'];
    fileNameConSp = [fileName '_Concatsp_GoodTr_Run' num2str(onlyRun) '.mat'];
    
%     if(onlyRun == 1)
        fullPath = [path fileNameConSp];
        if(exist(fullPath) == 0)
            disp(['Concatenating spike file does not exist. Please run function',...
                  ' ConcatenateSpikes first']);
            return;
        end
        load(fullPath);
%     else
%         disp('CCG is only calculated for the running period');
%     end
   
    %%%%%%%%% initialize constants
    GlobalConst;
    fileNameInfo = [fileName '_Info.mat'];
        
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo_smTr(path,fileName);
    end
    load(fullPath);
    mazeSess = beh.mazeSessAll;
    
    % calculate CCG for each subsession
    disp('Calculate CCG for each subsession')
    CCGSessGoodTr = cell(1,length(mazeSess));
    CCGSessBadTr = cell(1,length(mazeSess));
    CCGSessOKTr = cell(1,length(mazeSess));
    HalfBins = 10000;     % total number of bins on each side
    samplePerBin = 1 * (sampleFq/1000);  % number of samples within each bin    
    for i = 1:length(mazeSess)
        disp(['Session ' num2str(i)]);
        
        disp('Good trials');
        [CCGSessGoodTr{i}.ccgVal, CCGSessGoodTr{i}.ccgT] = ...
            CCG_wang(spTrainGoodTr{i}, spCluGoodTr{i}, ...
                     samplePerBin, HalfBins, sampleFq);
        CCGSessGoodTr{i}.ccgT = CCGSessGoodTr{i}.ccgT'; 
        CCGSessGoodTr{i}.totLfpInd = totLfpIndGoodTr{i};
        
        disp('Bad trials');
        [CCGSessBadTr{i}.ccgVal, CCGSessBadTr{i}.ccgT] = ...
            CCG_wang(spTrainBadTr{i}, spCluBadTr{i}, ...
                     samplePerBin, HalfBins, sampleFq);
        CCGSessBadTr{i}.ccgT = CCGSessBadTr{i}.ccgT'; 
        CCGSessBadTr{i}.totLfpInd = totLfpIndBadTr{i};
        
        disp('OK trials');
        [CCGSessOKTr{i}.ccgVal, CCGSessOKTr{i}.ccgT] = ...
            CCG_wang(spTrainOKTr{i}, spCluOKTr{i}, ...
                     samplePerBin, HalfBins, sampleFq);
        CCGSessOKTr{i}.ccgT = CCGSessOKTr{i}.ccgT'; 
        CCGSessOKTr{i}.totLfpInd = totLfpIndOKTr{i};
    end
    
    fullPath = [path fileNameCCG];
    save(fullPath, 'CCGSessGoodTr','CCGSessBadTr','CCGSessOKTr','-v7.3');
                   
    clear mydata;

end
