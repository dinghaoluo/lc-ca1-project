function CalCCG_CtrlOnly(path,fileName,onlyRun,mazeSess)
% Calculate CCG for each subsession
% 
% by Yingxue 8/25/2017

    %%%%%%%%% check arguments
    if nargin<2
        disp('At least three arguments are needed for this function.');
        return;
    elseif nargin == 2
        onlyRun = 1;
    elseif nargin > 4
        disp('Too many input arguments');        
        return;
    end
        
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    fileNameCCG = [fileName '_CCG_Ctrl_Run' num2str(onlyRun) '_mazeSess' num2str(mazeSess) '.mat'];
    fileNameConSp = [fileName '_Concatsp_Ctrl_Run' num2str(onlyRun) '_mazeSess' num2str(mazeSess) '.mat'];
    
    fullPath = [path fileNameConSp];
    if(exist(fullPath) == 0)
        disp(['Concatenating spike file does not exist. Please run function',...
              ' ConcatenateSpikes_smTr_CtrlOnly first']);
        return;
    end
    load(fullPath);

   
    %%%%%%%%% initialize constants
    GlobalConst;
    fileNameInfo = [fileName '_Info.mat'];
        
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo_smTr(path,fileName);
    end
    load(fullPath);
    
    % calculate CCG for each subsession
    disp('Calculate CCG for ctrl trials')
    HalfBins = 10000;     % total number of bins on each side
    samplePerBin = 1 * (sampleFq/1000);  % number of samples within each bin 
    [CCGSessCtrl.ccgVal, CCGSessCtrl.ccgT] = ...
            CCG_wang(spTrainCtrl, spCluCtrl, ...
                     samplePerBin, HalfBins, sampleFq);
    CCGSessCtrl.ccgT = CCGSessCtrl.ccgT'; 
    CCGSessCtrl.totLfpInd = totLfpIndCtrl;
    
    fullPath = [path fileNameCCG];
    save(fullPath, 'CCGSessCtrl','-v7.3');
                   
    clear mydata;

end