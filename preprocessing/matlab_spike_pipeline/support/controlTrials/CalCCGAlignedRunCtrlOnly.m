function CalCCGAlignedRunCtrlOnly(path,fileName,onlyRun,mazeSess)
% Calculate CCG for each subsession
% 
% by Yingxue 8/25/2017

    %%%%%%%%% check arguments
    if nargin<3
        disp('At least four arguments are needed for this function.');
        return;
    elseif nargin > 4
        disp('Too many input arguments');        
        return;
    end
        
    %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    fileNameCCG = [fileName '_CCGAlignedRunCtrl_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    fileNameConSp = [fileName '_ConcatspAlignedRunCtrl_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    
    fullPath = [path fileNameConSp];
    if(exist(fullPath) == 0)
        disp(['Concatenating spike file does not exist. Please run function',...
              ' ConcatenateSpikesAlignedRunCtrl first']);
        return;
    end
    load(fullPath);
    
    %%%%%%%%% initialize constants
    GlobalConst;
    
    % calculate CCG for non stimulated trials'
    disp('Calculate CCG for non stimulated trials')
    HalfBins = 10000;     % total number of bins on each side
    samplePerBin = 1 * (sampleFq/1000);      
            % number of samples within each bin
    [CCGNonStim.ccgVal, CCGNonStim.ccgT] = ...
        CCG_wang(spTrainNonStim, spCluNonStim, ...
                 samplePerBin, HalfBins, sampleFq);
    CCGNonStim.ccgT = CCGNonStim.ccgT'; 
    CCGNonStim.totLfpInd = totLfpIndNonStim;
    fullPath = [path fileNameCCG];
    save(fullPath, 'CCGNonStim','-v7.3');
    clear CCGNonStim
        
    clear mydata;

end
