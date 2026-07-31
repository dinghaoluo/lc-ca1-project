function CalCCGAlignedRun(path,fileName,onlyRun,mazeSess)
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
    fileNameCCG = [fileName '_CCGAlignedRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    fileNameConSp = [fileName '_ConcatspAlignedRun_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    
    fullPath = [path fileNameConSp];
    if(exist(fullPath) == 0)
        disp(['Concatenating spike file does not exist. Please run function',...
              ' ConcatenateSpikesAlignedRun first']);
        return;
    end
    load(fullPath);
    
    fullPath = [path fileName '_PeakFRAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];    
    if(exist(fullPath) == 0)
        disp('The _PeakFRAligned file does not exist');
        return;
    end
    load(fullPath,'trialNoNonStimGood','trialNoNonStimBad','trialNoStim','trialNoStimCtrl','pulseMeth');

    %%%%%%%%% initialize constants
    GlobalConst;
    
    % calculate CCG for non stimulated good trials'
    disp('Calculate CCG for non stimulated good trials')
    HalfBins = 10000;     % total number of bins on each side
    samplePerBin = 1 * (sampleFq/1000);      
            % number of samples within each bin
    [CCGNonStimGood.ccgVal, CCGNonStimGood.ccgT] = ...
        CCG_wang(spTrainNonStimGood, spCluNonStimGood, ...
                 samplePerBin, HalfBins, sampleFq);
    CCGNonStimGood.ccgT = CCGNonStimGood.ccgT'; 
    CCGNonStimGood.totLfpInd = totLfpIndNonStimGood;
    fullPath = [path fileNameCCG];
    save(fullPath, 'CCGNonStimGood','-v7.3');
    clear CCGNonStimGood
    
    % calculate CCG for non stimulated bad trials'
    disp('Calculate CCG for non stimulated bad trials')
    [CCGNonStimBad.ccgVal, CCGNonStimBad.ccgT] = ...
        CCG_wang(spTrainNonStimBad, spCluNonStimBad, ...
                 samplePerBin, HalfBins, sampleFq);
    CCGNonStimBad.ccgT = CCGNonStimBad.ccgT'; 
    CCGNonStimBad.totLfpInd = totLfpIndNonStimBad;     
    save(fullPath, 'CCGNonStimBad','-v7.3','-append');
    
    %% added by Yingxue on 2/15/2021
    CCGStim = [];
    CCGStimCtrl = [];
    save(fullPath, 'CCGStim','CCGStimCtrl','-v7.3','-append');
    for i = 1:length(pulseMeth)
        % calculate CCG for stimulated trials'
        disp('Calculate CCG for stimulated trials')
        [CCGStim{i}.ccgVal, CCGStim{i}.ccgT] = ...
            CCG_wang(spTrainStim{i}, spCluStim{i}, ...
                     samplePerBin, HalfBins, sampleFq);
        CCGStim{i}.ccgT = CCGStim{i}.ccgT'; 
        CCGStim{i}.totLfpInd = totLfpIndStim{i};
        save(fullPath, 'CCGStim','-v7.3','-append');
        
        % calculate CCG for non stimulated bad trials'
        disp('Calculate CCG for control trials during stimulation')
        [CCGStimCtrl{i}.ccgVal, CCGStimCtrl{i}.ccgT] = ...
            CCG_wang(spTrainStimCtrl{i}, spCluStimCtrl{i}, ...
                     samplePerBin, HalfBins, sampleFq);
        CCGStimCtrl{i}.ccgT = CCGStimCtrl{i}.ccgT'; 
        CCGStimCtrl{i}.totLfpInd = totLfpIndStimCtrl{i};
        save(fullPath, 'CCGStimCtrl','-v7.3','-append');
    end
                   
    clear mydata;

end
