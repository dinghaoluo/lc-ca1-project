function detectNeuWithPeak(path, fileName, onlyRun, mazeSess)
% detect neurons with significant peaks using shuffling method

     %%%%%%%%% load recording file
    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    fileNamePeakFR = [fileName '_PeakFR_msess' num2str(mazeSess) '_RunOnset' num2str(onlyRun) '.mat'];        
    fileNameConv = [fileName '_convSpikesAligned_msess' num2str(mazeSess) '_BefRun' num2str(onlyRun) '.mat'];
    
    fullPath = [path fileNameConv];
    if(exist(fullPath,'file') == 0)
        disp(['The firing profile file does not exist. Try to run the',...
                    ' "ConvSpikeTrain_AlignedRunOnset" function first.']);
        return;
    end
    load(fullPath);
        
    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The _behPar file does not exist');
        return;
    end
    load(fullPath);
    
    fullPath = [path fileName '_PeakFRAligned_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];    
    if(exist(fullPath) == 0)
        disp('The _PeakFRAligned file does not exist');
        return;
    end
    load(fullPath,'trialNoNonStimGood','trialNoNonStimBad');
            
    % significance probability
    p = [99.9,99];
    
    % median trial length
    medTrLen = prctile(behPar.numSamplesRun,75)/1250;
    indTime = find(timeStepRun/1250 >= medTrLen,1);
    if(isempty(indTime))
        indTime = length(timeStepRun);
    end
    isPeakNeuArr = zeros(length(p),length(filteredSpikeArrayRunOnSet));
    isPeakNeuArrNoStim = zeros(length(p),length(filteredSpikeArrayRunOnSet));
    
    % detect neurons with peaks based on all the trials in the session
    for i = 1:length(filteredSpikeArrayRunOnSet)
        isPeakNeuArr(:,i) = neuPeakDetection(filteredSpikeArrayRunOnSet{i}(:,1:indTime),paramC,i,p);
    end
    
    indPeakNeu = find(isPeakNeuArr(1,:) == 1)
        
    % detect neurons with peaks based on good trials in the session
    trialNoNoStim = [trialNoNonStimGood; trialNoNonStimBad];
    for i = 1:length(filteredSpikeArrayRunOnSet)
        isPeakNeuArrNoStim(:,i) = neuPeakDetection(filteredSpikeArrayRunOnSet{i}(trialNoNoStim,1:indTime),paramC,i,p);
    end
    
    fullpath = [path fileName '_NeuronWPeak_msess' num2str(mazeSess) '_RunOnset' num2str(onlyRun) '.mat'];
    save(fullpath, 'isPeakNeuArr','isPeakNeuArrNoStim');
end