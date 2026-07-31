function getTimePerDistBin(path,fileName,onlyRun,mazeSess)
% get the time spend on each distance bin

    fullPath = [path fileName '_alignRun_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to run file does not exist');
        return;
    end
    load(fullPath,'trialsRun');
    
    fullPath = [path fileName '_alignRew_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to reward file does not exist');
        return;
    end
    load(fullPath,'trialsRew');
    
    fullPath = [path fileName '_alignCue_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to cue file does not exist');
        return;
    end
    load(fullPath,'trialsCue');
    
    fileNameTpD = [fileName '_timePerDistBin_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    
    GlobalConst;
    
    trialNo = length(trialsRun.xMM);
    for j = 1:trialNo
        if(~isempty(trialsRun.xMM{j}))
            spaceSteps = 0:spaceMergeBin:max(trialsRun.xMM{j});
            timePerDistBinRun.spaceSteps{j} = spaceSteps;
            if(onlyRun == 1)
                xMM = trialsRun.xMM{j}(trialsRun.speed_MMsec{j} > minSpeed); % changed from >= to > on 11/14/2021
            else
                xMM = trialsRun.xMM{j};
            end
            timePerDistBinRun.hist{j} = spikeTime2Dist(xMM, spaceSteps);
        end
        
        if(~isempty(trialsRew.xMM{j}))
            spaceSteps = 0:spaceMergeBin:max(trialsRew.xMM{j});
            timePerDistBinRew.spaceSteps{j} = spaceSteps;
            if(onlyRun == 1)
                xMM = trialsRew.xMM{j}(trialsRew.speed_MMsec{j} > minSpeed); % changed from >= to > on 11/14/2021
            else
                xMM = trialsRew.xMM{j};
            end
            timePerDistBinRew.hist{j} = spikeTime2Dist(xMM, spaceSteps);
        end
        
        if(~isempty(trialsCue.xMM{j}))
            spaceSteps = 0:spaceMergeBin:max(trialsCue.xMM{j});
            timePerDistBinCue.spaceSteps{j} = spaceSteps;
            if(onlyRun == 1)
                xMM = trialsCue.xMM{j}(trialsCue.speed_MMsec{j} > minSpeed); % changed from >= to > on 11/14/2021
            else
                xMM = trialsCue.xMM{j};
            end
            timePerDistBinCue.hist{j} = spikeTime2Dist(xMM, spaceSteps);
        end
    end
    
    save([path fileNameTpD], 'timePerDistBinRun', 'timePerDistBinRew',...
        'timePerDistBinCue');    
end
    
function timePerBin = spikeTime2Dist(xMM, spaceSteps)
                        
    numBins = length(spaceSteps);
    step = spaceSteps(2) - spaceSteps(1);
    timePerBin = zeros(1,numBins);
    for i = 1:numBins
        ind = find(xMM >= spaceSteps(i)-step/2 & xMM < spaceSteps(i)+step/2);
        if(~isempty(ind))
            timePerBin(i) = length(ind);
        else
            timePerBin(i) = 1;
        end
    end 
end
