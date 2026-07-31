function getTimePerDistBin(path,fileName,onlyRun,mazeSess)
% get the time spend on each distance bin

    fullPath = [path fileName '_alignRun_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to run file does not exist');
        return;
    end
    load(fullPath,'trialsRun');
    
    fileNameTpD = [fileName '_timePerDistBin_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'];
    
    GlobalConst2P;
    
    trialNo = length(trialsRun.xMM);
    for j = 1:trialNo
        if(~isempty(trialsRun.xMM{j}))
            spaceSteps = 0:spaceMergeBin:max(trialsRun.xMM{j});
            timePerDistBinRun.spaceSteps{j} = spaceSteps;
            if(onlyRun == 1)
                ind = trialsRun.speed_MMsec{j} > minSpeed;
                spikes = trialsRun.spikes{j}(ind,:);
                xMM = trialsRun.xMM{j}(ind);
            else
                spikes = trialsRun.spikes{j};
                xMM = trialsRun.xMM{j};
            end
            [timePerDistBinRun.spikesHist{j},timePerDistBinRun.hist{j}] = ...
                spikeTime2Dist(spikes, xMM, spaceSteps);
        end
    end
    
    save([path fileNameTpD], 'timePerDistBinRun','-v7.3');   
    clear trialsRun
    
    fullPath = [path fileName '_alignRew_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to reward file does not exist');
        return;
    end
    load(fullPath,'trialsRew');
    
    for j = 1:trialNo
        if(~isempty(trialsRew.xMM{j}))
            spaceSteps = 0:spaceMergeBin:max(trialsRew.xMM{j});
            timePerDistBinRew.spaceSteps{j} = spaceSteps;
            if(onlyRun == 1)
                ind = trialsRew.speed_MMsec{j} > minSpeed;
                spikes = trialsRew.spikes{j}(ind,:);
                xMM = trialsRew.xMM{j}(ind);
            else
                spikes = trialsRew.spikes{j};
                xMM = trialsRew.xMM{j};
            end
            [timePerDistBinRew.spikesHist{j},timePerDistBinRew.hist{j}] = ...
                spikeTime2Dist(spikes, xMM, spaceSteps);
        end
    end
    
    save([path fileNameTpD], 'timePerDistBinRew','-append');   
    clear trialsRew
    
    fullPath = [path fileName '_alignCue_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to cue file does not exist');
        return;
    end
    load(fullPath,'trialsCue');
    
    for j = 1:trialNo
        if(~isempty(trialsCue.xMM{j}))
            spaceSteps = 0:spaceMergeBin:max(trialsCue.xMM{j});
            timePerDistBinCue.spaceSteps{j} = spaceSteps;
            if(onlyRun == 1)
                ind = trialsCue.speed_MMsec{j} > minSpeed;
                spikes = trialsCue.spikes{j}(ind,:);
                xMM = trialsCue.xMM{j}(ind);
            else
                spikes = trialsCue.spikes{j};
                xMM = trialsCue.xMM{j};
            end
            [timePerDistBinCue.spikesHist{j},timePerDistBinCue.hist{j}] = ...
                spikeTime2Dist(spikes, xMM, spaceSteps);
        end
    end
    
    save([path fileNameTpD], 'timePerDistBinCue','-append');  
    
    clear all;
end
    
% function timePerBin = spikeTime2Dist(xMM, spaceSteps)
%                         
%     numBins = length(spaceSteps);
%     step = spaceSteps(2) - spaceSteps(1);
%     timePerBin = zeros(1,numBins);
%     for i = 1:numBins
%         ind = find(xMM >= spaceSteps(i)-step/2 & xMM < spaceSteps(i)+step/2);
%         if(~isempty(ind))
%             timePerBin(i) = length(ind);
%         else
%             timePerBin(i) = 1;
%         end
%     end 
% end

function [spikesPerBin,timePerBin] = ...
                            spikeTime2Dist(spikes,dist,spaceSteps)
                        
    numBins = length(spaceSteps);
    step = spaceSteps(2) - spaceSteps(1);

    timePerBin = hist(dist(dist<=spaceSteps(end)+step/2),spaceSteps);

    numNeurons = size(spikes,2);
    spikesPerBin = zeros(numNeurons,numBins);
    accumTime = 1;
    for i = 1:numBins
        if(timePerBin(i) == 0)
            continue;
        end
        spikePerBinTmp = spikes(accumTime:accumTime+timePerBin(i)-1,:)';
        spikesPerBin(:,i) = sum(spikePerBinTmp,2);
        accumTime = accumTime + timePerBin(i);
    end
    
end
