function averageRunningProfile(path, fileName)
    
    fullPath = [path fileName '.mat']; 
    if(exist(fullPath) == 0)
        disp('The recording file does not exist');
        return;
    end
    load(fullPath,'cluList','lap');
    neuronNo = 1:length(cluList.all);
    trialNo = 1:length(lap.trackLen);
    
    fullPath = [path fileName '_alignRun.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to run file does not exist');
        return;
    end
    load(fullPath,'trialsRun');
    
    GlobalConst;
    
    figure
    hold on;
    numTr = 0;
    trialLen = trialLenT*sampleFq;
    speedAvg = zeros(1,trialLenT*sampleFq); 
    for i = trialNo
        if(~isempty(trialsRun.speed_MMsec{i}))
            speedSM = smooth(trialsRun.speed_MMsec{i},smoothSpan);
            lenTr = length(speedSM);
            if(lenTr > trialLen)
                speedAvg = speedAvg + speedSM(1:trialLen)';
            else
                speedAvg(1:lenTr) = speedAvg(1:lenTr) + speedSM';
            end
            numTr = numTr + 1;
            plot([1:length(speedSM)]/sampleFq,speedSM);
        end
    end  
    speedAvg = speedAvg/numTr;
    h = plot([1:length(speedAvg)]/sampleFq,speedAvg,'k-');
    set(h,'LineWidth',8);
    xlabel('Time (s)')
    ylabel('Speed (mm/s)')
end
