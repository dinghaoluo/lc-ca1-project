function alignToRewardBeh(path, fileName, mazeSess)
% align each trial to the reward time of the last trial

    fileNameRew = [fileName '_alignRew_msess' num2str(mazeSess) '.mat'];                    
 
    fileNameBehElec = [fileName '_BehavElectrDataLFP.mat'];            
    fullPath = [path fileNameBehElec];
    if(exist(fullPath) == 0)
        disp('The recording BehElectDataLFP file does not exist.');
        disp(fullPath)
        return;
    end
    load(fullPath);
    
    GlobalConst;
    
    trialsRew = [];
    indTr = find(Laps.mazeSess == mazeSess);
    for i = 1:length(indTr)-1
        trialsRew.pumpLfpInd{i} = Laps.pumpLfpInd{indTr(i)};
        if(~isempty(trialsRew.pumpLfpInd{i}))
            trialsRew.goodTrial(i+1) = 1;
            trialsRew.startLfpInd(i+1) = Laps.pumpLfpInd{indTr(i)}(1); 
        else
            trialsRew.goodTrial(i+1) = -1; % trials without reward
            trialsRew.startLfpInd(i+1) = Laps.startLfpInd(indTr(i+1)); 
        end
        trialsRew.endLfpInd(i+1) = Laps.endLfpInd(indTr(i+1)); 
        trialsRew.numSamples(i+1) = trialsRew.endLfpInd(i+1) - ...
            trialsRew.startLfpInd(i+1)+1;

        indEeg = trialsRew.startLfpInd(i+1):trialsRew.endLfpInd(i+1);
        xMM = Track.xMMAll(indEeg);
        xMM = xMMContinuous(xMM);
        trialsRew.xMM{i+1} = xMM;
        [trialsRew.speed_MMsec{i+1}, trialsRew.accel_MMsecSq{i+1}] = ...
            MazeSpeedAccel(xMM,sampleFq);
        
        licks = [Laps.lickLfpInd{indTr(i)}; Laps.lickLfpInd{indTr(i+1)}];
        indLick = licks > trialsRew.startLfpInd(i+1);
        trialsRew.lickLfpInd{i+1} = licks(indLick)-trialsRew.startLfpInd(i+1)+1;
    end
    
    save([path fileNameRew],'trialsRew');
end

function xMM1 = xMMContinuous(xMM)
    indResetXMM = find(diff(xMM) < 0);
    indResetXMM = [indResetXMM+1; length(xMM)+1];
    for n = 1:length(indResetXMM)-1
        xMM(indResetXMM(n):indResetXMM(n+1)-1) = ...
            xMM(indResetXMM(n):indResetXMM(n+1)-1) + xMM(indResetXMM(n)-1);
    end
    xMM1 = xMM-xMM(1);  
end
