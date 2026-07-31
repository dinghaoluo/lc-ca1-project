function alignToCueBeh(path,fileName,mazeSess)
% align each trial to the cue onset

    fileNameCue = [fileName '_alignCue_msess' num2str(mazeSess) '.mat'];   
    
    fileNameBehElec = [fileName '_BehavElectrDataLFP.mat'];
            
    fullPath = [path fileNameBehElec];
    if(exist(fullPath) == 0)
        disp('The recording BehElectDataLFP file does not exist.');
        disp(fullPath)
        return;
    end
    load(fullPath);
    
    GlobalConst;
    
    trialsCue = [];
    indTr = find(Laps.mazeSess == mazeSess);
    for i = 1:length(indTr)
        trialsCue.pumpLfpInd{i} = Laps.pumpLfpInd{indTr(i)};
        trialsCue.startLfpInd(i) = Laps.startLfpInd(indTr(i)); 
        trialsCue.endLfpInd(i) = Laps.endLfpInd(indTr(i)); 
        trialsCue.numSamples(i) = trialsCue.endLfpInd(i) - ...
            trialsCue.startLfpInd(i)+1;
        
        indEeg = trialsCue.startLfpInd(i):trialsCue.endLfpInd(i);
        xMM = Track.xMMAll(indEeg);
        xMM = xMMContinuous(xMM);
        trialsCue.xMM{i} = xMM;
        [trialsCue.speed_MMsec{i}, trialsCue.accel_MMsecSq{i}] = ...
            MazeSpeedAccel(xMM,sampleFq);
        
        indLick = Laps.lickLfpInd{indTr(i)} > trialsCue.startLfpInd(i);
        trialsCue.lickLfpInd{i} = Laps.lickLfpInd{indTr(i)}(indLick) ...
            -trialsCue.startLfpInd(i)+1;
        
    end
    
    save([path fileNameCue],'trialsCue');
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
