function alignToCue2P(path, fileName, mazeSess)
% align the spikes in each trial to the cue onset of the last trial

    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    
    fileNameCue = [fileName '_alignCue_msess' num2str(mazeSess) '.mat'];  
  
    indRecName = strfind(fileName, '_');
    fileNameBehElec = [fileName(1:indRecName) 'Behav2PDataLFP.mat'];
            
    fullPath = [path fileNameBehElec];
    if(exist(fullPath) == 0)
        disp('The recording Behav2PDataLFP file does not exist.');
        disp(fullPath)
        return;
    end
    load(fullPath,'Laps','Track');
    
    fileNameRun = [fileName '_alignRun_msess' num2str(mazeSess) '.mat'];  
    fullPath = [path fileNameRun];
    if(exist(fullPath) == 0)
        disp('The recording _alignRun_msess file does not exist.');
        disp(fullPath)
        return;
    end
    load(fullPath);
    
    GlobalConst2P;
    
    trialsCue = [];
    indTr = find(Laps.mazeSess == mazeSess);
    for i = 2:length(indTr)
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
            MazeSpeedAccel2P(xMM,sampleFq);
        trialsCue.dFF{i} = Track.dFFSM(indEeg,:);
        trialsCue.dFFGF{i} = Track.dFFGF(indEeg,:); % added on 9/6/2023, gaussian filtered dFF
        trialsCue.F{i} = Track.F(indEeg,:);
        trialsCue.Fneu{i} = Track.Fneu(indEeg,:);
        trialsCue.spikes{i} = Track.spikesSM(indEeg,:);
        
        indLick = Laps.lickLfpInd{indTr(i)} > trialsCue.startLfpInd(i);
        trialsCue.lickLfpInd{i} = Laps.lickLfpInd{indTr(i)}(indLick);
        
        %% added by Yingxue on 9/7/2023
        indEegBef = trialsCue.startLfpInd(i)-nSampBef:trialsCue.startLfpInd(i)-1;
        trialsCue.dFFBef{i} = Track.dFFSM(indEegBef,:);
        trialsCue.dFFBefGF{i} = Track.dFFGF(indEegBef,:); % added on 9/6/2023, gaussian filtered dFF
        trialsCue.FBef{i} = Track.F(indEegBef,:);
        trialsCue.FneuBef{i} = Track.Fneu(indEegBef,:);
        trialsCue.spikesBef{i} = Track.spikesSM(indEegBef,:);
        trialsCue.speed_MMsecBef{i} = Track.speed_MMsec(indEegBef);
        
        licks = [Laps.lickLfpInd{indTr(i-1)}; Laps.lickLfpInd{indTr(i)}];
        indLick = licks >= trialsCue.startLfpInd(i)-nSampBef & ... 
            licks <= trialsCue.startLfpInd(i)-1;
        trialsCue.lickLfpIndBef{i} = licks(indLick);
    end
    
    save([path fileNameCue],'trialsCue','-v7.3');
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
