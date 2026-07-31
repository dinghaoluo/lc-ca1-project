function alignToReward2P(path, fileName, mazeSess)
% align the spikes in each trial to the reward time of the last trial

    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
   
    fileNameRew = [fileName '_alignRew_msess' num2str(mazeSess) '.mat'];                    
   
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
    
    trialsRew = [];
    indTr = find(Laps.mazeSess == mazeSess);
    for i = 1:length(indTr)-1
        if(i == 1)
            trialsRew.pumpLfpInd{i} = Laps.pumpLfpInd{indTr(i)};
        end
        trialsRew.pumpLfpInd{i+1} = Laps.pumpLfpInd{indTr(i+1)};
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
            MazeSpeedAccel2P(xMM,sampleFq);
        trialsRew.dFF{i+1} = Track.dFFSM(indEeg,:);
        trialsRew.dFFGF{i+1} = Track.dFFGF(indEeg,:); % added on 9/6/2023, gaussian filtered dFF
        trialsRew.F{i+1} = Track.F(indEeg,:);
        trialsRew.Fneu{i+1} = Track.Fneu(indEeg,:);
        trialsRew.spikes{i+1} = Track.spikesSM(indEeg,:);
        
        licks = [Laps.lickLfpInd{indTr(i)}; Laps.lickLfpInd{indTr(i+1)}];
        indLick = licks > trialsRew.startLfpInd(i+1);
        trialsRew.lickLfpInd{i+1} = licks(indLick);
        
        %% added by Yingxue on 9/7/2023
        indEegBef = trialsRew.startLfpInd(i+1)-nSampBef:trialsRew.startLfpInd(i+1)-1;
        trialsRew.dFFBef{i+1} = Track.dFFSM(indEegBef,:);
        trialsRew.dFFBefGF{i+1} = Track.dFFGF(indEegBef,:); % added on 9/6/2023, gaussian filtered dFF
        trialsRew.FBef{i+1} = Track.F(indEegBef,:);
        trialsRew.FneuBef{i+1} = Track.Fneu(indEegBef,:);
        trialsRew.spikesBef{i+1} = Track.spikesSM(indEegBef,:);
        trialsRew.speed_MMsecBef{i+1} = Track.speed_MMsec(indEegBef);
        
        licks = [Laps.lickLfpInd{indTr(i)}; Laps.lickLfpInd{indTr(i+1)}];
        indLick = licks >= trialsRew.startLfpInd(i+1)-nSampBef & ... 
            licks <= trialsRew.startLfpInd(i+1)-1;
        trialsRew.lickLfpIndBef{i+1} = licks(indLick);
    end
    
    save([path fileNameRew],'trialsRew','-v7.3');
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
