function alignToReward1(path, fileName, mazeSess)
% align the spikes in each trial to the reward time of the last trial
% added the spikes before reward on 3/3/2022, as compared to alignedReward

    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
   
    fileNameRew = [fileName '_alignRew_msess' num2str(mazeSess) '.mat'];                    
   
    indRecName = strfind(fileName, '_');
    fileNameBehElec = [fileName(1:indRecName) 'BehavElectrDataLFP.mat'];
            
    fullPath = [path fileNameBehElec];
    if(exist(fullPath) == 0)
        disp('The recording BehElectDataLFP file does not exist.');
        disp(fullPath)
        return;
    end
    load(fullPath);
    
    fileNameRun = [fileName '_alignRun_msess' num2str(mazeSess) '.mat'];  
    fullPath = [path fileNameRun];
    if(exist(fullPath) == 0)
        disp('The recording _alignRun_msess file does not exist.');
        disp(fullPath)
        return;
    end
    load(fullPath);
    
    GlobalConst;
    
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
        trialsRew.eeg{i+1} = Track.eeg(indEeg);
        xMM = Track.xMMAll(indEeg);
        xMM = xMMContinuous(xMM);
        trialsRew.xMM{i+1} = xMM;
        [trialsRew.speed_MMsec{i+1}, trialsRew.accel_MMsecSq{i+1}] = ...
            MazeSpeedAccel(xMM,sampleFq);
        trialsRew.thetaPhHilb{i+1} = Track.thetaPhHilb(indEeg);
        trialsRew.thetaPhLinInterp{i+1} = Track.thetaPhLinInterp(indEeg);

        indSpikes = Spike.res >= trialsRew.startLfpInd(i+1) & ...
            Spike.res < trialsRew.endLfpInd(i+1);
        trialsRew.res{i+1} = Spike.res(indSpikes) - trialsRew.startLfpInd(i+1) + 1;       
        trialsRew.totclu{i+1} = Spike.totclu(indSpikes);
        trialsRew.thPhaseHilbSpike{i+1} = Spike.thPhaseHilb(indSpikes);
        trialsRew.thPhaseInterpSpike{i+1} = Spike.thPhaseInterp(indSpikes);            
        trialsRew.xMMSpike{i+1} = xMM(trialsRew.res{i+1});
        trialsRew.speed_MMsecSpike{i+1} = ...
            trialsRew.speed_MMsec{i+1}(trialsRew.res{i+1});
        trialsRew.accel_MMsecSqSpike{i+1} = ...
            trialsRew.accel_MMsecSq{i+1}(trialsRew.res{i+1});

        licks = [Laps.lickLfpInd{indTr(i)}; Laps.lickLfpInd{indTr(i+1)}];
        indLick = licks > trialsRew.startLfpInd(i+1);
        trialsRew.lickLfpInd{i+1} = licks(indLick);
        
        % added by Dinghao on 3/3/2022, changed by Yingxue on 3/4/2022
        indSpikes = find(Spike.res < trialsRew.startLfpInd(i+1) & ...
            Spike.res >= trialsRew.startLfpInd(i+1) - nSampBefRew);
        trialsRew.resBef{i+1} = Spike.res(indSpikes)-trialsRew.startLfpInd(i+1)+1;
        trialsRew.totcluBef{i+1} = Spike.totclu(indSpikes);
        trialsRew.speed_MMsecSpikeBef{i+1} = Spike.speed_MMsec(indSpikes);
        trialsRew.thPhaseHilbSpikeBef{i+1} = Spike.thPhaseHilb(indSpikes);
        trialsRew.thPhaseInterpSpikeBef{i+1} = Spike.thPhaseInterp(indSpikes); 
        
        %% added on 3/12/2020 by Yingxue, including all the spikes from the 
        % end of last trial to the end of the current trial
        if(i > 0)
            if(i == 34)
                a = 1;
            end
            trackMM = Track.xMM(Laps.startLfpInd(indTr(i)):Laps.endLfpInd(indTr(i)));
            indEndTrackLastTr = find(trackMM > Laps.trackLen(indTr(i)),1);
            %% added on 7/22/2020 (because sometimes Track.xMM(indEndTrackLastTr) ~= Track.xMMAll(indEndTrackLastTr))
            if(indEndTrackLastTr ~= trialsRun.startLfpInd_LastTr(i+1) ...
                    - Laps.startLfpInd(indTr(i)) + 1)
                indEndTrackLastTr = trialsRun.startLfpInd_LastTr(i+1) ...
                    - Laps.startLfpInd(indTr(i)) + 1;
            end
                
            indSpikesLastToCurTrEnd = find(Spike.res >= ...
                indEndTrackLastTr+Laps.startLfpInd(indTr(i))-1 & ...
                Spike.res < trialsRew.endLfpInd(i+1)); 
            trialsRew.res_LasttoCurTr{i+1} = Spike.res(indSpikesLastToCurTrEnd)-...
                trialsRew.startLfpInd(i+1)+1; 
            trialsRew.totclu_LasttoCurTr{i+1} = Spike.totclu(indSpikesLastToCurTrEnd);
        end

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
