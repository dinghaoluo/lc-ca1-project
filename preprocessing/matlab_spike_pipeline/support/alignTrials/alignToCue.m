function alignToCue(path, fileName, mazeSess)
% align the spikes in each trial to the cue onset of the last trial

    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    
    fileNameCue = [fileName '_alignCue_msess' num2str(mazeSess) '.mat'];  
  
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
    
    trialsCue = [];
    indTr = find(Laps.mazeSess == mazeSess);
    for i = 1:length(indTr)
        trialsCue.pumpLfpInd{i} = Laps.pumpLfpInd{indTr(i)};
        trialsCue.startLfpInd(i) = Laps.startLfpInd(indTr(i)); 
        trialsCue.endLfpInd(i) = Laps.endLfpInd(indTr(i)); 
        trialsCue.numSamples(i) = trialsCue.endLfpInd(i) - ...
            trialsCue.startLfpInd(i)+1;
        
        indEeg = trialsCue.startLfpInd(i):trialsCue.endLfpInd(i);
        trialsCue.eeg{i} = Track.eeg(indEeg);
        xMM = Track.xMMAll(indEeg);
        xMM = xMMContinuous(xMM);
        trialsCue.xMM{i} = xMM;
        [trialsCue.speed_MMsec{i}, trialsCue.accel_MMsecSq{i}] = ...
            MazeSpeedAccel(xMM,sampleFq);
        trialsCue.thetaPhHilb{i} = Track.thetaPhHilb(indEeg);
        trialsCue.thetaPhLinInterp{i} = Track.thetaPhLinInterp(indEeg);
        
        indSpikes = Spike.res >= trialsCue.startLfpInd(i) & ...
            Spike.res < trialsCue.endLfpInd(i);
        trialsCue.res{i} = Spike.res(indSpikes)-trialsCue.startLfpInd(i)+1;    
        trialsCue.totclu{i} = Spike.totclu(indSpikes);
        trialsCue.thPhaseHilbSpike{i} = Spike.thPhaseHilb(indSpikes);
        trialsCue.thPhaseInterpSpike{i} = Spike.thPhaseInterp(indSpikes);
        trialsCue.xMMSpike{i} = xMM(trialsCue.res{i});
        trialsCue.speed_MMsecSpike{i} = ...
                trialsCue.speed_MMsec{i}(trialsCue.res{i});
        trialsCue.accel_MMsecSqSpike{i} = ...
                trialsCue.accel_MMsecSq{i}(trialsCue.res{i});

        indLick = Laps.lickLfpInd{indTr(i)} > trialsCue.startLfpInd(i);
        trialsCue.lickLfpInd{i} = Laps.lickLfpInd{indTr(i)}(indLick);
        
        %% added 20 Dec 2023, Dinghao, for monitoring pre-cue sequence
        indSpikes = find(Spike.res < trialsCue.startLfpInd(i) & ...
            Spike.res >= trialsCue.startLfpInd(i) - nSampBef);
        trialsCue.resBef{i} = Spike.res(indSpikes)-trialsCue.startLfpInd(i)+1;
        trialsCue.totcluBef{i} = Spike.totclu(indSpikes);
        trialsCue.speed_MMsecSpikeBef{i} = Spike.speed_MMsec(indSpikes);
        trialsCue.thPhaseHilbSpikeBef{i} = Spike.thPhaseHilb(indSpikes);
        trialsCue.thPhaseInterpSpikeBef{i} = Spike.thPhaseInterp(indSpikes); 
        
        %% added on 3/12/2020 by Yingxue, including all the spikes from the 
        % end of last trial to the end of the current trial
        if(i > 1)
            indEndTrackLastTr = find(Track.xMM(Laps.startLfpInd(indTr(i-1)):Laps.endLfpInd(indTr(i-1)))...
                        > Laps.trackLen(indTr(i-1)),1);
            %% added on 7/22/2020 (because sometimes Track.xMM(indEndTrackLastTr) ~= Track.xMMAll(indEndTrackLastTr))
            if(indEndTrackLastTr ~= trialsRun.startLfpInd_LastTr(i) ...
                    - Laps.startLfpInd(indTr(i-1)) + 1)
                indEndTrackLastTr = trialsRun.startLfpInd_LastTr(i) ...
                    - Laps.startLfpInd(indTr(i-1)) + 1;
            end
            indSpikesLastToCurTrEnd = find(Spike.res >= ...
                indEndTrackLastTr+Laps.startLfpInd(indTr(i-1))-1 & ...
                Spike.res < trialsCue.endLfpInd(i)); 
            trialsCue.res_LasttoCurTr{i} = Spike.res(indSpikesLastToCurTrEnd)-...
                trialsCue.startLfpInd(i)+1; 
            trialsCue.totclu_LasttoCurTr{i} = Spike.totclu(indSpikesLastToCurTrEnd);
        end
        
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
