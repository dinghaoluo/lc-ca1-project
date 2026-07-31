function alignToCueOff(path, fileName, mazeSess)
% align the spikes in each trial to the running onset

    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    
    fileNameCueOff = [fileName '_alignCueOff_msess' num2str(mazeSess) '.mat'];  
    
    fileNameTheta = [fileName '_thetaPower.mat'];  
       
    fullPath = [path fileNameTheta];
    if(exist(fullPath) == 0)
        disp('The _thetaPower file does not exist.');
        disp(fullPath)
        return;
    end
    load(fullPath);
    
    indRecName = strfind(fileName, '_');
    fileNameBehElec = [fileName(1:indRecName) 'BehavElectrDataLFP.mat'];
            
    fullPath = [path fileNameBehElec];
    if(exist(fullPath) == 0)
        disp('The recording BehElectDataLFP file does not exist.');
        disp(fullPath)
        return;
    end
    load(fullPath,'Laps','Track','Spike');
    
    fileNameRun = [fileName '_alignRun_msess' num2str(mazeSess) '.mat'];  
    fullPath = [path fileNameRun];
    if(exist(fullPath) == 0)
        disp('The recording _alignRun_msess file does not exist.');
        disp(fullPath)
        return;
    end
    load(fullPath);
    
    GlobalConst;
    
    minSpeed1 = 100;
    minSpeed = 10;
    smoothSpan = 100;
    befTime = 3; % take 3 sec before the cue offset
    trialsCueOff = [];
    indTr = find(Laps.mazeSess == mazeSess);
    for i = 1:length(indTr)
        trialsCueOff.pumpLfpInd{i} = Laps.pumpLfpInd{indTr(i)};
        if(i == 1)
            indStart = -1;
            trialsCueOff.goodTrial(i) = indStart;
            continue;
        end
        if(length(Laps.movieOnLfpInd{indTr(i)}) == 3)
            indCueOff = Laps.movieOnLfpInd{indTr(i)}(2);
        elseif(isfield(Laps,'movieLocation') && length(Laps.movieOnLfpInd{indTr(i)}) == 4)
            indCueOff = Laps.movieOnLfpInd{indTr(i)}(3);    
        else
            disp(['trial ' num2str(indTr(i)) ' does not have a cue off signal']);
            trialsCueOff.startLfpInd(i) = -1;
            continue;
        end
            
        xMM = Track.xMMAll(indCueOff-befTime*sampleFq:Laps.endLfpInd(indTr(i)));    
        xMM1 = xMMContinuous(xMM);
        [speed,ac] = MazeSpeedAccel(xMM1,sampleFq);
        speedSM = smooth(speed,smoothSpan);
        indTrueStart = befTime*sampleFq;
        
        trialsCueOff.startLfpInd(i) = indCueOff;

        trialsCueOff.endLfpInd(i) = Laps.endLfpInd(indTr(i)); 
        trialsCueOff.numSamples(i) = trialsCueOff.endLfpInd(i) - ...
            trialsCueOff.startLfpInd(i)+1;
        indSpikes = find(Spike.res >= trialsCueOff.startLfpInd(i) & ...
            Spike.res < trialsCueOff.endLfpInd(i));
        trialsCueOff.res{i} = Spike.res(indSpikes)-trialsCueOff.startLfpInd(i)+1;
        trialsCueOff.totclu{i} = Spike.totclu(indSpikes);
        trialsCueOff.thPhaseHilbSpike{i} = Spike.thPhaseHilb(indSpikes);
        trialsCueOff.thPhaseInterpSpike{i} = Spike.thPhaseInterp(indSpikes);
        indXMMSpikes = trialsCueOff.res{i};
        trialsCueOff.xMMSpike{i} = xMM1(indXMMSpikes+indTrueStart+1)-xMM1(indTrueStart+1);
        trialsCueOff.speed_MMsecSpike{i} = speed(indXMMSpikes+indTrueStart+1);
        trialsCueOff.accel_MMsecSqSpike{i} = ac(indXMMSpikes+indTrueStart+1);  
        indEeg = trialsCueOff.startLfpInd(i):trialsCueOff.endLfpInd(i);
        trialsCueOff.eeg{i} = Track.eeg(indEeg);
        trialsCueOff.xMM{i} = xMM1(indTrueStart+1:end)-xMM1(indTrueStart+1);
        trialsCueOff.speed_MMsec{i} = speed(indTrueStart+1:end);
        trialsCueOff.accel_MMsecSq{i} = ac(indTrueStart+1:end);
        trialsCueOff.thetaPhHilb{i} = Track.thetaPhHilb(indEeg);
        trialsCueOff.thetaPhLinInterp{i} = Track.thetaPhLinInterp(indEeg);
        trialsCueOff.startCueTime(i) = Laps.startLfpInd(indTr(i));
        
        indEegBef = trialsCueOff.startLfpInd(i)-befTime*sampleFq:trialsCueOff.startLfpInd(i)-1;
        trialsCueOff.eegBef{i} = Track.eeg(indEegBef);
        trialsCueOff.thetaPhHilbBef{i} = Track.thetaPhHilb(indEegBef);
        trialsCueOff.thetaPhLinInterpBef{i} = Track.thetaPhLinInterp(indEegBef);
        
        trialsCueOff.ThetaAmp{i} = thetaPower.ThetaAmp(indEeg);
        trialsCueOff.ThetaFreq{i} = thetaPower.ThetaFreq(indEeg);
        
        trialsCueOff.ThetaAmpBef{i} =  thetaPower.ThetaAmp(indEegBef);
        trialsCueOff.ThetaFreqBef{i} = thetaPower.ThetaFreq(indEegBef);
        
        trialsCueOff.speed_MMsecBef{i} = Track.speed_MMsec(indEegBef);
        
        licks = Laps.lickLfpInd{indTr(i)};
        indLick = licks > trialsCueOff.startLfpInd(i);
        trialsCueOff.lickLfpInd{i} = licks(indLick);
        
        indSpikes = find(Spike.res < trialsCueOff.startLfpInd(i) & ...
            Spike.res >= trialsCueOff.startLfpInd(i) - befTime*sampleFq);
        trialsCueOff.resBef{i} = Spike.res(indSpikes)-trialsCueOff.startLfpInd(i)+1;
        trialsCueOff.totcluBef{i} = Spike.totclu(indSpikes);
        trialsCueOff.speed_MMsecSpikeBef{i} = Spike.speed_MMsec(indSpikes);
        trialsCueOff.thPhaseHilbSpikeBef{i} = Spike.thPhaseHilb(indSpikes);
        trialsCueOff.thPhaseInterpSpikeBef{i} = Spike.thPhaseInterp(indSpikes); 
        
        trialsCueOff.stimOnLfpInd{i} = [];
        if(isfield(Laps,'stimOnLfpInd') & ~isempty(Laps.stimOnLfpInd{i}))
            trialsCueOff.stimOnLfpInd{i} = ...
                Laps.stimOnLfpInd{i}-trialsCueOff.startLfpInd(i)+1;
        end

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
                Spike.res < trialsCueOff.endLfpInd(i)); 
            trialsCueOff.res_LasttoCurTr{i} = Spike.res(indSpikesLastToCurTrEnd)-...
                trialsCueOff.startLfpInd(i)+1; 
            trialsCueOff.totclu_LasttoCurTr{i} = Spike.totclu(indSpikesLastToCurTrEnd);
        end
        
%         startCueTime = Laps.startLfpInd(indTr(i)) - ...
%             (indEndTrackLastTr+Laps.startLfpInd(indTr(i-1))-1);
%         disp(['TR ' num2str(i)]);
%         plot(speedAll);
%         hold on;
%         plot([indTrueStart indTrueStart],[0 600],'r');
%         plot([startCueTime startCueTime],[0 600], 'g');
%         hold off;
%         pause;
    end
    
    save([path fileNameCueOff],'trialsCueOff');
end

function [data,data1] = numOfConsecutiveOnes(arr)
    data = [];
    data1 = [];
    s = sprintf('%d', arr);
    %Reading the consequences of 1's from the string by using 0's as delimiters
    t1=textscan(s,'%s','delimiter','0','multipleDelimsAsOne',1);
    % Converting cell array of cell into a single cell array
    d = t1{:};
    % Computing the length of each run by going through the array and assigning
    % it into 
    for k = 1:length(d)
          data(k) = length(d{k});
    end
    
    t2=textscan(s,'%s','delimiter','1','multipleDelimsAsOne',1);
    % Converting cell array of cell into a single cell array
    f = t2{:};
    % Computing the length of each run by going through the array and assigning
    % it into 
    for k = 1:length(f)
          data1(k) = length(f{k});
    end
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