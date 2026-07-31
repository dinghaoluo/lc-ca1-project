function alignToCueOff2P(path, fileName, mazeSess)
% align the spikes in each trial to the running onset

    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    
    fileNameCueOff = [fileName '_alignCueOff_msess' num2str(mazeSess) '.mat'];  
           
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
        [speed,ac] = MazeSpeedAccel2P(xMM1,sampleFq);
        indTrueStart = befTime*sampleFq;
        
        trialsCueOff.startLfpInd(i) = indCueOff;

        trialsCueOff.endLfpInd(i) = Laps.endLfpInd(indTr(i)); 
        trialsCueOff.numSamples(i) = trialsCueOff.endLfpInd(i) - ...
            trialsCueOff.startLfpInd(i)+1;
         
        indEeg = trialsCueOff.startLfpInd(i):trialsCueOff.endLfpInd(i);
        trialsCueOff.xMM{i} = xMM1(indTrueStart+1:end)-xMM1(indTrueStart+1);
        trialsCueOff.speed_MMsec{i} = speed(indTrueStart+1:end);
        trialsCueOff.accel_MMsecSq{i} = ac(indTrueStart+1:end);
        trialsCueOff.startCueTime(i) = Laps.startLfpInd(indTr(i));
        
        trialsCueOff.dFF{i} = Track.dFFSM(indEeg,:);
        trialsCueOff.dFFGF{i} = Track.dFFGF(indEeg,:); % added on 9/6/2023, gaussian filtered dFF
        trialsCueOff.F{i} = Track.F(indEeg,:);
        trialsCueOff.Fneu{i} = Track.Fneu(indEeg,:);
        trialsCueOff.spikes{i} = Track.spikesSM(indEeg,:);
        
        indEegBef = trialsCueOff.startLfpInd(i)-befTime*sampleFq:trialsCueOff.startLfpInd(i)-1;
        trialsCueOff.speed_MMsecBef{i} = Track.speed_MMsec(indEegBef);
        trialsCueOff.dFFBef{i} = Track.dFFSM(indEeg,:);
        trialsCueOff.dFFBefGF{i} = Track.dFFGF(indEeg,:); % added on 9/6/2023, gaussian filtered dFF
        trialsCueOff.FBef{i} = Track.F(indEeg,:);
        trialsCueOff.FneuBef{i} = Track.Fneu(indEeg,:);
        trialsCueOff.spikesBef{i} = Track.spikesSM(indEeg,:);
        
        licks = Laps.lickLfpInd{indTr(i)};
        indLick = licks > trialsCueOff.startLfpInd(i);
        trialsCueOff.lickLfpInd{i} = licks(indLick);
        
        trialsCueOff.stimOnLfpInd{i} = [];
        if(isfield(Laps,'stimOnLfpInd') & ~isempty(Laps.stimOnLfpInd{i}))
            trialsCueOff.stimOnLfpInd{i} = ...
                Laps.stimOnLfpInd{i}-trialsCueOff.startLfpInd(i)+1;
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
    
    save([path fileNameCueOff],'trialsCueOff','-v7.3');
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