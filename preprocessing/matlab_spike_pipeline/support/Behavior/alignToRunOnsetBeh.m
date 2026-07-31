function alignToRunOnset(path, fileName, mazeSess)
% align each trial to the running onset

    fileNameRun = [fileName '_alignRun_msess' num2str(mazeSess) '.mat'];  
  
    fileNameBehElec = [fileName '_BehavElectrDataLFP.mat'];           
    fullPath = [path fileNameBehElec];
    if(exist(fullPath) == 0)
        disp('The recording BehElectDataLFP file does not exist.');
        disp(fullPath)
        return;
    end
    load(fullPath);
    
    GlobalConst;
    
    minSpeed1 = 100;
    minSpeed = 10;
    smoothSpan = 100;
    nSampBef = 3*sampleFq;
    trialsRun = [];
    indTr = find(Laps.mazeSess == mazeSess);
    disp(['number of trials = ' num2str(length(indTr))]);
    for i = 1:length(indTr)
        trialsRun.pumpLfpInd{i} = Laps.pumpLfpInd{indTr(i)};
        if(i == 1)
            indStart = -1;
            trialsRun.goodTrial(i) = indStart;
            continue;
        end
        xMM = Track.xMMAll(Laps.startLfpInd(indTr(i)):Laps.endLfpInd(indTr(i)));    
        xMM1 = xMMContinuous(xMM);
        speed = MazeSpeedAccel(xMM1,sampleFq);
        speedSM = smooth(speed,smoothSpan);
        
        % changed by Yingxue on 7/9/2020
        trackLastTr = Track.xMMAll(Laps.startLfpInd(indTr(i-1)):Laps.endLfpInd(indTr(i-1)));
        indStartLastTr = find(diff(trackLastTr) < -Laps.trackLen(indTr(i-1))/2);
         % check the reset of distance from last trial
        if(length(indStartLastTr) > 1)
            indEndTrackLastTr = find(trackLastTr(indStartLastTr(end-1)+1:indStartLastTr(end))...
                    > Laps.trackLen(indTr(i-1)),1)+indStartLastTr(end-1);
        elseif(length(indStartLastTr) == 1)
            indEndTrackLastTr = find(trackLastTr(indStartLastTr(end)+1:end)...
                    > Laps.trackLen(indTr(i-1)),1)+indStartLastTr(end);
            if(isempty(indEndTrackLastTr)) 
            % if the reset only occurs at the end of the last trial, there 
            % is no reset at the beginning of the last trial
                indEndTrackLastTr = find(trackLastTr...
                    > Laps.trackLen(indTr(i-1)),1);
            end
        else
            indEndTrackLastTr = find(trackLastTr...
                    > Laps.trackLen(indTr(i-1)),1);
        end
        
        xMMAll = Track.xMMAll(indEndTrackLastTr+Laps.startLfpInd(indTr(i-1))-1:Laps.endLfpInd(indTr(i)));
        xMM1All = xMMContinuous(xMMAll);
        [speedAll,acAll] = MazeSpeedAccel(xMM1All,sampleFq);
        speedAllSM = smooth(speedAll,smoothSpan);
        accAllSM = smooth(acAll,smoothSpan);
        lenDiff = length(xMMAll) - length(xMM);
        
        indSpeed = speedSM > minSpeed1; % 100 mm/s min speed
        [continuousRun,stopRun] = numOfConsecutiveOnes(indSpeed);
        indFirstRun = find(indSpeed == 1,1);
        indContiRun = find(continuousRun > sampleFq*0.3,1);
  
        if(indFirstRun > 1)
            if(indContiRun > 1)
                indStart = sum(continuousRun(1:indContiRun-1)) + ...
                    sum(stopRun(1:indContiRun));
                indStart = indStart(end)+1+lenDiff;
            else
                indStart = indFirstRun+lenDiff;
            end
        else
            if(indContiRun > 1)
                indStart = sum(continuousRun(1:indContiRun-1)) + ...
                    sum(stopRun(1:indContiRun-1));
                indStart = indStart(end)+1+lenDiff;
            else
                indSpeed = speedAllSM(1:lenDiff) > minSpeed1;
                indLastRunStart = find(indSpeed == 0,1,'last');
                if(isempty(indLastRunStart))
                    indStart = -2; % continuously run during the reward and blackout
                else
                    indStart = indLastRunStart;
                end  
            end
        end
        
        if(indStart > 0)    
            indTrueStart = find(speedAllSM(1:indStart-1) <= minSpeed,1,'last');
            if(isempty(indTrueStart))
                indStart = -3; % never fully stopped during the reward and blackout 
                indTrueStart = 0;
            else
                indStart = indTrueStart;
            end
        else
            indTrueStart = 0;
        end
        
        trialsRun.goodTrial(i) = indStart;
        if(indStart > 0)
            trialsRun.startLfpInd(i) = indTrueStart+1+Laps.startLfpInd(indTr(i-1))+ ...
                indEndTrackLastTr-1; %% shifted by +1 on 11/27/2019
        else
            trialsRun.startLfpInd(i) = Laps.startLfpInd(indTr(i-1))+ ...
                indEndTrackLastTr-1;
        end
        trialsRun.endLfpInd(i) = Laps.endLfpInd(indTr(i)); 
        trialsRun.numSamples(i) = trialsRun.endLfpInd(i) - ...
            trialsRun.startLfpInd(i)+1;
        trialsRun.xMM{i} = xMM1All(indTrueStart+1:end)-xMM1All(indTrueStart+1);
        trialsRun.speed_MMsec{i} = speedAll(indTrueStart+1:end);
        trialsRun.speed_MMsec_SM{i} = speedAllSM(indTrueStart+1:end);
        trialsRun.accel_MMsecSq{i} = acAll(indTrueStart+1:end);
        trialsRun.accel_MMsecSq_SM{i} = accAllSM(indTrueStart+1:end);
        trialsRun.startCueTime(i) = Laps.startLfpInd(indTr(i));
        
        licks = [Laps.lickLfpInd{indTr(i-1)}; Laps.lickLfpInd{indTr(i)}];
        indLick = licks > trialsRun.startLfpInd(i);
        trialsRun.lickLfpInd{i} = licks(indLick)-trialsRun.startLfpInd(i)+1;
        
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
    
    save([path fileNameRun],'trialsRun');
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