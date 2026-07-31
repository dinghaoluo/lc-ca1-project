function alignToRunOnset2P(path, fileName, mazeSess)
% align the spikes in each trial to the running onset

    indexFileName = findstr(fileName, '.mat');
    if(~isempty(indexFileName))
        fileName = fileName(1:indexFileName(end)-1);
    end
    
    fileNameRun = [fileName '_alignRun_msess' num2str(mazeSess) '.mat'];  
        
    indRecName = strfind(fileName, '_');
    fileNameBehElec = [fileName(1:indRecName) 'Behav2PDataLFP.mat'];
            
    fullPath = [path fileNameBehElec];
    if(exist(fullPath) == 0)
        disp('The recording Behav2PDataLFP file does not exist.');
        disp(fullPath)
        return;
    end
    load(fullPath,'Laps','Track');
    
    GlobalConst2P;
    
    minSpeed = 10;
    smoothSpan = 100;
    trialsRun = [];
    indTr = find(Laps.mazeSess == mazeSess);
    for i = 1:length(indTr)
        disp(num2str(indTr(i)))
        if(i == 34)
            a = 1;
        end
        trialsRun.pumpLfpInd{i} = Laps.pumpLfpInd{indTr(i)};
        if(i == 1)
            indStart = -1;
            trialsRun.goodTrial(i) = indStart;
            continue;
        end
        xMM = Track.xMMAll(Laps.startLfpInd(indTr(i)):Laps.endLfpInd(indTr(i)));    
        xMM1 = xMMContinuous(xMM);
        speed = MazeSpeedAccel2P(xMM1,sampleFq);
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
        [speedAll,acAll] = MazeSpeedAccel2P(xMM1All,sampleFq);
        speedAllSM = smooth(speedAll,smoothSpan);
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
                %% changed by Yingxue on 1/23/2022, calculate the index with lowest speed
                % instead of using reward location as run start
                [~,indTmp] = min(speedAllSM(1:indStart-1));
                indTrueStart = indTmp(end); 
                indStart = -3; % no full stop
            else
                indStart = indTrueStart;
            end
        else % no stop
            %% changed by Yingxue on 1/23/2022, calculate the index with lowest speed
            % instead of using reward location as run start
            [~,indTmp] = min(speedAllSM(1:lenDiff));
            indTrueStart = indTmp(end);
        end
       
        %% changed by Yingxue on 1/23/2022, calculate the index with lowest speed
        % instead of using reward location as run start
        trialsRun.goodTrial(i) = indStart;
        trialsRun.startLfpInd(i) = indTrueStart+Laps.startLfpInd(indTr(i-1))+ ...
                indEndTrackLastTr-1; %% shifted by +1 on 11/27/2019
        
        figure(1);
        plot(speedAllSM);
        hold on;
        plot(indTrueStart,speedAllSM(indTrueStart),'ro');
        hold off;
        %%
        
        trialsRun.endLfpInd(i) = Laps.endLfpInd(indTr(i)); 
        trialsRun.numSamples(i) = trialsRun.endLfpInd(i) - ...
            trialsRun.startLfpInd(i)+1;
        trialsRun.dFF{i} = Track.dFFSM(trialsRun.startLfpInd(i):trialsRun.endLfpInd(i),:);
        trialsRun.dFFGF{i} = Track.dFFGF(trialsRun.startLfpInd(i):trialsRun.endLfpInd(i),:); % added on 9/6/2023, gaussian filtered dFF
        trialsRun.F{i} = Track.F(trialsRun.startLfpInd(i):trialsRun.endLfpInd(i),:);
        trialsRun.Fneu{i} = Track.Fneu(trialsRun.startLfpInd(i):trialsRun.endLfpInd(i),:);
        trialsRun.spikes{i} = Track.spikesSM(trialsRun.startLfpInd(i):trialsRun.endLfpInd(i),:);
        trialsRun.xMM{i} = xMM1All(indTrueStart+1:end)-xMM1All(indTrueStart+1);
        trialsRun.speed_MMsec{i} = speedAll(indTrueStart+1:end);
        trialsRun.accel_MMsecSq{i} = acAll(indTrueStart+1:end);
        trialsRun.startCueTime(i) = Laps.startLfpInd(indTr(i));
        
        indEegBef = trialsRun.startLfpInd(i)-nSampBef:trialsRun.startLfpInd(i)-1;
        trialsRun.dFFBef{i} = Track.dFFSM(indEegBef,:);
        trialsRun.dFFBefGF{i} = Track.dFFGF(indEegBef,:); % added on 9/6/2023, gaussian filtered dFF
        trialsRun.FBef{i} = Track.F(indEegBef,:);
        trialsRun.FneuBef{i} = Track.Fneu(indEegBef,:);
        trialsRun.spikesBef{i} = Track.spikesSM(indEegBef,:);
        trialsRun.speed_MMsecBef{i} = Track.speed_MMsec(indEegBef);
        
        licks = [Laps.lickLfpInd{indTr(i-1)}; Laps.lickLfpInd{indTr(i)}];
        indLick = licks > trialsRun.startLfpInd(i);
        trialsRun.lickLfpInd{i} = licks(indLick);
                
        trialsRun.stimOnLfpInd{i} = [];
        if(isfield(Laps,'stimOnLfpInd') & ~isempty(Laps.stimOnLfpInd{i}))
            trialsRun.stimOnLfpInd{i} = ...
                Laps.stimOnLfpInd{i}-trialsRun.startLfpInd(i)+1;
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
    
    save([path fileNameRun],'trialsRun','-v7.3');
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