function getBehParametersBeh(path,fileName,mazeSess)
% get the behavior parameters for each trial
% e.g. getBehParameters('./','A011-20190218-01_DataStructure_mazeSection1_TrialType1',1)
    
    fullPath = [path fileName '_alignRun_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to run file does not exist');
        return;
    end
    load(fullPath,'trialsRun');
    
    fullPath = [path fileName '_alignRew_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to reward file does not exist');
        return;
    end
    load(fullPath,'trialsRew');
    
    fullPath = [path fileName '_alignCue_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to cue file does not exist');
        return;
    end
    load(fullPath,'trialsCue');
    
    GlobalConst;
    
    %% exclude the first few trials, changed on 7/2/2020  
    param.endTr = length(trialsRun.xMM);
    if(unique(trialsRun.mazeTypeMod)==0) %control trials
        param.startTr = param.endTr-40;
        if(param.startTr < 10)
            param.startTr = 10;
        end
    else
        param.startTr = 2; % exclude the first trial which is empty
    end
    trialNo = 1:length(trialsRun.xMM); % avoid the first trial
    smoothSpan = 100;  
    minSpeed = 10;
    minTotSample = 10*sampleFq;
    behParRun.speedProfile = zeros(length(trialNo),minTotSample);
    behParRun.lickProfile = zeros(length(trialNo),minTotSample);
    behParRun.mazeType = trialsRun.mazeType;
    behParRun.mazeTypeMod = trialsRun.mazeTypeMod;
    behParRun.stimOnLfpInd = zeros(1,length(trialNo));
    behParRun.stimOffLfpInd = zeros(1,length(trialNo));
    behParRun.stimAcrossTrialRun = zeros(1,length(trialNo));
    behParRun.stimAcrossTrialCue = zeros(1,length(trialNo));
    behParRun.stimAcrossTrialRew = zeros(1,length(trialNo));
    for i = trialNo
        disp(['Trial ' num2str(i)]);
        if(~isempty(trialsRun.xMM{i}))
%             disp(['Tr ' num2str(i)])
            if(i == 52)
                a = 0;
            end
            %% aligned to run onset
            behParRun.numSamples(i) = length(trialsRun.speed_MMsec{i});
            speedSM = trialsRun.speed_MMsec_SM{i};
            behParRun.lowSpeed{i} = speedSM < minSpeed;
            if(~isempty(trialsRun.pumpLfpInd{i}))
                indSpeed = speedSM(1:trialsRun.pumpLfpInd{i}(1) - ...
                    trialsRun.startLfpInd(i)+1) >= minSpeed;    
                behParRun.pumpLfpInd(i) = trialsRun.pumpLfpInd{i}(1) - ...
                    trialsRun.startLfpInd(i)+1;
                %% changed by Yingxue on 2/22/2022
                if(behParRun.pumpLfpInd(i) > 0) 
                    behParRun.pumpMM(i) = trialsRun.xMM{i}(behParRun.pumpLfpInd(i));
                else
                    behParRun.pumpMM(i) = -1;
                end
            else
                indSpeed = speedSM >= minSpeed;
                behParRun.pumpLfpInd(i) = -1;
                behParRun.pumpMM(i) = -1;
            end
            
            %% changed by Yingxue on 2/22/2022
            behParRun.speedSM{i} = speedSM;
            behParRun.maxSpeed(i) = max(speedSM);
            behParRun.meanSpeed(i) = mean(speedSM);
            behParRun.stdSpeed(i) = std(speedSM);
            accSM = trialsRun.accel_MMsecSq_SM{i}(indSpeed);
            behParRun.accSM{i} = accSM;
            if(isempty(indSpeed))
                continuousRun = 0;
                stopRun = 0;
                
                behParRun.maxSpeedRun(i) = 0;
                behParRun.meanSpeedRun(i) = 0;
                behParRun.stdSpeedRun(i) = 0;
                
                behParRun.maxRunLenT(i) = 0;
                behParRun.totRunLenT(i) = 0;
                behParRun.numRun(i) = 0; % number of runs that are longer than 50 ms
                
                behParRun.maxAcc(i) = 0;
                behParRun.meanAcc(i) = 0;
                behParRun.stdAcc(i) = 0;
            else
                [continuousRun,stopRun] = numOfConsecutiveOnes(indSpeed);
                
                behParRun.maxSpeedRun(i) = max(speedSM(indSpeed));
                behParRun.meanSpeedRun(i) = mean(speedSM(indSpeed));
                behParRun.stdSpeedRun(i) = std(speedSM(indSpeed));
                
                behParRun.maxRunLenT(i) = max(continuousRun)/sampleFq;
                behParRun.totRunLenT(i) = sum(continuousRun)/sampleFq;
                behParRun.numRun(i) = sum(continuousRun > 62.5); % number of runs that are longer than 50 ms
                
                behParRun.maxAcc(i) = max(accSM);
                behParRun.meanAcc(i) = mean(accSM);
                behParRun.stdAcc(i) = std(accSM);
            end
            %%
            
            behParRun.totStopLenT(i) = sum(stopRun)/sampleFq;
            
            behParRun.startCueToRun(i) = (trialsRun.startLfpInd(i) - ...
                trialsCue.startLfpInd(i))/sampleFq; 
            behParRun.rewardToRun(i) = (trialsRun.startLfpInd(i) - ...
                trialsRew.startLfpInd(i))/sampleFq; 
            
            if(behParRun.pumpMM(i) == -1) % changed by Yingxue 02/22/2022
                behParRun.rewarded(i) = -1;
            else
                behParRun.rewarded(i) = 1;
            end
            %% changed by Yingxue on 1/23/2022
            if(trialsRun.goodTrial(i) == -2)
                behParRun.nonStop(i) = 1;
            else
                behParRun.nonStop(i) = 0;
            end
            if(trialsRun.goodTrial(i) == -3)
                behParRun.nonFullStop(i) = 1;
            else
                behParRun.nonFullStop(i) = 0;
            end
            %%
            if(length(trialsRun.lickLfpInd{i}) >= 5)
                firstFiveLicks = trialsRun.lickLfpInd{i}(1:5);
            else
                firstFiveLicks = trialsRun.lickLfpInd{i};
            end
            behParRun.med1stFiveLickDist(i) = ...
                median(trialsRun.xMM{i}(firstFiveLicks));
            behParRun.medLickDist(i) = ...
                median(trialsRun.xMM{i}(trialsRun.lickLfpInd{i}));
            if(~isempty(trialsRun.pumpLfpInd{i}))
                indLick = trialsRun.lickLfpInd{i} <= behParRun.pumpLfpInd(i);
            else
                indLick = trialsRun.lickLfpInd{i} <= behParRun.numSamples(i);
            end
            behParRun.numLicksBefRew(i) = sum(indLick);
            behParRun.numLicksRew(i) = length(trialsRun.lickLfpInd{i}) - ...
                sum(indLick);
            behParRun.lickLfpIndBefRew{i} = trialsRun.lickLfpInd{i}(indLick);
            behParRun.lickLfpIndRew{i} = trialsRun.lickLfpInd{i}(~indLick);
            behParRun.lickMMBefRew{i} = trialsRun.xMM{i}(trialsRun.lickLfpInd{i}(indLick));
            behParRun.lickMMRew{i} = trialsRun.xMM{i}(trialsRun.lickLfpInd{i}(~indLick));
            
            ind = find(indLick == 1);
            if(sum(indLick) >= 5)
                indTmp = ind(1:5);
            else
                indTmp = ind;
            end
            behParRun.med1stFiveLickDistBefRew(i) = ...
                median(trialsRun.xMM{i}(trialsRun.lickLfpInd{i}(indTmp)));
            behParRun.medLickDistBefRew(i) = ...
                median(trialsRun.xMM{i}(trialsRun.lickLfpInd{i}(indLick)));
            if(i >= param.startTr & i <= param.endTr)
                if(behParRun.numSamples(i) >= minTotSample)
                    behParRun.speedProfile(i,:) = ...
                        behParRun.speedSM{i}(1:minTotSample);
                else
                    behParRun.speedProfile(i,1:behParRun.numSamples(i)) = ...
                        behParRun.speedSM{i}';
                end
                indTmp = find(trialsRun.lickLfpInd{i} <= minTotSample);
                behParRun.lickProfile(i,trialsRun.lickLfpInd{i}(indTmp)) = 1;
            end
            
            %% added by Yingxue on 1/25/2022
            if(behParRun.mazeType(i) ~= 0)
                behParRun.stimOnLfpInd(i) = trialsRun.stimOnLfpInd{i}(1);
                behParRun.stimOffLfpInd(i) = trialsRun.stimOffLfpIndExact(i);
                if(i+1 <= trialNo(end))
                    if(trialsRun.stimOffLfpIndExact(i) > trialsRun.startLfpInd(i+1))
                        behParRun.stimAcrossTrialRun(i) = 1;
                    end
                    if(trialsRun.stimOffLfpIndExact(i) > trialsCue.startLfpInd(i+1))
                        behParRun.stimAcrossTrialCue(i) = 1;
                    end
                    if(trialsRun.stimOffLfpIndExact(i) > trialsRew.startLfpInd(i+1))
                        behParRun.stimAcrossTrialRew(i) = 1;
                    end
                end
            end
            
            %% aligned to reward 
            behParRew.numSamples(i) = length(trialsRew.speed_MMsec{i});
            speed = trialsRew.speed_MMsec{i};
            speedSM = smooth(speed,smoothSpan);
            indSpeed = speedSM >= minSpeed;
            [continuousRun,stopRun] = numOfConsecutiveOnes(indSpeed);
            behParRew.maxSpeedRun(i) = max(speedSM(indSpeed));
            behParRew.meanSpeedRun(i) = mean(speedSM(indSpeed));
            behParRew.stdSpeedRun(i) = std(speedSM(indSpeed));
            behParRew.maxSpeed(i) = max(speedSM);
            behParRew.meanSpeed(i) = mean(speedSM);
            behParRew.stdSpeed(i) = std(speedSM);
            behParRew.maxRunLenT(i) = max(continuousRun)/sampleFq;
            behParRew.totRunLenT(i) = sum(continuousRun)/sampleFq;
            behParRew.numRun(i) = sum(continuousRun > 62.5); % number of runs that are longer than 50 ms
            
            acc = trialsRew.accel_MMsecSq{i};
            accSM = smooth(acc,smoothSpan);
            behParRew.maxAcc(i) = max(accSM);
            behParRew.meanAcc(i) = mean(accSM);
            behParRew.stdAcc(i) = std(accSM);
            
            behParRew.totStopLenTRew(i) = sum(stopRun)/sampleFq;
            
            %% aligned to cue onset
            behParCue.numSamples(i) = length(trialsCue.speed_MMsec{i});
            speed = trialsCue.speed_MMsec{i};
            speedSM = smooth(speed,smoothSpan);
            indSpeed = speedSM >= minSpeed;
            [continuousRun,stopRun] = numOfConsecutiveOnes(indSpeed);
            behParCue.maxSpeedRun(i) = max(speedSM(indSpeed));
            behParCue.meanSpeedRun(i) = mean(speedSM(indSpeed));
            behParCue.stdSpeedRun(i) = std(speedSM(indSpeed));
            behParCue.maxSpeed(i) = max(speedSM);
            behParCue.meanSpeed(i) = mean(speedSM);
            behParCue.stdSpeed(i) = std(speedSM);
            behParCue.maxRunLenT(i) = max(continuousRun)/sampleFq;
            behParCue.totRunLenT(i) = sum(continuousRun)/sampleFq;
            behParCue.numRun(i) = sum(continuousRun > 62.5); % number of runs that are longer than 50 ms
            
            acc = trialsRew.accel_MMsecSq{i};
            accSM = smooth(acc,smoothSpan);
            behParCue.maxAcc(i) = max(accSM);
            behParCue.meanAcc(i) = mean(accSM);
            behParCue.stdAcc(i) = std(accSM);
            
            behParCue.totStopLenT(i) = sum(stopRun)/sampleFq;
            
            if(length(trialsCue.lickLfpInd{i}) >= 5)
                firstFiveLicks = trialsCue.lickLfpInd{i}(1:5);
            else
                firstFiveLicks = trialsCue.lickLfpInd{i};
            end
            behParCue.med1stFiveLickDist(i) = ...
                median(trialsCue.xMM{i}(firstFiveLicks));
            
            behParCue.medLickDist(i) = ...
                median(trialsCue.xMM{i}(trialsCue.lickLfpInd{i}));
            
            if(~isempty(trialsCue.pumpLfpInd{i}))
                indLick = trialsCue.lickLfpInd{i} <= trialsCue.pumpLfpInd{i}(1)...
                    - trialsCue.startLfpInd(i)+1;
            else
                indLick = trialsCue.lickLfpInd{i} <= behParCue.numSamples(i);
            end
            behParCue.numLicksBefRew(i) = sum(indLick);
            behParCue.numLicksRew(i) = length(trialsCue.lickLfpInd{i}) - ...
                sum(indLick);
            behParCue.lickLfpIndBefRew{i} = trialsCue.lickLfpInd{i}(indLick);
            behParCue.lickLfpIndRew{i} = trialsCue.lickLfpInd{i}(~indLick);
            behParCue.lickMMBefRew{i} = trialsCue.xMM{i}(trialsCue.lickLfpInd{i}(indLick));
            behParCue.lickMMRew{i} = trialsCue.xMM{i}(trialsCue.lickLfpInd{i}(~indLick));
            
            behParCue.medLickDistBefRew(i) = ...
                median(trialsCue.xMM{i}(trialsCue.lickLfpInd{i}(indLick)));
            ind = find(indLick == 1);
            if(sum(indLick) >= 5)
                indTmp = ind(1:5);
            else
                indTmp = ind;
            end
            behParCue.med1stFiveLickDistBefRew(i) = ...
                median(trialsCue.xMM{i}(trialsCue.lickLfpInd{i}(indTmp)));
            
        end
    end
%     figure;
%     timeBin = -3:0.5:50;
%     hist(behParRun.startCueToRun,timeBin);
%     title(fileName)
    behParRun.indTrBadBeh = behParRun.totRunLenT > 13 | behParRun.numRun > 10 | ...
        behParRun.totStopLenT > 2 | behParRun.rewarded == -1 | behParRun.nonStop == 1 ...
        | behParRun.nonFullStop == 1;
    
    save([path fileName '_behPar_msess' num2str(mazeSess) '.mat'],...
        'behParRun','behParRew','behParCue','param');
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
