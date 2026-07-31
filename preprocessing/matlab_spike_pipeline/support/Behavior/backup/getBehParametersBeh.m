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
    
    %% exclude the first 5 and last 5 trials
    param.startTr = 5;
    param.endTr = length(trialsRun.xMM)-5;
    if(param.endTr < param.startTr)
        param.startTr = 2;
        param.endTr = length(trialsRun.xMM);
    end
    trialNo = 1:length(trialsRun.xMM);
    smoothSpan = 100;  
    minSpeed = 10;
    minTotSample = 10*sampleFq;
    behParRun.speedProfile = zeros(length(trialNo),minTotSample);
    behParRun.lickProfile = zeros(length(trialNo),minTotSample);
    for i = trialNo
        if(~isempty(trialsRun.xMM{i}))
%             disp(['Tr ' num2str(i)])
%             if(i == 66)
%                 a = 0;
%             end
            %% aligned to run onset
            behParRun.numSamples(i) = length(trialsRun.speed_MMsec{i});
            speedSM = trialsRun.speed_MMsec_SM{i};
            if(~isempty(trialsRun.pumpLfpInd{i}))
                indSpeed = speedSM(1:trialsRun.pumpLfpInd{i}(1) - ...
                    trialsRun.startLfpInd(i)+1) >= minSpeed;
                behParRun.pumpLfpInd(i) = trialsRun.pumpLfpInd{i}(1) - ...
                    trialsRun.startLfpInd(i)+1;
                behParRun.pumpMM(i) = trialsRun.xMM{i}(behParRun.pumpLfpInd(i));
            else
                indSpeed = speedSM >= minSpeed;
                behParRun.pumpLfpInd(i) = -1;
                behParRun.pumpMM(i) = -1;
            end
            
            [continuousRun,stopRun] = numOfConsecutiveOnes(indSpeed);
            behParRun.speedSM{i} = speedSM;
            behParRun.maxSpeed(i) = max(speedSM(indSpeed));
            behParRun.meanSpeed(i) = mean(speedSM(indSpeed));
            behParRun.stdSpeed(i) = std(speedSM(indSpeed));
            behParRun.maxRunLenT(i) = max(continuousRun)/sampleFq;
            behParRun.totRunLenT(i) = sum(continuousRun)/sampleFq;
            behParRun.numRun(i) = sum(continuousRun > 62.5); % number of runs that are longer than 50 ms
            
            accSM = trialsRun.accel_MMsecSq_SM{i}(indSpeed);
            behParRun.accSM{i} = accSM;
            behParRun.maxAcc(i) = max(accSM);
            behParRun.meanAcc(i) = mean(accSM);
            behParRun.stdAcc(i) = std(accSM);
            
            behParRun.totStopLenT(i) = sum(stopRun)/sampleFq;
            
            behParRun.startCueToRun(i) = (trialsRun.startLfpInd(i) - ...
                trialsCue.startLfpInd(i))/sampleFq; 
            behParRun.rewardToRun(i) = (trialsRun.startLfpInd(i) - ...
                trialsRew.startLfpInd(i))/sampleFq; 
            behParRun.rewarded(i) = trialsRew.goodTrial(i);
            if(trialsRun.goodTrial(i) < -1)
                behParRun.nonStop(i) = 1;
            else
                behParRun.nonStop(i) = 0;
            end
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
            
            %% aligned to reward 
            behParRew.numSamples(i) = length(trialsRew.speed_MMsec{i});
            speed = trialsRew.speed_MMsec{i};
            speedSM = smooth(speed,smoothSpan);
            indSpeed = speedSM >= minSpeed;
            [continuousRun,stopRun] = numOfConsecutiveOnes(indSpeed);
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
        behParRun.totStopLenT > 2 | behParRun.rewarded == -1 | behParRun.nonStop == 1;
    
    save([path fileName '_behPar.mat'],'behParRun','behParRew','behParCue','param');
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
