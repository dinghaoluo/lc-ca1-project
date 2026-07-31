function getBehParameters2P(path,fileName,mazeSess)
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
    
    fullPath = [path fileName '_alignCueOff_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to cue file does not exist');
        return;
    end
    load(fullPath,'trialsCueOff');
    
    GlobalConst2P;
    
    trialNo = 1:length(trialsRun.xMM);
    smoothSpan = 100;  
    minSpeed = 10;
    for i = trialNo
        if(~isempty(trialsRun.xMM{i}))
            behPar.numSamplesRun(i) = length(trialsRun.speed_MMsec{i});
            behPar.numSamplesRew(i) = length(trialsRew.speed_MMsec{i});
            behPar.numSamplesCue(i) = length(trialsCue.speed_MMsec{i});
            if(trialsCueOff.startLfpInd(i) ~= -1) % in case there is no cue off signal
                behPar.numSamplesCueOff(i) = length(trialsCueOff.speed_MMsec{i});
            else
                behPar.numSamplesCueOff(i) = 0;
            end
            speed = trialsRun.speed_MMsec{i};
            speedSM = smooth(speed,smoothSpan);
            indSpeed = speedSM >= minSpeed;
            [continuousRun,stopRun] = numOfConsecutiveOnes(indSpeed);
            behPar.maxSpeed(i) = max(speedSM);
            behPar.meanSpeed(i) = mean(speedSM);
            behPar.maxRunLenT(i) = max(continuousRun)/sampleFq;
            behPar.totRunLenT(i) = sum(continuousRun)/sampleFq;
            behPar.numRun(i) = sum(continuousRun > 62.5); % number of runs that are longer than 50 ms
            
            acc = trialsRun.accel_MMsecSq{i};
            accSM = smooth(acc,smoothSpan);
            behPar.maxAcc(i) = max(accSM);
            behPar.meanAcc(i) = mean(accSM);
            
            behPar.totStopLenT(i) = sum(stopRun)/sampleFq;
            
            speed = trialsRew.speed_MMsec{i};
            speedSM = smooth(speed,smoothSpan);
            indSpeed = speedSM >= minSpeed;
            [continuousRun,stopRun] = numOfConsecutiveOnes(indSpeed);
                        
            behPar.startCueToRun(i) = (trialsRun.startLfpInd(i) - ...
                trialsCue.startLfpInd(i))/sampleFq; 
            behPar.rewardToRun(i) = (trialsRun.startLfpInd(i) - ...
                trialsRew.startLfpInd(i))/sampleFq; 
            %% changed by Yingxue on 02/22/2022
            if(~isempty(trialsRun.pumpLfpInd{i}))
                behPar.rewarded(i) = 1;
            else
                behPar.rewarded(i) = -1;
            end
            %% changed by Yingxue on 1/23/2022
            if(trialsRun.goodTrial(i) == -2)
                behPar.nonStop(i) = 1;
            else
                behPar.nonStop(i) = 0;
            end
            if(trialsRun.goodTrial(i) == -3)
                behPar.nonFullStop(i) = 1;
            else
                behPar.nonFullStop(i) = 0;
            end
            %%
            if(length(trialsRun.lickLfpInd{i}) >= 5)
                firstFiveLicks = trialsRun.lickLfpInd{i}(1:5);
            else
                firstFiveLicks = trialsRun.lickLfpInd{i};
            end
            behPar.med1stFiveLickDist(i) = ...
                median(trialsRun.xMM{i}(firstFiveLicks-trialsRun.startLfpInd(i)));
            behPar.medLickDist(i) = ...
                median(trialsRun.xMM{i}(trialsRun.lickLfpInd{i}-trialsRun.startLfpInd(i)));
            
            if(isempty(trialsRun.stimOnLfpInd{i}))
                behPar.stimOn(i) = 0;
            else
                behPar.stimOn(i) = 1;
            end
        end
    end
    figure;
    timeBin = -3:0.5:50;
    hist(behPar.startCueToRun,timeBin);
    title(fileName)
    behPar.indTrBadBeh = behPar.totRunLenT > 13 | behPar.numRun > 10 | ...
        behPar.totStopLenT > 2 | behPar.rewarded == -1 | behPar.nonStop == 1;
    
    save([path fileName '_behPar_msess' num2str(mazeSess) '.mat'],'behPar');
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
