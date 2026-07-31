function RunSpeedOverTimeBeh_opto(path, fileName, onlyRun, mazeSess)

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
    
    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat'];
    if(exist(fullPath) == 0)
        disp('The behPar file does not exist');
        return;
    end
    load(fullPath,'param','behParCue');
    
    GlobalConst;
    paramC.trialLenT = 20; %sec
    paramC.timeIncBef = ((0:trialLenT*sampleFq-1) - nSampBef)/sampleFq;
    paramC.time = (1:(trialLenT*sampleFq - nSampBef))/sampleFq;
    paramC.timeRew = (1:(trialLenT*sampleFq))/sampleFq;
        
    nTrials = length(trialsRun.goodTrial);
    
    trInd = param.startTr:param.endTr;
    %% removing the grooming trial. Added by Yingxue on 4/22/2021
    indGroomTr = max(behParCue.numSamples) > 30*sampleFq; % a trial longer than 15 sec
    if(sum(indGroomTr) > 0)
        [~,indGroomTr] = max(behParCue.numSamples);
        trInd = setdiff(trInd,indGroomTr);
    end
    %%
    nTr = length(trInd);
    
    mazeT = unique(trialsRun.mazeType);
    mazeT = mazeT(mazeT ~= 0);

    %% new code
%     %%% changed by Yingxue on 12/15/2020
%     trIndOpt = [];
%     if(~isempty(mazeT))
%         if(mazeT == 4)
%             trIndOpt = find(trialsRun.mazeType ~= 0);
%             trIndOpt = trIndOpt+1;  
%             if(trIndOpt(end) > length(trialsRun.mazeType))
%                 trIndOpt = trIndOpt(1:end-1);
%             end
%             trIndOptCtrl = setdiff(trIndOpt(1):trIndOpt(end),trIndOpt);        
%         else
%             trIndOpt = find(trialsRun.mazeType ~= 0);
%             trIndOptCtrl = setdiff(trIndOpt(1):trIndOpt(end),trIndOpt); 
%         end
%     end
    
    %% old code
    trIndOpt = [];
    if(~isempty(mazeT))
        if(mazeT == 4)
            trIndOpt = trInd(trialsRun.mazeType(trInd) ~= 0)+1;   
            if(trIndOpt(end) > trInd(end))
                trIndOpt = trIndOpt(1:end-1);
            end
            trIndOptCtrl = setdiff(trInd,trIndOpt);        
        else
            trIndOpt = trInd(trialsRun.mazeType(trInd) ~= 0);
            trIndOptCtrl = trInd(trialsRun.mazeType(trInd) == 0);
        end
        %% added by Yingxue on1/25/2022
        trIndOptCtrl1 = trIndOpt+1;
        if(trIndOptCtrl1(end) > trInd(end))
            trIndOptCtrl1 = trIndOptCtrl1(1:end-1);
        end
    end    
    
    speedOverTime.runSpeedAlignedRun = ...
        speedAlignedRun(trialsRun,1:nTrials,paramC.trialLenT*sampleFq);
    
    speedOverTime.trInd = trInd;    
    speedOverTime.meanRun = mean(speedOverTime.runSpeedAlignedRun(trInd,:));
    speedOverTime.stdRun = std(speedOverTime.runSpeedAlignedRun(trInd,:));
    speedOverTime.SEMRun = std(speedOverTime.runSpeedAlignedRun(trInd,:))/sqrt(nTr);
    
    timeStep = paramC.timeIncBef;
    speedOverTime.indBaseline = find(timeStep >= -3 & timeStep < -2);
    speedOverTime.meanRunSpeedBL = mean(speedOverTime.runSpeedAlignedRun(trInd,speedOverTime.indBaseline),2);
    speedOverTime.indBefRun = find(timeStep >= -0.5 & timeStep < 0);
    speedOverTime.meanRunSpeedBefRun = mean(speedOverTime.runSpeedAlignedRun(trInd,speedOverTime.indBefRun),2);
    speedOverTime.ind0to1 = find(timeStep >= 0 & timeStep < 1);
    speedOverTime.meanRunSpeed0to1 = mean(speedOverTime.runSpeedAlignedRun(trInd,speedOverTime.ind0to1),2);
    speedOverTime.ind3to5 = find(timeStep >= 3 & timeStep < 5);
    speedOverTime.meanRunSpeed3to5 = mean(speedOverTime.runSpeedAlignedRun(trInd,speedOverTime.ind3to5),2);
    
    speedOverTime.runSpeedAlignedRew = ...
        speedAlignedRew(trialsRew,1:nTrials,paramC.trialLenT*sampleFq);
       
    speedOverTime.meanRew = mean(speedOverTime.runSpeedAlignedRew(trInd,:));
    speedOverTime.stdRew = std(speedOverTime.runSpeedAlignedRew(trInd,:));
    speedOverTime.SEMRew = std(speedOverTime.runSpeedAlignedRew(trInd,:))/sqrt(nTr);
    
    timeStep = paramC.timeRew;
    speedOverTime.ind0to1Rew = find(timeStep > 0 & timeStep <= 1);
    speedOverTime.meanRewSpeed0to1 = mean(speedOverTime.runSpeedAlignedRew(trInd,speedOverTime.ind0to1Rew),2);
    speedOverTime.ind1to2Rew = find(timeStep > 1 & timeStep <= 2);
    speedOverTime.meanRewSpeed1to2 = mean(speedOverTime.runSpeedAlignedRew(trInd,speedOverTime.ind1to2Rew),2);
    speedOverTime.ind2to3Rew = find(timeStep > 2 & timeStep <= 3);
    speedOverTime.meanRewSpeed2to3 = mean(speedOverTime.runSpeedAlignedRew(trInd,speedOverTime.ind2to3Rew),2);
    speedOverTime.ind3to5Rew = find(timeStep > 3 & timeStep <= 5);
    speedOverTime.meanRewSpeed3to5 = mean(speedOverTime.runSpeedAlignedRew(trInd,speedOverTime.ind3to5Rew),2);
    
    if(~isempty(trIndOpt))
        nTrIndOpt = length(trIndOpt); %% changed by Yingxue on 1/25/2022 (SEM was caculated using nTr before)
        speedOverTimeOpt.trInd = trIndOpt;
        speedOverTimeOpt.meanRun = mean(speedOverTime.runSpeedAlignedRun(trIndOpt,:));
        speedOverTimeOpt.stdRun = std(speedOverTime.runSpeedAlignedRun(trIndOpt,:));
        speedOverTimeOpt.SEMRun = std(speedOverTime.runSpeedAlignedRun(trIndOpt,:))/sqrt(nTrIndOpt);
    
        speedOverTimeOpt.meanRunSpeedBL = mean(speedOverTime.runSpeedAlignedRun(trIndOpt,speedOverTime.indBaseline),2);
        speedOverTimeOpt.meanRunSpeedBefRun = mean(speedOverTime.runSpeedAlignedRun(trIndOpt,speedOverTime.indBefRun),2);
        speedOverTimeOpt.meanRunSpeed0to1 = mean(speedOverTime.runSpeedAlignedRun(trIndOpt,speedOverTime.ind0to1),2);
        speedOverTimeOpt.meanRunSpeed3to5 = mean(speedOverTime.runSpeedAlignedRun(trIndOpt,speedOverTime.ind3to5),2);
        
        nTrIndOptCtrl = length(trIndOptCtrl); %% changed by Yingxue on 1/25/2022 (SEM was caculated using nTr before)
        speedOverTimeOptCtrl.trInd = trIndOptCtrl;
        speedOverTimeOptCtrl.meanRun = mean(speedOverTime.runSpeedAlignedRun(trIndOptCtrl,:));
        speedOverTimeOptCtrl.stdRun = std(speedOverTime.runSpeedAlignedRun(trIndOptCtrl,:));
        speedOverTimeOptCtrl.SEMRun = std(speedOverTime.runSpeedAlignedRun(trIndOptCtrl,:))/sqrt(nTrIndOptCtrl);
        
        speedOverTimeOptCtrl.meanRunSpeedBL = mean(speedOverTime.runSpeedAlignedRun(trIndOptCtrl,speedOverTime.indBaseline),2);
        speedOverTimeOptCtrl.meanRunSpeedBefRun = mean(speedOverTime.runSpeedAlignedRun(trIndOptCtrl,speedOverTime.indBefRun),2);
        speedOverTimeOptCtrl.meanRunSpeed0to1 = mean(speedOverTime.runSpeedAlignedRun(trIndOptCtrl,speedOverTime.ind0to1),2);
        speedOverTimeOptCtrl.meanRunSpeed3to5 = mean(speedOverTime.runSpeedAlignedRun(trIndOptCtrl,speedOverTime.ind3to5),2);
        
        %% added by Yingxue on 1/25/2022
        nTrIndOptCtrl1 = length(trIndOptCtrl1); 
        speedOverTimeOptCtrl1.trInd = trIndOptCtrl1;
        speedOverTimeOptCtrl1.meanRun = mean(speedOverTime.runSpeedAlignedRun(trIndOptCtrl1,:));
        speedOverTimeOptCtrl1.stdRun = std(speedOverTime.runSpeedAlignedRun(trIndOptCtrl1,:));
        speedOverTimeOptCtrl1.SEMRun = std(speedOverTime.runSpeedAlignedRun(trIndOptCtrl1,:))/sqrt(nTrIndOptCtrl1);
        
        speedOverTimeOptCtrl1.meanRunSpeedBL = mean(speedOverTime.runSpeedAlignedRun(trIndOptCtrl1,speedOverTime.indBaseline),2);
        speedOverTimeOptCtrl1.meanRunSpeedBefRun = mean(speedOverTime.runSpeedAlignedRun(trIndOptCtrl1,speedOverTime.indBefRun),2);
        speedOverTimeOptCtrl1.meanRunSpeed0to1 = mean(speedOverTime.runSpeedAlignedRun(trIndOptCtrl1,speedOverTime.ind0to1),2);
        speedOverTimeOptCtrl1.meanRunSpeed3to5 = mean(speedOverTime.runSpeedAlignedRun(trIndOptCtrl1,speedOverTime.ind3to5),2);
    else
        speedOverTimeOpt.trInd = [];
        speedOverTimeOpt.meanRun = [];
        speedOverTimeOpt.stdRun = [];
        speedOverTimeOpt.SEMRun = [];
        
        speedOverTimeOpt.meanRunSpeedBL = [];
        speedOverTimeOpt.meanRunSpeedBefRun = [];
        speedOverTimeOpt.meanRunSpeed0to1 = [];
        speedOverTimeOpt.meanRunSpeed3to5 = [];
        
        speedOverTimeOptCtrl.trInd = [];
        speedOverTimeOptCtrl.meanRun = [];
        speedOverTimeOptCtrl.stdRun = [];
        speedOverTimeOptCtrl.SEMRun = [];
        
        speedOverTimeOptCtrl.meanRunSpeedBL = [];
        speedOverTimeOptCtrl.meanRunSpeedBefRun = [];
        speedOverTimeOptCtrl.meanRunSpeed0to1 = [];
        speedOverTimeOptCtrl.meanRunSpeed3to5 = [];
        
        %% added by Yingxue on 1/25/2022
        speedOverTimeOptCtrl1.trInd = [];
        speedOverTimeOptCtrl1.meanRun = [];
        speedOverTimeOptCtrl1.stdRun = [];
        speedOverTimeOptCtrl1.SEMRun = [];
        
        speedOverTimeOptCtrl1.meanRunSpeedBL = [];
        speedOverTimeOptCtrl1.meanRunSpeedBefRun = [];
        speedOverTimeOptCtrl1.meanRunSpeed0to1 = [];
        speedOverTimeOptCtrl1.meanRunSpeed3to5 = [];
    end
            
    save([path fileName '_runSpeedTime_Run' num2str(onlyRun) '_msess' num2str(mazeSess) '.mat'],...
        'speedOverTime','speedOverTimeOpt','speedOverTimeOptCtrl','speedOverTimeOptCtrl1','paramC');
end

function runSpeed = speedAlignedRun(trialsRun,trialNo,trialLen)
    runSpeed = zeros(length(trialNo),trialLen);
    for i = 1:length(trialNo)
        tmpSpeed = [trialsRun.speed_MMsecBef{trialNo(i)}',...
            trialsRun.speed_MMsec{trialNo(i)}'];
        tmpSpeed(tmpSpeed < 0) = 0; % changed 9/12/2020 by Yingxue
        if(length(tmpSpeed) > trialLen)
            runSpeed(i,:) = tmpSpeed(1:trialLen);
        else
            runSpeed(i,1:length(tmpSpeed)) = tmpSpeed;
        end
    end
end

function runSpeed = speedAlignedRew(trialsRew,trialNo,trialLen)
    runSpeed = zeros(length(trialNo),trialLen);
    for i = 1:length(trialNo)
        tmpSpeed = trialsRew.speed_MMsec{trialNo(i)};
        tmpSpeed(tmpSpeed < 0) = 0;
        if(length(tmpSpeed) > trialLen)
            runSpeed(i,:) = tmpSpeed(1:trialLen);
        else
            runSpeed(i,1:length(tmpSpeed)) = tmpSpeed;
        end
    end
end