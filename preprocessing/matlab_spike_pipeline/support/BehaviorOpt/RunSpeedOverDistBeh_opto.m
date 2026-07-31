function RunSpeedOverDistBeh_opto(path, fileName, onlyRun, mazeSess)

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
    
    fullPath = [path fileName '_behPar_msess' num2str(mazeSess) '.mat'];
    if(exist(fullPath) == 0)
        disp('The behPar file does not exist');
        return;
    end
    load(fullPath,'param','behParCue','behParRun');
        
    GlobalConst;
    tracks = 1800;
    if(spaceMergeBin ~= 0)
        param.spaceSteps = [0:spaceMergeBin:tracks];
    else
        param.spaceSteps = [0:tracks];
    end
    
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
    
%     %% new code
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
    
%%    old code   
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
    
    numBins = length(param.spaceSteps);
    step = spaceMergeBin;
    timePerBin.Run = zeros(nTrials,numBins);
    speedOverDist.Run = zeros(nTrials,numBins);
    timePerBin.Rew = zeros(nTrials,numBins);
    speedOverDist.Rew = zeros(nTrials,numBins);
    timePerBin.Cue = zeros(nTrials,numBins);
    speedOverDist.Cue = zeros(nTrials,numBins);
    for tr = 1:nTrials
        %% aligned to run onset
        speedTmp = zeros(1,numBins);
        timePerBinTmp = zeros(1,numBins);
        timePerBinLSpTmp = zeros(1,numBins);
        lowSpeedTmp = zeros(1,numBins);
        if(onlyRun == 1)
            indSpeed = find(trialsRun.speed_MMsec{tr} > minSpeed);
            distRun = trialsRun.xMM{tr}(indSpeed);
        else
            distRun = trialsRun.xMM{tr};
            indSpeed = 1:length(trialsRun.xMM{tr});
        end
        for i = 1:numBins
            ind = find(distRun >= param.spaceSteps(i)-step/2 & distRun < param.spaceSteps(i)+step/2);
            indOrig = indSpeed(ind);
            if(~isempty(ind))
                timePerBinTmp(i) = length(ind);
                speedTmp(i) = mean(trialsRun.speed_MMsec{tr}(indOrig));
            else
                timePerBinTmp(i) = 1;
                if(i > 1)
                    speedTmp(i) = speedTmp(i-1);
                end
            end
            
            %% added on 8/13/2022
            ind = find(trialsRun.xMM{tr} >= param.spaceSteps(i)-step/2 & trialsRun.xMM{tr} < param.spaceSteps(i)+step/2);
            indLowSpeed = behParRun.lowSpeed{tr}(ind);
            if(~isempty(ind))
                timePerBinLSpTmp(i) = length(ind);
                lowSpeedTmp(i) = sum(indLowSpeed);
            else
                timePerBinLSpTmp(i) = 1;
                if(i > 1)
                    lowSpeedTmp(i) = lowSpeedTmp(i-1);
                end
            end
        end
        speedOverDist.Run(tr,:) = speedTmp;
        timePerBin.Run(tr,:) = timePerBinTmp;
        speedOverDist.lowSpeed(tr,:) = lowSpeedTmp; %% added on 8/13/2022
        timePerBin.LowSpeed(tr,:) = timePerBinLSpTmp; %% added on 8/13/2022
        
        %% aligned to reward
        speedTmp = zeros(1,numBins);
        timePerBinTmp = zeros(1,numBins);
        if(onlyRun == 1)
            indSpeed = find(trialsRew.speed_MMsec{tr} > minSpeed);
            distRun = trialsRew.xMM{tr}(indSpeed);
        else
            distRun = trialsRew.xMM{tr};
            indSpeed = 1:length(trialsRew.xMM{tr});
        end
        for i = 1:numBins
            ind = find(distRun >= param.spaceSteps(i)-step/2 & distRun < param.spaceSteps(i)+step/2);
            indOrig = indSpeed(ind);
            if(~isempty(ind))
                timePerBinTmp(i) = length(ind);
                speedTmp(i) = mean(trialsRew.speed_MMsec{tr}(indOrig));
            else
                timePerBinTmp(i) = 1;
                if(i > 1)
                    speedTmp(i) = speedTmp(i-1);
                end
            end
        end
        speedOverDist.Rew(tr,:) = speedTmp;
        timePerBin.Rew(tr,:) = timePerBinTmp;
        
        %% aligned to cue
        speedTmp = zeros(1,numBins);
        timePerBinTmp = zeros(1,numBins);
        if(onlyRun == 1)
            indSpeed = find(trialsCue.speed_MMsec{tr} > minSpeed);
            distRun = trialsCue.xMM{tr}(indSpeed);
        else
            distRun = trialsCue.xMM{tr};
            indSpeed = 1:length(trialsCue.xMM{tr});
        end
        for i = 1:numBins
            ind = find(distRun >= param.spaceSteps(i)-step/2 & distRun < param.spaceSteps(i)+step/2);
            indOrig = indSpeed(ind);
            if(~isempty(ind))
                timePerBinTmp(i) = length(ind);
                speedTmp(i) = mean(trialsCue.speed_MMsec{tr}(indOrig));
            else
                timePerBinTmp(i) = 1;
                if(i > 1)
                    speedTmp(i) = speedTmp(i-1);
                end
            end
        end
        speedOverDist.Cue(tr,:) = speedTmp;
        timePerBin.Cue(tr,:) = timePerBinTmp;
    end
    
    speedOverDist.trInd = trInd;
    speedOverDist.meanRun = mean(speedOverDist.Run(trInd,:));
    speedOverDist.stdRun = std(speedOverDist.Run(trInd,:));
    speedOverDist.SEMRun = std(speedOverDist.Run(trInd,:))/sqrt(nTr);
    
    %% added on 8/13/2022
    speedOverDist.meanRunLowSpeed = mean(speedOverDist.lowSpeed(trInd,:));
    speedOverDist.stdRunLowSpeed = std(speedOverDist.lowSpeed(trInd,:));
    speedOverDist.SEMRunLowSpeed = std(speedOverDist.lowSpeed(trInd,:))/sqrt(nTr);
    %%
    
    speedOverDist.ind30to100 = find(param.spaceSteps >= 300 & param.spaceSteps < 1000);
    speedOverDist.meanRun30to100 = mean(speedOverDist.Run(trInd,speedOverDist.ind30to100),2);
    speedOverDist.stdRun30to100 = std(speedOverDist.Run(trInd,speedOverDist.ind30to100),[],2);
    speedOverDist.SEMRun30to100 = std(speedOverDist.Run(trInd,speedOverDist.ind30to100),[],2)/sqrt(nTr);
    
    speedOverDist.indAfter100 = find(param.spaceSteps >= 1000);
    speedOverDist.meanRunAfter100 = mean(speedOverDist.Run(trInd,speedOverDist.indAfter100),2);
    speedOverDist.stdRunAfter100 = std(speedOverDist.Run(trInd,speedOverDist.indAfter100),[],2);
    speedOverDist.SEMRunAfter100 = std(speedOverDist.Run(trInd,speedOverDist.indAfter100),[],2)/sqrt(nTr);
    
    if(~isempty(trIndOpt))
        nTrIndOpt = length(trIndOpt); %% changed by Yingxue on 1/25/2022 (SEM was caculated using nTr before)
        speedOverDistOpt.trInd = trIndOpt;
        speedOverDistOpt.meanRun = mean(speedOverDist.Run(trIndOpt,:));
        speedOverDistOpt.stdRun = std(speedOverDist.Run(trIndOpt,:));
        speedOverDistOpt.SEMRun = std(speedOverDist.Run(trIndOpt,:))/sqrt(nTrIndOpt);
        
        %% added on 8/13/2022
        speedOverDistOpt.meanRunLowSpeed = mean(speedOverDist.lowSpeed(trIndOpt,:));
        speedOverDistOpt.stdRunLowSpeed = std(speedOverDist.lowSpeed(trIndOpt,:));
        speedOverDistOpt.SEMRunLowSpeed = std(speedOverDist.lowSpeed(trIndOpt,:))/sqrt(nTrIndOpt);
        %%
        
        speedOverDistOpt.meanRun30to100 = mean(speedOverDist.Run(trIndOpt,speedOverDist.ind30to100),2);
        speedOverDistOpt.stdRun30to100 = std(speedOverDist.Run(trIndOpt,speedOverDist.ind30to100),[],2);
        speedOverDistOpt.SEMRun30to100 = std(speedOverDist.Run(trIndOpt,speedOverDist.ind30to100),[],2)/sqrt(nTrIndOpt);
        
        speedOverDistOpt.meanRunAfter100 = mean(speedOverDist.Run(trIndOpt,speedOverDist.indAfter100),2);
        speedOverDistOpt.stdRunAfter100 = std(speedOverDist.Run(trIndOpt,speedOverDist.indAfter100),[],2);
        speedOverDistOpt.SEMRunAfter100 = std(speedOverDist.Run(trIndOpt,speedOverDist.indAfter100),[],2)/sqrt(nTrIndOpt);
        
        nTrIndOptCtrl = length(trIndOptCtrl); %% changed by Yingxue on 1/25/2022 (SEM was caculated using nTr before)
        speedOverDistOptCtrl.trInd = trIndOptCtrl;        
        speedOverDistOptCtrl.meanRun = mean(speedOverDist.Run(trIndOptCtrl,:));
        speedOverDistOptCtrl.stdRun = std(speedOverDist.Run(trIndOptCtrl,:));
        speedOverDistOptCtrl.SEMRun = std(speedOverDist.Run(trIndOptCtrl,:))/sqrt(nTrIndOptCtrl);
        
        %% added on 8/13/2022
        speedOverDistOptCtrl.meanRunLowSpeed = mean(speedOverDist.lowSpeed(trIndOptCtrl,:));
        speedOverDistOptCtrl.stdRunLowSpeed = std(speedOverDist.lowSpeed(trIndOptCtrl,:));
        speedOverDistOptCtrl.SEMRunLowSpeed = std(speedOverDist.lowSpeed(trIndOptCtrl,:))/sqrt(nTrIndOptCtrl);
        %%
        
        speedOverDistOptCtrl.meanRun30to100 = mean(speedOverDist.Run(trIndOptCtrl,speedOverDist.ind30to100),2);
        speedOverDistOptCtrl.stdRun30to100 = std(speedOverDist.Run(trIndOptCtrl,speedOverDist.ind30to100),[],2);
        speedOverDistOptCtrl.SEMRun30to100 = std(speedOverDist.Run(trIndOptCtrl,speedOverDist.ind30to100),[],2)/sqrt(nTrIndOptCtrl);
        
        speedOverDistOptCtrl.meanRunAfter100 = mean(speedOverDist.Run(trIndOptCtrl,speedOverDist.indAfter100),2);
        speedOverDistOptCtrl.stdRunAfter100 = std(speedOverDist.Run(trIndOptCtrl,speedOverDist.indAfter100),[],2);
        speedOverDistOptCtrl.SEMRunAfter100 = std(speedOverDist.Run(trIndOptCtrl,speedOverDist.indAfter100),[],2)/sqrt(nTrIndOptCtrl);
        
        %% added by Yingxue on 1/25/2022
        nTrIndOptCtrl1 = length(trIndOptCtrl1); 
        speedOverDistOptCtrl1.trInd = trIndOptCtrl1;        
        speedOverDistOptCtrl1.meanRun = mean(speedOverDist.Run(trIndOptCtrl1,:));
        speedOverDistOptCtrl1.stdRun = std(speedOverDist.Run(trIndOptCtrl1,:));
        speedOverDistOptCtrl1.SEMRun = std(speedOverDist.Run(trIndOptCtrl1,:))/sqrt(nTrIndOptCtrl1);
        
        %% added on 8/13/2022
        speedOverDistOptCtrl1.meanRunLowSpeed = mean(speedOverDist.lowSpeed(trIndOptCtrl1,:));
        speedOverDistOptCtrl1.stdRunLowSpeed = std(speedOverDist.lowSpeed(trIndOptCtrl1,:));
        speedOverDistOptCtrl1.SEMRunLowSpeed = std(speedOverDist.lowSpeed(trIndOptCtrl1,:))/sqrt(nTrIndOptCtrl1);
        %%
        
        speedOverDistOptCtrl1.meanRun30to100 = mean(speedOverDist.Run(trIndOptCtrl1,speedOverDist.ind30to100),2);
        speedOverDistOptCtrl1.stdRun30to100 = std(speedOverDist.Run(trIndOptCtrl1,speedOverDist.ind30to100),[],2);
        speedOverDistOptCtrl1.SEMRun30to100 = std(speedOverDist.Run(trIndOptCtrl1,speedOverDist.ind30to100),[],2)/sqrt(nTrIndOptCtrl1);
        
        speedOverDistOptCtrl1.meanRunAfter100 = mean(speedOverDist.Run(trIndOptCtrl1,speedOverDist.indAfter100),2);
        speedOverDistOptCtrl1.stdRunAfter100 = std(speedOverDist.Run(trIndOptCtrl1,speedOverDist.indAfter100),[],2);
        speedOverDistOptCtrl1.SEMRunAfter100 = std(speedOverDist.Run(trIndOptCtrl1,speedOverDist.indAfter100),[],2)/sqrt(nTrIndOptCtrl1);
    else
        speedOverDistOpt.trInd = [];
        speedOverDistOpt.meanRun = [];
        speedOverDistOpt.stdRun = [];
        speedOverDistOpt.SEMRun = [];
        
        %% added on 8/13/2022
        speedOverDistOpt.meanRunLowSpeed = [];
        speedOverDistOpt.stdRunLowSpeed = [];
        speedOverDistOpt.SEMRunLowSpeed = [];
        
        speedOverDistOpt.meanRun30to100 = [];
        speedOverDistOpt.stdRun30to100 = [];
        speedOverDistOpt.SEMRun30to100 = [];
        
        speedOverDistOpt.meanRunAfter100 = [];
        speedOverDistOpt.stdRunAfter100 = [];
        speedOverDistOpt.SEMRunAfter100 = [];
        
        speedOverDistOptCtrl.trInd = [];        
        speedOverDistOptCtrl.meanRun = [];
        speedOverDistOptCtrl.stdRun = [];
        speedOverDistOptCtrl.SEMRun = [];
        
        %% added on 8/13/2022
        speedOverDistOptCtrl.meanRunLowSpeed = [];
        speedOverDistOptCtrl.stdRunLowSpeed = [];
        speedOverDistOptCtrl.SEMRunLowSpeed = [];
        
        speedOverDistOptCtrl.meanRun30to100 = [];
        speedOverDistOptCtrl.stdRun30to100 = [];
        speedOverDistOptCtrl.SEMRun30to100 = [];
        
        speedOverDistOptCtrl.meanRunAfter100 = [];
        speedOverDistOptCtrl.stdRunAfter100 = [];
        speedOverDistOptCtrl.SEMRunAfter100 = [];
        
        %% added by Yingxue on 1/25/2022
        speedOverDistOptCtrl1.trInd = [];        
        speedOverDistOptCtrl1.meanRun = [];
        speedOverDistOptCtrl1.stdRun = [];
        speedOverDistOptCtrl1.SEMRun = [];
        
        %% added on 8/13/2022
        speedOverDistOptCtrl1.meanRunLowSpeed = [];
        speedOverDistOptCtrl1.stdRunLowSpeed = [];
        speedOverDistOptCtrl1.SEMRunLowSpeed = [];
                
        speedOverDistOptCtrl1.meanRun30to100 = [];
        speedOverDistOptCtrl1.stdRun30to100 = [];
        speedOverDistOptCtrl1.SEMRun30to100 = [];
        
        speedOverDistOptCtrl1.meanRunAfter100 = [];
        speedOverDistOptCtrl1.stdRunAfter100 = [];
        speedOverDistOptCtrl1.SEMRunAfter100 = [];
    end
    
    speedOverDist.meanRew = mean(speedOverDist.Rew(trInd,:));
    speedOverDist.stdRew = std(speedOverDist.Rew(trInd,:));
    speedOverDist.SEMRew = std(speedOverDist.Rew(trInd,:))/sqrt(nTr);
    
    speedOverDist.meanCue = mean(speedOverDist.Cue(trInd,:));
    speedOverDist.stdCue = std(speedOverDist.Cue(trInd,:));
    speedOverDist.SEMCue = std(speedOverDist.Cue(trInd,:))/sqrt(nTr);
    
    save([path fileName '_runSpeedDist_Run' num2str(onlyRun) '_msess' num2str(mazeSess) '.mat'],...
        'speedOverDist','speedOverDistOpt','speedOverDistOptCtrl','speedOverDistOptCtrl1',...
        'timePerBin','param');
end