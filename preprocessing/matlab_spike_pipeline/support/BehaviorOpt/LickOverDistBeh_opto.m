function LickOverDistBeh_opto(path, fileName, mazeSess)

    fullPath = [path fileName '_alignRun_msess' num2str(mazeSess) '.mat']; 
    if(exist(fullPath) == 0)
        disp('The aligned to run file does not exist');
        disp(fullPath);
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
    load(fullPath,'param','behParCue');
        
    GlobalConst;
    tracks = 2200;
    spaceMergeBin = 10; %mm
    if(spaceMergeBin ~= 0)
        param.spaceSteps = [-spaceMergeBin/2:spaceMergeBin:tracks+spaceMergeBin/2];
    else
        param.spaceSteps = [0:tracks];
    end
    
    nTrials = length(trialsRun.goodTrial);
    trInd = param.startTr:param.endTr;
    %% removing the grooming trial. Added by Yingxue on 4/22/2021
    indGroomTr = max(behParCue.numSamples) > 30*sampleFq; % a trial longer than 30 sec
    if(sum(indGroomTr) > 0)
        [~,indGroomTr] = max(behParCue.numSamples);
        if(length(indGroomTr) > 1)
            a = 1;
        end
        trInd = setdiff(trInd,indGroomTr(1));
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

    numBins = length(param.spaceSteps);
    step = spaceMergeBin;
    lickOverDist.Run = zeros(nTrials,numBins);
    lickOverDist.Rew = zeros(nTrials,numBins);
    lickOverDist.Cue = zeros(nTrials,numBins);
    for tr = 1:nTrials
        %% aligned to run onset
        licks = trialsRun.lickLfpInd{tr};
        if(~isempty(licks))
            licks = trialsRun.xMM{tr}(licks);
        else
            continue;
        end
        lickTmp = hist(licks,param.spaceSteps);
        lickOverDist.Run(tr,:) = lickTmp;
        
        %% aligned to reward
        licks = trialsRew.lickLfpInd{tr};
        if(~isempty(licks))
            licks = trialsRew.xMM{tr}(licks);
        else
            continue;
        end
        lickTmp = hist(licks,param.spaceSteps);
        lickOverDist.Rew(tr,:) = lickTmp;
        
        %% aligned to cue
        licks = trialsCue.lickLfpInd{tr};
        if(~isempty(licks))
            licks = trialsCue.xMM{tr}(licks);
        else
            continue;
        end
        lickTmp = hist(licks,param.spaceSteps);
        lickOverDist.Cue(tr,:) = lickTmp;
    end
    
    %%
    lickOverDist.trInd = trInd;
    lickOverDist.meanRun = mean(lickOverDist.Run(trInd,:));
    lickOverDist.stdRun = std(lickOverDist.Run(trInd,:));
    lickOverDist.SEMRun = std(lickOverDist.Run(trInd,:))/sqrt(nTr);
    
    lickOverDist.ind30to100 = find(param.spaceSteps >= 300 & param.spaceSteps < 1000);
    lickOverDist.meanRun30to100 = mean(lickOverDist.Run(trInd,lickOverDist.ind30to100),2);
    lickOverDist.stdRun30to100 = std(lickOverDist.Run(trInd,lickOverDist.ind30to100),[],2);
    lickOverDist.SEMRun30to100 = std(lickOverDist.Run(trInd,lickOverDist.ind30to100),[],2)/sqrt(nTr);
    
    lickOverDist.ind30to120 = find(param.spaceSteps >= 300 & param.spaceSteps < 1200);
    lickOverDist.meanRun30to120 = mean(lickOverDist.Run(trInd,lickOverDist.ind30to120),2);
    lickOverDist.stdRun30to120 = std(lickOverDist.Run(trInd,lickOverDist.ind30to120),[],2);
    lickOverDist.SEMRun30to120 = std(lickOverDist.Run(trInd,lickOverDist.ind30to120),[],2)/sqrt(nTr);
    
    lickOverDist.ind100to150 = find(param.spaceSteps >= 1000 & param.spaceSteps < 1500);
    lickOverDist.meanRun100to150 = mean(lickOverDist.Run(trInd,lickOverDist.ind100to150),2);
    lickOverDist.stdRun100to150 = std(lickOverDist.Run(trInd,lickOverDist.ind100to150),[],2);
    lickOverDist.SEMRun100to150 = std(lickOverDist.Run(trInd,lickOverDist.ind100to150),[],2)/sqrt(nTr);
    
    lickOverDist.ind150to180 = find(param.spaceSteps >= 1500 & param.spaceSteps < 1800);
    lickOverDist.meanRun150to180 = mean(lickOverDist.Run(trInd,lickOverDist.ind150to180),2);
    lickOverDist.stdRun150to180 = std(lickOverDist.Run(trInd,lickOverDist.ind150to180),[],2);
    lickOverDist.SEMRun150to180 = std(lickOverDist.Run(trInd,lickOverDist.ind150to180),[],2)/sqrt(nTr);
    
    lickOverDist.ind180to210 = find(param.spaceSteps >= 1800 & param.spaceSteps < 2100);
    lickOverDist.meanRun180to210 = mean(lickOverDist.Run(trInd,lickOverDist.ind180to210),2);
    lickOverDist.stdRun180to210 = std(lickOverDist.Run(trInd,lickOverDist.ind180to210),[],2);
    lickOverDist.SEMRun180to210 = std(lickOverDist.Run(trInd,lickOverDist.ind180to210),[],2)/sqrt(nTr);
    
    if(~isempty(trIndOpt))
        nTrOpt = length(trIndOpt); %% changed by Yingxue on 1/25/2022 (SEM was caculated using nTr before)
        lickOverDistOpt.trInd = trIndOpt;
        lickOverDistOpt.meanRun = mean(lickOverDist.Run(trIndOpt,:));
        lickOverDistOpt.stdRun = std(lickOverDist.Run(trIndOpt,:));
        lickOverDistOpt.SEMRun = std(lickOverDist.Run(trIndOpt,:))/sqrt(nTrOpt);
        
        lickOverDistOpt.meanRun30to100 = mean(lickOverDist.Run(trIndOpt,lickOverDist.ind30to100),2);
        lickOverDistOpt.stdRun30to100 = std(lickOverDist.Run(trIndOpt,lickOverDist.ind30to100),[],2);
        lickOverDistOpt.SEMRun30to100 = std(lickOverDist.Run(trIndOpt,lickOverDist.ind30to100),[],2)/sqrt(nTrOpt);
        
        lickOverDistOpt.meanRun30to120 = mean(lickOverDist.Run(trIndOpt,lickOverDist.ind30to120),2);
        lickOverDistOpt.stdRun30to120 = std(lickOverDist.Run(trIndOpt,lickOverDist.ind30to120),[],2);
        lickOverDistOpt.SEMRun30to120 = std(lickOverDist.Run(trIndOpt,lickOverDist.ind30to120),[],2)/sqrt(nTrOpt);

        lickOverDistOpt.meanRun100to150 = mean(lickOverDist.Run(trIndOpt,lickOverDist.ind100to150),2);
        lickOverDistOpt.stdRun100to150 = std(lickOverDist.Run(trIndOpt,lickOverDist.ind100to150),[],2);
        lickOverDistOpt.SEMRun100to150 = std(lickOverDist.Run(trIndOpt,lickOverDist.ind100to150),[],2)/sqrt(nTrOpt);

        lickOverDistOpt.meanRun150to180 = mean(lickOverDist.Run(trIndOpt,lickOverDist.ind150to180),2);
        lickOverDistOpt.stdRun150to180 = std(lickOverDist.Run(trIndOpt,lickOverDist.ind150to180),[],2);
        lickOverDistOpt.SEMRun150to180 = std(lickOverDist.Run(trIndOpt,lickOverDist.ind150to180),[],2)/sqrt(nTrOpt);

        lickOverDistOpt.meanRun180to210 = mean(lickOverDist.Run(trIndOpt,lickOverDist.ind180to210),2);
        lickOverDistOpt.stdRun180to210 = std(lickOverDist.Run(trIndOpt,lickOverDist.ind180to210),[],2);
        lickOverDistOpt.SEMRun180to210 = std(lickOverDist.Run(trIndOpt,lickOverDist.ind180to210),[],2)/sqrt(nTrOpt);
        
        nTrOptCtrl = length(trIndOptCtrl); %% changed by Yingxue on 1/25/2022 (SEM was caculated using nTr before)
        lickOverDistOptCtrl.trInd = trIndOptCtrl;        
        lickOverDistOptCtrl.meanRun = mean(lickOverDist.Run(trIndOptCtrl,:));
        lickOverDistOptCtrl.stdRun = std(lickOverDist.Run(trIndOptCtrl,:));
        lickOverDistOptCtrl.SEMRun = std(lickOverDist.Run(trIndOptCtrl,:))/sqrt(nTrOptCtrl);
        
        lickOverDistOptCtrl.meanRun30to100 = mean(lickOverDist.Run(trIndOptCtrl,lickOverDist.ind30to100),2);
        lickOverDistOptCtrl.stdRun30to100 = std(lickOverDist.Run(trIndOptCtrl,lickOverDist.ind30to100),[],2);
        lickOverDistOptCtrl.SEMRun30to100 = std(lickOverDist.Run(trIndOptCtrl,lickOverDist.ind30to100),[],2)/sqrt(nTrOptCtrl);
        
        lickOverDistOptCtrl.meanRun30to120 = mean(lickOverDist.Run(trIndOptCtrl,lickOverDist.ind30to120),2);
        lickOverDistOptCtrl.stdRun30to120 = std(lickOverDist.Run(trIndOptCtrl,lickOverDist.ind30to120),[],2);
        lickOverDistOptCtrl.SEMRun30to120 = std(lickOverDist.Run(trIndOptCtrl,lickOverDist.ind30to120),[],2)/sqrt(nTrOptCtrl);

        lickOverDistOptCtrl.meanRun100to150 = mean(lickOverDist.Run(trIndOptCtrl,lickOverDist.ind100to150),2);
        lickOverDistOptCtrl.stdRun100to150 = std(lickOverDist.Run(trIndOptCtrl,lickOverDist.ind100to150),[],2);
        lickOverDistOptCtrl.SEMRun100to150 = std(lickOverDist.Run(trIndOptCtrl,lickOverDist.ind100to150),[],2)/sqrt(nTrOptCtrl);

        lickOverDistOptCtrl.meanRun150to180 = mean(lickOverDist.Run(trIndOptCtrl,lickOverDist.ind150to180),2);
        lickOverDistOptCtrl.stdRun150to180 = std(lickOverDist.Run(trIndOptCtrl,lickOverDist.ind150to180),[],2);
        lickOverDistOptCtrl.SEMRun150to180 = std(lickOverDist.Run(trIndOptCtrl,lickOverDist.ind150to180),[],2)/sqrt(nTrOptCtrl);

        lickOverDistOptCtrl.meanRun180to210 = mean(lickOverDist.Run(trIndOptCtrl,lickOverDist.ind180to210),2);
        lickOverDistOptCtrl.stdRun180to210 = std(lickOverDist.Run(trIndOptCtrl,lickOverDist.ind180to210),[],2);
        lickOverDistOptCtrl.SEMRun180to210 = std(lickOverDist.Run(trIndOptCtrl,lickOverDist.ind180to210),[],2)/sqrt(nTrOptCtrl);
        
        %% added by Yingxue on 1/25/2022
        nTrOptCtrl1 = length(trIndOptCtrl1);
        lickOverDistOptCtrl1.trInd = trIndOptCtrl1;        
        lickOverDistOptCtrl1.meanRun = mean(lickOverDist.Run(trIndOptCtrl1,:));
        lickOverDistOptCtrl1.stdRun = std(lickOverDist.Run(trIndOptCtrl1,:));
        lickOverDistOptCtrl1.SEMRun = std(lickOverDist.Run(trIndOptCtrl1,:))/sqrt(nTrOptCtrl1);
        
        lickOverDistOptCtrl1.meanRun30to100 = mean(lickOverDist.Run(trIndOptCtrl1,lickOverDist.ind30to100),2);
        lickOverDistOptCtrl1.stdRun30to100 = std(lickOverDist.Run(trIndOptCtrl1,lickOverDist.ind30to100),[],2);
        lickOverDistOptCtrl1.SEMRun30to100 = std(lickOverDist.Run(trIndOptCtrl1,lickOverDist.ind30to100),[],2)/sqrt(nTrOptCtrl1);
        
        lickOverDistOptCtrl1.meanRun30to120 = mean(lickOverDist.Run(trIndOptCtrl1,lickOverDist.ind30to120),2);
        lickOverDistOptCtrl1.stdRun30to120 = std(lickOverDist.Run(trIndOptCtrl1,lickOverDist.ind30to120),[],2);
        lickOverDistOptCtrl1.SEMRun30to120 = std(lickOverDist.Run(trIndOptCtrl1,lickOverDist.ind30to120),[],2)/sqrt(nTrOptCtrl1);

        lickOverDistOptCtrl1.meanRun100to150 = mean(lickOverDist.Run(trIndOptCtrl1,lickOverDist.ind100to150),2);
        lickOverDistOptCtrl1.stdRun100to150 = std(lickOverDist.Run(trIndOptCtrl1,lickOverDist.ind100to150),[],2);
        lickOverDistOptCtrl1.SEMRun100to150 = std(lickOverDist.Run(trIndOptCtrl1,lickOverDist.ind100to150),[],2)/sqrt(nTrOptCtrl1);

        lickOverDistOptCtrl1.meanRun150to180 = mean(lickOverDist.Run(trIndOptCtrl1,lickOverDist.ind150to180),2);
        lickOverDistOptCtrl1.stdRun150to180 = std(lickOverDist.Run(trIndOptCtrl1,lickOverDist.ind150to180),[],2);
        lickOverDistOptCtrl1.SEMRun150to180 = std(lickOverDist.Run(trIndOptCtrl1,lickOverDist.ind150to180),[],2)/sqrt(nTrOptCtrl1);

        lickOverDistOptCtrl1.meanRun180to210 = mean(lickOverDist.Run(trIndOptCtrl1,lickOverDist.ind180to210),2);
        lickOverDistOptCtrl1.stdRun180to210 = std(lickOverDist.Run(trIndOptCtrl1,lickOverDist.ind180to210),[],2);
        lickOverDistOptCtrl1.SEMRun180to210 = std(lickOverDist.Run(trIndOptCtrl1,lickOverDist.ind180to210),[],2)/sqrt(nTrOptCtrl1);
    else
        lickOverDistOpt.trInd = [];
        lickOverDistOpt.meanRun = [];
        lickOverDistOpt.stdRun = [];
        lickOverDistOpt.SEMRun = [];
        
        lickOverDistOpt.meanRun30to100 = [];
        lickOverDistOpt.stdRun30to100 = [];
        lickOverDistOpt.SEMRun30to100 = [];
        
        lickOverDistOpt.meanRun30to120 = [];
        lickOverDistOpt.stdRun30to120 = [];
        lickOverDistOpt.SEMRun30to120 = [];

        lickOverDistOpt.meanRun100to150 = [];
        lickOverDistOpt.stdRun100to150 = [];
        lickOverDistOpt.SEMRun100to150 = [];

        lickOverDistOpt.meanRun150to180 = [];
        lickOverDistOpt.stdRun150to180 = [];
        lickOverDistOpt.SEMRun150to180 = [];

        lickOverDistOpt.meanRun180to210 = [];
        lickOverDistOpt.stdRun180to210 = [];
        lickOverDistOpt.SEMRun180to210 = [];
        
        lickOverDistOptCtrl.trInd = [];        
        lickOverDistOptCtrl.meanRun = [];
        lickOverDistOptCtrl.stdRun = [];
        lickOverDistOptCtrl.SEMRun = [];
        
        lickOverDistOptCtrl.meanRun30to100 = [];
        lickOverDistOptCtrl.stdRun30to100 = [];
        lickOverDistOptCtrl.SEMRun30to100 = [];
        
        lickOverDistOptCtrl.meanRun30to120 = [];
        lickOverDistOptCtrl.stdRun30to120 = [];
        lickOverDistOptCtrl.SEMRun30to120 = [];

        lickOverDistOptCtrl.meanRun100to150 = [];
        lickOverDistOptCtrl.stdRun100to150 = [];
        lickOverDistOptCtrl.SEMRun100to150 = [];

        lickOverDistOptCtrl.meanRun150to180 = [];
        lickOverDistOptCtrl.stdRun150to180 = [];
        lickOverDistOptCtrl.SEMRun150to180 = [];

        lickOverDistOptCtrl.meanRun180to210 = [];
        lickOverDistOptCtrl.stdRun180to210 = [];
        lickOverDistOptCtrl.SEMRun180to210 = [];
        
        %% added by Yingxue on 1/25/2022
        lickOverDistOptCtrl1.trInd = [];        
        lickOverDistOptCtrl1.meanRun = [];
        lickOverDistOptCtrl1.stdRun = [];
        lickOverDistOptCtrl1.SEMRun = [];
        
        lickOverDistOptCtrl1.meanRun30to100 = [];
        lickOverDistOptCtrl1.stdRun30to100 = [];
        lickOverDistOptCtrl1.SEMRun30to100 = [];
        
        lickOverDistOptCtrl1.meanRun30to120 = [];
        lickOverDistOptCtrl1.stdRun30to120 = [];
        lickOverDistOptCtrl1.SEMRun30to120 = [];

        lickOverDistOptCtrl1.meanRun100to150 = [];
        lickOverDistOptCtrl1.stdRun100to150 = [];
        lickOverDistOptCtrl1.SEMRun100to150 = [];

        lickOverDistOptCtrl1.meanRun150to180 = [];
        lickOverDistOptCtrl1.stdRun150to180 = [];
        lickOverDistOptCtrl1.SEMRun150to180 = [];

        lickOverDistOptCtrl1.meanRun180to210 = [];
        lickOverDistOptCtrl1.stdRun180to210 = [];
        lickOverDistOptCtrl1.SEMRun180to210 = [];
    end
        
    lickOverDist.meanRew = mean(lickOverDist.Rew(trInd,:));
    lickOverDist.stdRew = std(lickOverDist.Rew(trInd,:));
    lickOverDist.SEMRew = std(lickOverDist.Rew(trInd,:))/sqrt(nTr);
    
    lickOverDist.meanCue = mean(lickOverDist.Cue(trInd,:));
    lickOverDist.stdCue = std(lickOverDist.Cue(trInd,:));
    lickOverDist.SEMCue = std(lickOverDist.Cue(trInd,:))/sqrt(nTr);
    
    save([path fileName '_lickDist_msess' num2str(mazeSess) '.mat'],...
        'lickOverDist','lickOverDistOpt','lickOverDistOptCtrl',...
        'lickOverDistOptCtrl1','param');
end