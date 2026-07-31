function RunSpeedOverDistSim_opto(path,fileName,onlyRun,mazeSess)
% trial by trial correlation of licking 
    
    fullPath = [path fileName '_runSpeedDist_Run' num2str(onlyRun) '_msess' num2str(mazeSess) '.mat'];
    if(exist(fullPath) == 0)
        disp('The speed over distance file does not exist');
        disp(fullPath);
        return;
    end
    load(fullPath,'speedOverDist','speedOverDistOpt','speedOverDistOptCtrl','param');
    
    speedOverDistSim.Run = squareform(pdist(speedOverDist.Run,'cosine'));
    speedOverDistSim.Rew = squareform(pdist(speedOverDist.Rew,'cosine'));
    speedOverDistSim.Cue = squareform(pdist(speedOverDist.Cue,'cosine'));
    
    speedOverDistSim.EucRun = squareform(pdist(speedOverDist.Run,'euclidean'));
    speedOverDistSim.EucRew = squareform(pdist(speedOverDist.Rew,'euclidean'));
    speedOverDistSim.EucCue = squareform(pdist(speedOverDist.Cue,'euclidean'));
    
    %% similarity
    speedOverDistSim.trInd = speedOverDist.trInd;
    nTrials = length(speedOverDist.trInd);
    nElem = nTrials*nTrials-nTrials;
    simArr = triu(speedOverDistSim.Run(speedOverDist.trInd,speedOverDist.trInd),1);
    simArr = simArr(:);
    speedOverDistSim.meanRun = sum(simArr(isnan(simArr) == 0))/nElem;
    speedOverDistSim.semRun = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
    
    simArr = triu(speedOverDistSim.Rew(speedOverDist.trInd,speedOverDist.trInd),1);
    simArr = simArr(:);
    speedOverDistSim.meanRew = sum(simArr(isnan(simArr) == 0))/nElem;
    speedOverDistSim.semRew = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
    
    simArr = triu(speedOverDistSim.Cue(speedOverDist.trInd,speedOverDist.trInd),1);
    simArr = simArr(:);
    speedOverDistSim.meanCue = sum(simArr(isnan(simArr) == 0))/nElem;
    speedOverDistSim.semCue = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
               
    if(~isempty(speedOverDistOpt.trInd))
        speedOverDistSimOpt.trInd = speedOverDistOpt.trInd;
        nTrials = length(speedOverDistOpt.trInd);
        nElem = nTrials*nTrials-nTrials;
        simArr = triu(speedOverDistSim.Run(speedOverDistOpt.trInd,speedOverDistOpt.trInd),1);
        simArr = simArr(:);
        speedOverDistSimOpt.meanRun = sum(simArr(isnan(simArr) == 0))/nElem;
        speedOverDistSimOpt.semRun = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
        
        simArr = triu(speedOverDistSim.Rew(speedOverDistOpt.trInd,speedOverDistOpt.trInd),1);
        simArr = simArr(:);
        speedOverDistSimOpt.meanRew = sum(simArr(isnan(simArr) == 0))/nElem;
        speedOverDistSimOpt.semRew = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
        
        simArr = triu(speedOverDistSim.Cue(speedOverDistOpt.trInd,speedOverDistOpt.trInd),1);
        simArr = simArr(:);
        speedOverDistSimOpt.meanCue = sum(simArr(isnan(simArr) == 0))/nElem;
        speedOverDistSimOpt.semCue = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
        
        speedOverDistSimOptCtrl.trInd = speedOverDistOptCtrl.trInd;
        nTrials = length(speedOverDistOptCtrl.trInd);
        nElem = nTrials*nTrials-nTrials;
        simArr = triu(speedOverDistSim.Run(speedOverDistOptCtrl.trInd,speedOverDistOptCtrl.trInd),1);
        simArr = simArr(:);
        speedOverDistSimOptCtrl.meanRun = sum(simArr(isnan(simArr) == 0))/nElem;
        speedOverDistSimOptCtrl.semRun = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
        
        simArr = triu(speedOverDistSim.Rew(speedOverDistOptCtrl.trInd,speedOverDistOptCtrl.trInd),1);
        simArr = simArr(:);
        speedOverDistSimOptCtrl.meanRew = sum(simArr(isnan(simArr) == 0))/nElem;
        speedOverDistSimOptCtrl.semRew = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
        
        simArr = triu(speedOverDistSim.Cue(speedOverDistOptCtrl.trInd,speedOverDistOptCtrl.trInd),1);
        simArr = simArr(:);
        speedOverDistSimOptCtrl.meanCue = sum(simArr(isnan(simArr) == 0))/nElem;
        speedOverDistSimOptCtrl.semCue = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
    else
        speedOverDistSimOpt.trInd = [];
        speedOverDistSimOpt.meanRun = [];
        speedOverDistSimOpt.meanRew = [];
        speedOverDistSimOpt.meanCue = [];
        speedOverDistSimOpt.semRun = [];
        speedOverDistSimOpt.semRew = [];
        speedOverDistSimOpt.semCue = [];
        
        speedOverDistSimOptCtrl.trInd = [];
        speedOverDistSimOptCtrl.meanRun = [];
        speedOverDistSimOptCtrl.meanRew = [];
        speedOverDistSimOptCtrl.meanCue = [];
        speedOverDistSimOptCtrl.semRun = [];
        speedOverDistSimOptCtrl.semRew = [];
        speedOverDistSimOptCtrl.semCue = [];
    end
    
    %% euclidean distance
    nTrials = length(speedOverDist.trInd);
    nElem = nTrials*nTrials-nTrials;
    simArr = triu(speedOverDistSim.EucRun(speedOverDist.trInd,speedOverDist.trInd),1);
    simArr = simArr(:);
    speedOverDistSim.meanEucRun = sum(simArr(isnan(simArr) == 0))/nElem;
    speedOverDistSim.semEucRun = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
    
    simArr = triu(speedOverDistSim.EucRew(speedOverDist.trInd,speedOverDist.trInd),1);
    simArr = simArr(:);
    speedOverDistSim.meanEucRew = sum(simArr(isnan(simArr) == 0))/nElem;
    speedOverDistSim.semEucRew = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
    
    simArr = triu(speedOverDistSim.EucCue(speedOverDist.trInd,speedOverDist.trInd),1);
    simArr = simArr(:);
    speedOverDistSim.meanEucCue = sum(simArr(isnan(simArr) == 0))/nElem;
    speedOverDistSim.semEucCue = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
               
    if(~isempty(speedOverDistOpt.trInd))
        nTrials = length(speedOverDistOpt.trInd);
        nElem = nTrials*nTrials-nTrials;
        simArr = triu(speedOverDistSim.EucRun(speedOverDistOpt.trInd,speedOverDistOpt.trInd),1);
        simArr = simArr(:);
        speedOverDistSimOpt.meanEucRun = sum(simArr(isnan(simArr) == 0))/nElem;
        speedOverDistSimOpt.semEucRun = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
        
        simArr = triu(speedOverDistSim.EucRew(speedOverDistOpt.trInd,speedOverDistOpt.trInd),1);
        simArr = simArr(:);
        speedOverDistSimOpt.meanEucRew = sum(simArr(isnan(simArr) == 0))/nElem;
        speedOverDistSimOpt.semEucRew = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
        
        simArr = triu(speedOverDistSim.EucCue(speedOverDistOpt.trInd,speedOverDistOpt.trInd),1);
        simArr = simArr(:);
        speedOverDistSimOpt.meanEucCue = sum(simArr(isnan(simArr) == 0))/nElem;
        speedOverDistSimOpt.semEucCue = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
        
        nTrials = length(speedOverDistOptCtrl.trInd);
        nElem = nTrials*nTrials-nTrials;
        simArr = triu(speedOverDistSim.EucRun(speedOverDistOptCtrl.trInd,speedOverDistOptCtrl.trInd),1);
        simArr = simArr(:);
        speedOverDistSimOptCtrl.meanEucRun = sum(simArr(isnan(simArr) == 0))/nElem;
        speedOverDistSimOptCtrl.semEucRun = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
        
        simArr = triu(speedOverDistSim.EucRew(speedOverDistOptCtrl.trInd,speedOverDistOptCtrl.trInd),1);
        simArr = simArr(:);
        speedOverDistSimOptCtrl.meanEucRew = sum(simArr(isnan(simArr) == 0))/nElem;
        speedOverDistSimOptCtrl.semEucRew = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
        
        simArr = triu(speedOverDistSim.EucCue(speedOverDistOptCtrl.trInd,speedOverDistOptCtrl.trInd),1);
        simArr = simArr(:);
        speedOverDistSimOptCtrl.meanEucCue = sum(simArr(isnan(simArr) == 0))/nElem;
        speedOverDistSimOptCtrl.semEucCue = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
    else
        speedOverDistSimOpt.meanEucRun = [];
        speedOverDistSimOpt.meanEucRew = [];
        speedOverDistSimOpt.meanEucCue = [];
        speedOverDistSimOpt.semEucRun = [];
        speedOverDistSimOpt.semEucRew = [];
        speedOverDistSimOpt.semEucCue = [];
        
        speedOverDistSimOptCtrl.meanEucRun = [];
        speedOverDistSimOptCtrl.meanEucRew = [];
        speedOverDistSimOptCtrl.meanEucCue = [];
        speedOverDistSimOptCtrl.semEucRun = [];
        speedOverDistSimOptCtrl.semEucRew = [];
        speedOverDistSimOptCtrl.semEucCue = [];
    end
            
    save([path fileName '_speedDistSim_msess' num2str(mazeSess) '_Run' num2str(onlyRun) '.mat'],...
        'speedOverDistSim','speedOverDistSimOpt','speedOverDistSimOptCtrl','param');
end
