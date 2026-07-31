function LickOverDistSim_opto(path,fileName,mazeSess)
% trial by trial correlation of licking 
    
    fullPath = [path fileName '_lickDist_msess' num2str(mazeSess) '.mat'];
    if(exist(fullPath) == 0)
        disp('The lick over distance file does not exist');
        disp(fullPath);
        return;
    end
    load(fullPath,'lickOverDist','lickOverDistOpt','lickOverDistOptCtrl','param');
    
    lickOverDistSim.Run = squareform(pdist(lickOverDist.Run,'cosine'));
    lickOverDistSim.Rew = squareform(pdist(lickOverDist.Rew,'cosine'));
    lickOverDistSim.Cue = squareform(pdist(lickOverDist.Cue,'cosine'));
    
    lickOverDistSim.EucRun = squareform(pdist(lickOverDist.Run,'euclidean'));
    lickOverDistSim.EucRew = squareform(pdist(lickOverDist.Rew,'euclidean'));
    lickOverDistSim.EucCue = squareform(pdist(lickOverDist.Cue,'euclidean'));
    
    %% Similarity
    lickOverDistSim.trInd = lickOverDist.trInd;
    nTrials = length(lickOverDist.trInd);
    nElem = nTrials*nTrials-nTrials;
    simArr = triu(lickOverDistSim.Run(lickOverDist.trInd,lickOverDist.trInd),1);
    simArr = simArr(:);
    lickOverDistSim.meanRun = sum(simArr(isnan(simArr) == 0))/nElem;
    lickOverDistSim.semRun = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
    
    simArr = triu(lickOverDistSim.Rew(lickOverDist.trInd,lickOverDist.trInd),1);
    simArr = simArr(:);
    lickOverDistSim.meanRew = sum(simArr(isnan(simArr) == 0))/nElem;
    lickOverDistSim.semRew = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
    
    simArr = triu(lickOverDistSim.Cue(lickOverDist.trInd,lickOverDist.trInd),1);
    simArr = simArr(:);
    lickOverDistSim.meanCue = sum(simArr(isnan(simArr) == 0))/nElem;
    lickOverDistSim.semCue = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
               
    if(~isempty(lickOverDistOpt.trInd))
        lickOverDistSimOpt.trInd = lickOverDistOpt.trInd;
        nTrials = length(lickOverDistOpt.trInd);
        nElem = nTrials*nTrials-nTrials;
        simArr = triu(lickOverDistSim.Run(lickOverDistOpt.trInd,lickOverDistOpt.trInd),1);
        simArr = simArr(:);
        lickOverDistSimOpt.meanRun = sum(simArr(isnan(simArr) == 0))/nElem;
        lickOverDistSimOpt.semRun = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
        
        simArr = triu(lickOverDistSim.Rew(lickOverDistOpt.trInd,lickOverDistOpt.trInd),1);
        simArr = simArr(:);
        lickOverDistSimOpt.meanRew = sum(simArr(isnan(simArr) == 0))/nElem;
        lickOverDistSimOpt.semRew = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
        
        simArr = triu(lickOverDistSim.Cue(lickOverDistOpt.trInd,lickOverDistOpt.trInd),1);
        simArr = simArr(:);
        lickOverDistSimOpt.meanCue = sum(simArr(isnan(simArr) == 0))/nElem;
        lickOverDistSimOpt.semCue = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
        
        lickOverDistSimOptCtrl.trInd = lickOverDistOptCtrl.trInd;
        nTrials = length(lickOverDistOptCtrl.trInd);
        nElem = nTrials*nTrials-nTrials;
        simArr = triu(lickOverDistSim.Run(lickOverDistOptCtrl.trInd,lickOverDistOptCtrl.trInd),1);
        simArr = simArr(:);
        lickOverDistSimOptCtrl.meanRun = sum(simArr(isnan(simArr) == 0))/nElem;
        lickOverDistSimOptCtrl.semRun = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
        
        simArr = triu(lickOverDistSim.Rew(lickOverDistOptCtrl.trInd,lickOverDistOptCtrl.trInd),1);
        simArr = simArr(:);
        lickOverDistSimOptCtrl.meanRew = sum(simArr(isnan(simArr) == 0))/nElem;
        lickOverDistSimOptCtrl.semRew = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
        
        simArr = triu(lickOverDistSim.Cue(lickOverDistOptCtrl.trInd,lickOverDistOptCtrl.trInd),1);
        simArr = simArr(:);
        lickOverDistSimOptCtrl.meanCue = sum(simArr(isnan(simArr) == 0))/nElem;
        lickOverDistSimOptCtrl.semCue = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
    else
        lickOverDistSimOpt.trInd = [];
        lickOverDistSimOpt.meanRun = [];
        lickOverDistSimOpt.meanRew = [];
        lickOverDistSimOpt.meanCue = [];
        lickOverDistSimOpt.semRun = [];
        lickOverDistSimOpt.semRew = [];
        lickOverDistSimOpt.semCue = [];
        
        lickOverDistSimOptCtrl.trInd = [];
        lickOverDistSimOptCtrl.meanRun = [];
        lickOverDistSimOptCtrl.meanRew = [];
        lickOverDistSimOptCtrl.meanCue = [];
        lickOverDistSimOptCtrl.semRun = [];
        lickOverDistSimOptCtrl.semRew = [];
        lickOverDistSimOptCtrl.semCue = [];
    end
    
    
    %% Euclidean distance
    nTrials = length(lickOverDist.trInd);
    nElem = nTrials*nTrials-nTrials;
    simArr = triu(lickOverDistSim.EucRun(lickOverDist.trInd,lickOverDist.trInd),1);
    simArr = simArr(:);
    lickOverDistSim.meanEucRun = sum(simArr(isnan(simArr) == 0))/nElem;
    lickOverDistSim.semEucRun = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
    
    simArr = triu(lickOverDistSim.EucRew(lickOverDist.trInd,lickOverDist.trInd),1);
    simArr = simArr(:);
    lickOverDistSim.meanEucRew = sum(simArr(isnan(simArr) == 0))/nElem;
    lickOverDistSim.semEucRew = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
    
    simArr = triu(lickOverDistSim.EucCue(lickOverDist.trInd,lickOverDist.trInd),1);
    simArr = simArr(:);
    lickOverDistSim.meanEucCue = sum(simArr(isnan(simArr) == 0))/nElem;
    lickOverDistSim.semEucCue = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
               
    if(~isempty(lickOverDistOpt.trInd))
        nTrials = length(lickOverDistOpt.trInd);
        nElem = nTrials*nTrials-nTrials;
        simArr = triu(lickOverDistSim.EucRun(lickOverDistOpt.trInd,lickOverDistOpt.trInd),1);
        simArr = simArr(:);
        lickOverDistSimOpt.meanEucRun = sum(simArr(isnan(simArr) == 0))/nElem;
        lickOverDistSimOpt.semEucRun = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
        
        simArr = triu(lickOverDistSim.EucRew(lickOverDistOpt.trInd,lickOverDistOpt.trInd),1);
        simArr = simArr(:);
        lickOverDistSimOpt.meanEucRew = sum(simArr(isnan(simArr) == 0))/nElem;
        lickOverDistSimOpt.semEucRew = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
        
        simArr = triu(lickOverDistSim.EucCue(lickOverDistOpt.trInd,lickOverDistOpt.trInd),1);
        simArr = simArr(:);
        lickOverDistSimOpt.meanEucCue = sum(simArr(isnan(simArr) == 0))/nElem;
        lickOverDistSimOpt.semEucCue = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
        
        nTrials = length(lickOverDistOptCtrl.trInd);
        nElem = nTrials*nTrials-nTrials;
        simArr = triu(lickOverDistSim.EucRun(lickOverDistOptCtrl.trInd,lickOverDistOptCtrl.trInd),1);
        simArr = simArr(:);
        lickOverDistSimOptCtrl.meanEucRun = sum(simArr(isnan(simArr) == 0))/nElem;
        lickOverDistSimOptCtrl.semEucRun = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
        
        simArr = triu(lickOverDistSim.EucRew(lickOverDistOptCtrl.trInd,lickOverDistOptCtrl.trInd),1);
        simArr = simArr(:);
        lickOverDistSimOptCtrl.meanEucRew = sum(simArr(isnan(simArr) == 0))/nElem;
        lickOverDistSimOptCtrl.semEucRew = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
        
        simArr = triu(lickOverDistSim.EucCue(lickOverDistOptCtrl.trInd,lickOverDistOptCtrl.trInd),1);
        simArr = simArr(:);
        lickOverDistSimOptCtrl.meanEucCue = sum(simArr(isnan(simArr) == 0))/nElem;
        lickOverDistSimOptCtrl.semEucCue = std(simArr(isnan(simArr) == 0))/sqrt(nElem);
    else
        lickOverDistSimOpt.meanEucRun = [];
        lickOverDistSimOpt.meanEucRew = [];
        lickOverDistSimOpt.meanEucCue = [];
        lickOverDistSimOpt.semEucRun = [];
        lickOverDistSimOpt.semEucRew = [];
        lickOverDistSimOpt.semEucCue = [];
        
        lickOverDistSimOptCtrl.meanEucRun = [];
        lickOverDistSimOptCtrl.meanEucRew = [];
        lickOverDistSimOptCtrl.meanEucCue = [];
        lickOverDistSimOptCtrl.semEucRun = [];
        lickOverDistSimOptCtrl.semEucRew = [];
        lickOverDistSimOptCtrl.semEucCue = [];
    end
            
    save([path fileName '_lickDistSim_msess' num2str(mazeSess) '.mat'],...
        'lickOverDistSim','lickOverDistSimOpt','lickOverDistSimOptCtrl','param');
end
