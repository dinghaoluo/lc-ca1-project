function LickOverDistNoCue(path, fileName, mazeSess)

    
    %%%%%%%%% initialize constants
    fileNameInfo = [fileName '_Info.mat'];     
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo(path,fileName);
    end
    load(fullPath);
    indTr = find(beh.mazeSess == mazeSess);
    nTrials = length(indTr);
    
    totNTrial = length(beh.mazeSess);
    
    ind = strfind(fileName, '_');
    fileNameGen = [fileName(1:ind(1)) 'BehavElectrDataLFP.mat'];
    load([path fileNameGen],'Track','Laps');
    
    fileNameBeh = [fileName(1:ind(1)-1) 'BTDT.mat'];
    load([path fileNameBeh]);
%     
%     fileName = [fileName '.mat'];
%     fullPath = [path fileName];
%     if(exist(fullPath) == 0)
%         disp('Recording file does not exist.');
%         return;
%     end
%     load(fullPath,'trials');
        
    GlobalConst;
    tracks = 2200;
    trackStart = 300; % only for the aligned to pump case
    trackEnd = -1500;
    spaceMergeBin = 10; %mm
    if(spaceMergeBin ~= 0)
        param.spaceSteps = [-spaceMergeBin/2:spaceMergeBin:tracks+spaceMergeBin/2];
        param.spaceStepsAligned = [-spaceMergeBin/2-1500:spaceMergeBin...
            :spaceMergeBin/2+trackStart];
    else
        param.spaceSteps = [0:tracks];
        param.spaceStepsAligned = [trackStart:tracks+trackStart];
    end
    
    numBins = length(param.spaceSteps);
    LickOverDist = zeros(nTrials,numBins);
    lickAlignedOverDist = zeros(nTrials,length(param.spaceStepsAligned));
    
    for tr = 1:nTrials
        if(indTr(end) == totNTrial)
            continue;
        end
            
        indTrCur = indTr(tr);
        indStart = Laps.startLfpInd(indTrCur);
        indEnd = Laps.startLfpInd(indTrCur+1)-1;
        indNextEnd = Laps.endLfpInd(indTrCur+1);
        indLick = behEventsTdt.lick(:,2) >= indStart & ...
            behEventsTdt.lick(:,2) <= indEnd;
        licks = behEventsTdt.lick(indLick,2);
        indLickAll =  behEventsTdt.lick(:,2) >= indStart & ...
            behEventsTdt.lick(:,2) <= indNextEnd;
        licksAll = behEventsTdt.lick(indLickAll,2) - indStart + 1;
               
        if(~isempty(licks))
            licksXMM = Track.xMMAll(licks);
        else
            continue;
        end
        lickTmp = hist(licksXMM,param.spaceSteps);
        lickOverDist(tr,:) = lickTmp;
        
        %% align licks to the pump on
        pump = Laps.pumpLfpInd{indTrCur}(1) - indStart + 1;
        if(~isempty(pump))
            xMMAll = Track.xMMAll(indStart:indNextEnd); 
            xMMAll1 = xMMContinuous(xMMAll);
            xMMLicks = xMMAll1(licksAll) - xMMAll1(pump(1));
            indLicks = xMMLicks >= trackEnd & xMMLicks <= trackStart;
            xMMLicks = xMMLicks(indLicks);
        else
            continue;
        end
              
        lickTmp = hist(xMMLicks,param.spaceStepsAligned);
        lickAlignedOverDist(tr,:) = lickTmp;        
    end
    
    lickOverDistMean = mean(lickOverDist);
    lickOverDistStd = std(lickOverDist);
    lickOverDistSEM = std(lickOverDist)/sqrt(nTrials);
    
    lickAlignedOverDistMean = mean(lickAlignedOverDist);
    lickAlignedOverDistStd = std(lickAlignedOverDist);
    lickAlignedOverDistSEM = std(lickAlignedOverDist)/sqrt(nTrials);
    
    save([path fileName(1:end) '_lickDist.mat'],'lickOverDist',...
        'lickOverDistMean','lickOverDistStd','lickOverDistSEM',...
        'lickAlignedOverDist','lickAlignedOverDistMean',...
        'lickAlignedOverDistStd','lickAlignedOverDistSEM','param');
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
