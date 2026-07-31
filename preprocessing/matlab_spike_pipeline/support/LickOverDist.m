function LickOverDist(path, fileName, mazeSess)

    
    %%%%%%%%% initialize constants
    fileNameInfo = [fileName '_Info.mat'];     
    fullPath = [path fileNameInfo];
    if(exist(fullPath) == 0)
        BasicInfo_smTr(path,fileName);
    end
    load(fullPath);
    indTr = beh.indTrCtrl(beh.mazeSess(beh.indTrCtrl) == mazeSess);
    nTrials = length(indTr);
    
    totNTrial = length(beh.mazeSess);
    
    fileName = [fileName '.mat'];
    fullPath = [path fileName];
    if(exist(fullPath) == 0)
        disp('Recording file does not exist.');
        return;
    end
    load(fullPath,'trials');
        
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
        indTrCur = indTr(tr);
        licks = trials{indTrCur}.lickLfpInd;        
        if(~isempty(licks))
            licksXMM = trials{indTrCur}.xMMAll(licks);
        else
            continue;
        end
        lickTmp = hist(licksXMM,param.spaceSteps);
        lickOverDist(tr,:) = lickTmp;
        
        %% align licks to the pump on
        pump = trials{indTrCur}.pumpLfpInd;
        if(~isempty(pump))
            if(indTrCur + 1 <= totNTrial)
                xMMAll = [trials{indTrCur}.xMMAll; trials{indTrCur+1}.xMMAll]; 
                xMMAll1 = xMMContinuous(xMMAll);
                licksAll = [trials{indTrCur}.lickLfpInd; ...
                    trials{indTrCur+1}.lickLfpInd + trials{indTrCur}.Nsamples];
                xMMLicks = xMMAll1(licksAll) - xMMAll1(pump(1));
                indLicks = xMMLicks >= trackEnd & xMMLicks <= trackStart;
                xMMLicks = xMMLicks(indLicks);
            else
                continue;
            end
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
    
    save([path fileName(1:end-4) '_lickDist.mat'],'lickOverDist',...
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