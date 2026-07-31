function [recDataRunPre,recDataRunManip,recDataRunPost] = ...
                accumRecDataRun(path,recSess,manipSess,mazeSess)
    nRec = size(path,1);
    tmp.numSamplesMean = [];
    tmp.maxSpeedMean = [];
    tmp.meanSpeedMean = [];
    tmp.maxRunLenTMean = [];
    tmp.totRunLenTMean = [];
    tmp.numRunMean = [];
    tmp.maxAccMean = [];
    tmp.meanAccMean = [];
    tmp.totStopLenTMean = [];
    tmp.startCueToRunMean = [];
    tmp.numLicksBefRewMean = [];
    tmp.numLicksRewMean = [];
    tmp.med1stFiveLickDistMean = [];
    tmp.medLickDistMean = [];
    tmp.med1stFiveLickDistBefRewMean = [];
    tmp.medLickDistBefRewMean = [];
    tmp.percRewarded = [];
    tmp.percNonStop = [];
%     tmp.speedProfile = [];
%     tmp.lickProfile = [];
    tmp.pumpLfpIndMean = [];
    tmp.pumpMMMean = [];
    tmp.tDiffInjStart = [];
    tmp.tDiffInjEnd = [];
    tmp.speedProfile = [];
    tmp.lickProfile = [];
    tmp.speedProfileRndSel = [];
    tmp.lickProfileRndSel = [];
    tmp.indRec = [];
    tmp.indSess = [];
    
    recDataRunPre = tmp;
    recDataRunManip = tmp;
    recDataRunPost = tmp;
    
    for i = 1:nRec
        if(i == 18)
            a = 1;
        end
        ind = findstr(path(i,:),'\');
        recName = path(i,ind(end-1)+1:ind(end)-1);
        disp(recName);
        fullPath = [path(i,:) recName '_compSess.mat'];
        if(~exist(fullPath,'file'))
            disp([fullPath ' does not exist.']);
            return;
        end
        load(fullPath,'sessDataRun','sessDataRecTime','sessDataLick','sessDataSpeed');
        
        indPre = find(recSess{i} < manipSess{i}(1));
        if(isempty(indPre))
            indPre = NaN;
        end
        indPost = find(recSess{i} > manipSess{i}(end));
        if(isempty(indPost))
            indPost = NaN;
        end
        [~,indManip] = intersect(recSess{i}, manipSess{i});
        nSessManip = length(manipSess{i});
        recDataRunPre = concatData(recDataRunPre,sessDataRun,sessDataRecTime,...
                    sessDataLick,sessDataSpeed,indPre*ones(1,nSessManip));
        recDataRunPre.indRec = [recDataRunPre.indRec i*ones(1,nSessManip)];
        if(isnan(indPre))
            recDataRunPre.indSess = [recDataRunPre.indSess indPre*ones(1,nSessManip)];
        else
            recDataRunPre.indSess = [recDataRunPre.indSess recSess{i}(indPre)*ones(1,nSessManip)];
        end
        
        recDataRunPost = concatData(recDataRunPost,sessDataRun,sessDataRecTime,...
            sessDataLick,sessDataSpeed,indPost*ones(1,nSessManip));
        recDataRunPost.indRec = [recDataRunPost.indRec i*ones(1,nSessManip)];
        if(isnan(indPost))
            recDataRunPost.indSess = [recDataRunPost.indSess indPost*ones(1,nSessManip)];
        else
            recDataRunPost.indSess = [recDataRunPost.indSess recSess{i}(indPost)*ones(1,nSessManip)];
        end
        
        recDataRunManip = concatData(recDataRunManip,sessDataRun,sessDataRecTime,...
            sessDataLick,sessDataSpeed,indManip);
        recDataRunManip.indRec = [recDataRunManip.indRec i*ones(1,nSessManip)];
        recDataRunManip.indSess = [recDataRunManip.indSess recSess{i}(indManip')];
    end

end

function recData = concatData(recData,sessData,time,sessDataLick,sessDataSpeed,indSess)
    randTrNo = 20;
    nSess = length(indSess);
    if(sum(isnan(indSess)) >= 1)
        recData.numSamplesMean = [recData.numSamplesMean NaN*ones(1,nSess)];
        recData.maxSpeedMean = [recData.maxSpeedMean NaN*ones(1,nSess)];
        recData.meanSpeedMean = [recData.meanSpeedMean NaN*ones(1,nSess)];
        recData.maxRunLenTMean = [recData.maxRunLenTMean NaN*ones(1,nSess)];
        recData.totRunLenTMean = [recData.totRunLenTMean NaN*ones(1,nSess)];
        recData.numRunMean = [recData.numRunMean NaN*ones(1,nSess)];
        recData.maxAccMean = [recData.maxAccMean NaN*ones(1,nSess)];
        recData.meanAccMean = [recData.meanAccMean NaN*ones(1,nSess)];
        recData.totStopLenTMean = [recData.totStopLenTMean NaN*ones(1,nSess)];
        recData.startCueToRunMean = [recData.startCueToRunMean NaN*ones(1,nSess)];
        recData.numLicksBefRewMean = [recData.numLicksBefRewMean NaN*ones(1,nSess)];
        recData.numLicksRewMean = [recData.numLicksRewMean NaN*ones(1,nSess)];
        recData.med1stFiveLickDistMean = [recData.med1stFiveLickDistMean NaN*ones(1,nSess)];
        recData.medLickDistMean = [recData.medLickDistMean NaN*ones(1,nSess)];
        recData.med1stFiveLickDistBefRewMean = [recData.med1stFiveLickDistBefRewMean NaN*ones(1,nSess)];
        recData.medLickDistBefRewMean = [recData.medLickDistBefRewMean NaN*ones(1,nSess)];
        recData.percRewarded = [recData.percRewarded NaN*ones(1,nSess)];
        recData.percNonStop = [recData.percNonStop NaN*ones(1,nSess)];
        recData.pumpLfpIndMean = [recData.pumpLfpIndMean NaN*ones(1,nSess)];
        recData.pumpMMMean = [recData.pumpMMMean NaN*ones(1,nSess)];
        recData.tDiffInjStart = [recData.tDiffInjStart NaN*ones(1,nSess)];
        recData.tDiffInjEnd = [recData.tDiffInjEnd NaN*ones(1,nSess)];
    else
        recData.numSamplesMean = [recData.numSamplesMean sessData.numSamplesMean(indSess)];
        recData.maxSpeedMean = [recData.maxSpeedMean sessData.maxSpeedMean(indSess)];
        recData.meanSpeedMean = [recData.meanSpeedMean sessData.meanSpeedMean(indSess)];
        recData.maxRunLenTMean = [recData.maxRunLenTMean sessData.maxRunLenTMean(indSess)];
        recData.totRunLenTMean = [recData.totRunLenTMean sessData.totRunLenTMean(indSess)];
        recData.numRunMean = [recData.numRunMean sessData.numRunMean(indSess)];
        recData.maxAccMean = [recData.maxAccMean sessData.maxAccMean(indSess)];
        recData.meanAccMean = [recData.meanAccMean sessData.meanAccMean(indSess)];
        recData.totStopLenTMean = [recData.totStopLenTMean sessData.totStopLenTMean(indSess)];
        recData.startCueToRunMean = [recData.startCueToRunMean sessData.startCueToRunMean(indSess)];
        recData.numLicksBefRewMean = [recData.numLicksBefRewMean sessData.numLicksBefRewMean(indSess)];
        recData.numLicksRewMean = [recData.numLicksRewMean sessData.numLicksRewMean(indSess)];
        recData.med1stFiveLickDistMean = [recData.med1stFiveLickDistMean sessData.med1stFiveLickDistMean(indSess)];
        recData.medLickDistMean = [recData.medLickDistMean sessData.medLickDistMean(indSess)];
        recData.med1stFiveLickDistBefRewMean = [recData.med1stFiveLickDistBefRewMean sessData.med1stFiveLickDistBefRewMean(indSess)];
        recData.medLickDistBefRewMean = [recData.medLickDistBefRewMean sessData.medLickDistBefRewMean(indSess)];
        recData.percRewarded = [recData.percRewarded sessData.percRewarded(indSess)];
        recData.percNonStop = [recData.percNonStop sessData.percNonStop(indSess)];
        recData.spaceStepsLick = sessDataLick.spaceSteps;
        recData.spaceStepsSpeed = sessDataSpeed.spaceSteps;
        sess = unique(indSess);
        for i = 1:length(sess)
            recData.speedProfile = [recData.speedProfile; sessDataSpeed.meanRun{sess(i)}];
            recData.lickProfile = [recData.lickProfile; sessDataLick.meanRun{sess(i)}];
            
            indTrR = randperm(size(sessDataSpeed.Run{sess(i)},1));
            if(length(indTrR) > randTrNo)
                recData.speedProfileRndSel = [recData.speedProfileRndSel;...
                    sessDataSpeed.Run{sess(i)}(indTrR(1:randTrNo),:)];
                recData.lickProfileRndSel = [recData.lickProfileRndSel;...
                    sessDataLick.Run{sess(i)}(indTrR(1:randTrNo),:)];
            else
                recData.speedProfileRndSel = [recData.speedProfileRndSel;...
                    sessDataSpeed.Run{sess(i)}];
                recData.lickProfileRndSel = [recData.lickProfileRndSel;...
                    sessDataLick.Run{sess(i)}];
            end
        end
        recData.pumpLfpIndMean = [recData.pumpLfpIndMean sessData.pumpLfpIndMean(indSess)];
        recData.pumpMMMean = [recData.pumpMMMean sessData.pumpMMMean(indSess)];
        recData.tDiffInjStart = [recData.tDiffInjStart time.tDiffInjStart(indSess)];
        recData.tDiffInjEnd = [recData.tDiffInjEnd time.tDiffInjEnd(indSess)];
    end
end